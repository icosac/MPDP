/**
 * @file dp.cu
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief GPU dynamic programming solver for multi-point paths.
 *
 * Structure of the algorithm, and where the parallelism is:
 *
 *   for each round (one coarse round + `nref` refinements)      <- sequential
 *       build the angle set of every point                      <- host, O(N*discr)
 *       for idx = N-1 .. 1                                      <- sequential
 *           for every (angle of point idx-1, angle of point idx) <- PARALLEL
 *               evaluate the point-to-point curve and add the
 *               already-known cost of the tail
 *           reduce over the second index                        <- PARALLEL
 *       walk the `next` pointers back to the angles             <- one thread
 *
 */

#ifndef CUDA_ON
#error "dp.cu must be compiled with CUDA_ON defined"
#endif

#include <dp.cuh>

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mpdp {
namespace gpu {

namespace {

	/*!
	 * Throws on a CUDA error, with the call site attached.
	 * @param err The status to check.
	 * @param what What was being attempted.
	 */
	inline void
	check (cudaError_t err, const char* what)
	{
		if (err != cudaSuccess)
		{
			std::ostringstream oss;
			oss << "CUDA error in " << what << ": " << cudaGetErrorString (err);
			throw std::runtime_error (oss.str());
		}
	}

//! Checks the status of the most recent kernel launch.
#define MPDP_CHECK_LAUNCH(what)                 \
	do {                                          \
		check (cudaGetLastError(), (what));         \
	} while (0)

	//////////////////////////////////////////////////////////////////////////////
	// Angle sampling - a faithful port of srcCC/dp.cc
	//////////////////////////////////////////////////////////////////////////////

	/*!
	 * Returns up to two circles of radius `r` through two points.
	 * Credit to Marco Frego & Paolo Bevilacqua.
	 * @param x1 Abscissa of the first point.
	 * @param y1 Ordinate of the first point.
	 * @param x2 Abscissa of the second point.
	 * @param y2 Ordinate of the second point.
	 * @param r The radius.
	 * @param XC Abscissas of the found centres.
	 * @param YC Ordinates of the found centres.
	 */
	void
	circles (
			double x1,
			double y1,
			double x2,
			double y2,
			double r,
			std::vector<double>& XC,
			std::vector<double>& YC)
	{
		const double TOL = 1e-8;

		const double q	= std::hypot (x2 - x1, y2 - y1);
		const double x3 = 0.5 * (x1 + x2);
		const double y3 = 0.5 * (y1 + y2);

		const double delta = r * r - q * q / 4.;

		XC.clear();
		YC.clear();

		if (delta < -TOL) { return; }

		if (delta < TOL)
		{
			XC.push_back (x3);
			YC.push_back (y3);
		}
		else
		{
			const double deltaS = std::sqrt (delta);
			XC.push_back (x3 + deltaS * (y1 - y2) / q);
			YC.push_back (y3 + deltaS * (x2 - x1) / q);
			XC.push_back (x3 - deltaS * (y1 - y2) / q);
			YC.push_back (y3 - deltaS * (x2 - x1) / q);
		}
	}

	/*!
	 * The "interesting" angles between two consecutive points: the heading that
	 * lines them up, and the tangents to the two circles of radius 1/Kmax
	 * through both. Credit to Marco Frego & Paolo Bevilacqua.
	 *
	 * @param i Index of the second of the two points.
	 * @param thPrev Filled with the angles suggested for point i-1.
	 * @param thCur Filled with the angles suggested for point i.
	 * @param points All the points.
	 * @param kmax The maximum curvature.
	 */
	void
	guessInitialAngles (
			const std::size_t i,
			std::vector<double>& thPrev,
			std::vector<double>& thCur,
			const std::vector<Configuration2>& points,
			double kmax)
	{
		thPrev.clear();
		thCur.clear();
		thPrev.reserve (5);
		thCur.reserve (5);

		const double m_pi = 3.14159265358979323846264338328;

		// aligned on straight line
		double th = std::atan2 (
				points[i].y() - points[i - 1].y(), points[i].x() - points[i - 1].x());
		thPrev.push_back (th);
		thCur.push_back (th);

		// aligned on circle
		std::vector<double> XC, YC;
		circles (
				points[i - 1].x(), points[i - 1].y(), points[i].x(), points[i].y(), 1. / kmax,
				XC, YC);
		for (std::size_t j = 0; j < XC.size(); ++j)
		{
			th = std::atan2 (points[i - 1].y() - YC[j], points[i - 1].x() - XC[j]);
			thPrev.push_back (th + m_pi / 2.);
			thPrev.push_back (th - m_pi / 2.);
			th = std::atan2 (points[i].y() - YC[j], points[i].x() - XC[j]);
			thCur.push_back (th + m_pi / 2.);
			thCur.push_back (th - m_pi / 2.);
		}
	}

	/*!
	 * Builds the angle set of the first round: `discr` samples spread over the
	 * whole circle around each point's current angle, plus the guess angles.
	 *
	 * @param discr Number of uniform samples.
	 * @param fixedAngles One flag per point.
	 * @param points The points.
	 * @param kmax The maximum curvature.
	 * @return One vector of angles per point.
	 */
	std::vector<std::vector<double>>
	samplingAnglesFirst (
			int discr,
			const std::vector<bool>& fixedAngles,
			const std::vector<Configuration2>& points,
			double kmax)
	{
		const std::size_t n = points.size();
		std::vector<std::vector<double>> rows (n);

		const double m_pi		 = 3.14159265358979323846264338328;
		const double dtheta	 = 2 * m_pi / discr;

		for (std::size_t i = 0; i < n; ++i)
		{
			rows[i].reserve (discr + 12);
			for (int j = 0; j < discr; ++j)
			{
				rows[i].push_back (points[i].th() + dtheta * j);
			}
			if (i > 0)
			{
				std::vector<double> thPrev, thCur;
				guessInitialAngles (i, thPrev, thCur, points, kmax);
				rows[i - 1].insert (rows[i - 1].end(), thPrev.begin(), thPrev.end());
				rows[i].insert (rows[i].end(), thCur.begin(), thCur.end());
			}
		}
		for (std::size_t i = 0; i < n; ++i)
		{
			if (fixedAngles[i]) { rows[i] = {points[i].th()}; }
		}
		return rows;
	}

	/*!
	 * Builds the angle set of a refinement round: a window of half-width
	 * `hrange` around each point's current angle, plus the guess angles.
	 *
	 * @param points The points, carrying the angles found by the previous round.
	 * @param fixedAngles One flag per point.
	 * @param hrange Half-width of the sampling window.
	 * @param hn Number of samples on each side of the centre.
	 * @param kmax The maximum curvature.
	 * @return One vector of angles per point.
	 */
	std::vector<std::vector<double>>
	samplingAnglesRefine (
			const std::vector<Configuration2>& points,
			const std::vector<bool>& fixedAngles,
			double hrange,
			int hn,
			double kmax)
	{
		const std::size_t n = points.size();
		std::vector<std::vector<double>> rows (n);

		const double dtheta = hrange / hn;
		for (std::size_t i = 0; i < n; ++i)
		{
			rows[i].reserve (2 * hn + 12);
			rows[i].push_back (points[i].th());
			for (int j = 1; j <= hn; ++j)
			{
				rows[i].push_back (dtheta * j + points[i].th());
				rows[i].push_back (-dtheta * j + points[i].th());
			}
			if (i >= 1)
			{
				std::vector<double> thPrev, thCur;
				guessInitialAngles (i, thPrev, thCur, points, kmax);
				rows[i - 1].insert (rows[i - 1].end(), thPrev.begin(), thPrev.end());
				rows[i].insert (rows[i].end(), thCur.begin(), thCur.end());
			}
		}
		for (std::size_t i = 0; i < n; ++i)
		{
			if (fixedAngles[i]) { rows[i] = {points[i].th()}; }
		}
		return rows;
	}

	//////////////////////////////////////////////////////////////////////////////
	// Kernels
	//////////////////////////////////////////////////////////////////////////////

	/*!
	 * Initialises the cost matrix: the last row costs nothing (it is where the
	 * path ends), every other cell starts unreachable.
	 * @param len The cost matrix, row-major, padded to `W`.
	 * @param next The successor matrix.
	 * @param row_len Number of valid cells per row.
	 * @param W Padded row width.
	 * @param n Number of rows.
	 */
	template <typename T>
	__global__ void
	dpInitKernel (T* len, int* next, const int* row_len, int W, int n)
	{
		const int tid		 = blockIdx.x * blockDim.x + threadIdx.x;
		const int stride = gridDim.x * blockDim.x;
		for (int c = tid; c < n * W; c += stride)
		{
			const int row = c / W;
			const int col = c - row * W;
			const bool valid_last = (row == n - 1) && (col < row_len[row]);
			len[c]	= valid_last ? T (0) : Scalar<T>::huge();
			next[c] = -1;
		}
	}

	/*!
	 * One stage of the backward sweep.
	 *
	 * Thread (x, y) of block (bx, by) handles previous-angle
	 * `i = bx*blockDim.y + y` and strides over next-angles starting at
	 * `by*32 + x`. The reduction over the next-angle index is done inside the
	 * warp, which is why `blockDim.x` is fixed to the warp size: every warp owns
	 * exactly one `i`, so no shared memory or barrier is needed.
	 *
	 * When the grid splits a row over several blocks in y, each block writes a
	 * partial best into `dst_*` and `dpStageCombineKernel` finishes the job;
	 * when it does not, `dst_*` are the destination row itself and there is
	 * nothing left to combine.
	 *
	 * @param x0 Abscissa of the previous point.
	 * @param y0 Ordinate of the previous point.
	 * @param x1 Abscissa of the current point.
	 * @param y1 Ordinate of the current point.
	 * @param th_prev Angles of the previous point.
	 * @param th_cur Angles of the current point.
	 * @param len_cur Already computed tail cost for each angle of the current point.
	 * @param params Curve parameters.
	 * @param dst_len Where to write the best cost.
	 * @param dst_idx Where to write the index achieving it.
	 * @param i_count Number of angles of the previous point.
	 * @param j_count Number of angles of the current point.
	 * @param dst_stride Number of partials stored per `i`.
	 */
	template <typename T, class CurveT>
	__global__ void
	dpStageKernel (
			T x0,
			T y0,
			T x1,
			T y1,
			const T* __restrict__ th_prev,
			const T* __restrict__ th_cur,
			const T* __restrict__ len_cur,
			const T* __restrict__ params,
			T* __restrict__ dst_len,
			int* __restrict__ dst_idx,
			int i_count,
			int j_count,
			int dst_stride)
	{
		const int i = blockIdx.x * blockDim.y + threadIdx.y;
		// `i` does not depend on threadIdx.x, so a whole warp leaves together and
		// the full-mask shuffles below stay legal.
		if (i >= i_count) { return; }

		const T th0 = th_prev[i];

		T best_l	 = Scalar<T>::huge();
		int best_j = -1;

		const int j_start	 = blockIdx.y * blockDim.x + threadIdx.x;
		const int j_stride = gridDim.y * blockDim.x;
		for (int j = j_start; j < j_count; j += j_stride)
		{
			const T l =
					CurveT::length (x0, y0, th0, x1, y1, th_cur[j], params) + len_cur[j];
			if (l < best_l || (l == best_l && j < best_j))
			{
				best_l = l;
				best_j = j;
			}
		}

		#pragma unroll
		for (int off = 16; off > 0; off >>= 1)
		{
			const T other_l	 = __shfl_down_sync (0xffffffffu, best_l, off);
			const int other_j = __shfl_down_sync (0xffffffffu, best_j, off);
			if (other_l < best_l)
			{
				best_l = other_l;
				best_j = other_j;
			}
			else if (other_l == best_l && other_j >= 0 && (best_j < 0 || other_j < best_j))
			{
				best_j = other_j;
			}
		}

		if (threadIdx.x == 0)
		{
			dst_len[i * dst_stride + blockIdx.y] = best_l;
			dst_idx[i * dst_stride + blockIdx.y] = best_j;
		}
	}

	/*!
	 * Merges the partial minima produced when a row was split over several
	 * blocks.
	 * @param part_len The partial costs.
	 * @param part_idx The partial indices.
	 * @param len_prev Destination cost row.
	 * @param next_prev Destination successor row.
	 * @param i_count Number of angles of the previous point.
	 * @param stride Number of partials stored per `i`.
	 */
	template <typename T>
	__global__ void
	dpStageCombineKernel (
			const T* __restrict__ part_len,
			const int* __restrict__ part_idx,
			T* __restrict__ len_prev,
			int* __restrict__ next_prev,
			int i_count,
			int stride)
	{
		const int i = blockIdx.x * blockDim.x + threadIdx.x;
		if (i >= i_count) { return; }

		T best_l	 = Scalar<T>::huge();
		int best_j = -1;
		for (int b = 0; b < stride; ++b)
		{
			const int j = part_idx[i * stride + b];
			if (j < 0) { continue; }
			const T l = part_len[i * stride + b];
			if (l < best_l || (l == best_l && j < best_j))
			{
				best_l = l;
				best_j = j;
			}
		}
		len_prev[i]	 = best_l;
		next_prev[i] = best_j;
	}

	/*!
	 * Picks the best cell of the first row and walks the successor pointers to
	 * recover the angles. Sequential by nature and only O(N), so it runs in a
	 * single thread rather than paying a round trip to the host.
	 * @param theta The angle matrix.
	 * @param len The cost matrix.
	 * @param next The successor matrix.
	 * @param row_len Number of valid cells per row.
	 * @param W Padded row width.
	 * @param n Number of rows.
	 * @param out_angles Receives one angle per point.
	 * @param out_len Receives the best cost, or a negative value if the path is
	 *        broken.
	 */
	template <typename T>
	__global__ void
	dpTraceKernel (
			const T* __restrict__ theta,
			const T* __restrict__ len,
			const int* __restrict__ next,
			const int* __restrict__ row_len,
			int W,
			int n,
			T* __restrict__ out_angles,
			T* __restrict__ out_len)
	{
		if (threadIdx.x != 0 || blockIdx.x != 0) { return; }

		T best_l	 = Scalar<T>::huge();
		int best_j = -1;
		for (int j = 0; j < row_len[0]; ++j)
		{
			if (len[j] < best_l)
			{
				best_l = len[j];
				best_j = j;
			}
		}

		if (best_j < 0)
		{
			out_len[0] = T (-1);
			return;
		}

		out_len[0]		= best_l;
		out_angles[0] = theta[best_j];

		int cur = best_j;
		for (int r = 0; r + 1 < n; ++r)
		{
			const int nx = next[r * W + cur];
			if (nx < 0 || nx >= row_len[r + 1])
			{
				out_len[0] = T (-2);
				return;
			}
			out_angles[r + 1] = theta[(r + 1) * W + nx];
			cur								= nx;
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	// Device memory held for the duration of one solve
	//////////////////////////////////////////////////////////////////////////////

	//! Owns every device allocation of a single `solveDP` call.
	template <typename T>
	struct DeviceWorkspace {
		T* theta			= nullptr;
		T* len				= nullptr;
		int* next			= nullptr;
		int* row_len	= nullptr;
		T* params			= nullptr;
		T* part_len		= nullptr;
		int* part_idx = nullptr;
		T* out_angles = nullptr;
		T* out_len		= nullptr;

		int n = 0;	///< Rows.
		int W = 0;	///< Padded row width.
		int max_split = 0;	///< Partials stored per `i`.

		//! Frees everything; safe to call twice.
		~DeviceWorkspace()
		{
			cudaFree (theta);
			cudaFree (len);
			cudaFree (next);
			cudaFree (row_len);
			cudaFree (params);
			cudaFree (part_len);
			cudaFree (part_idx);
			cudaFree (out_angles);
			cudaFree (out_len);
		}

		/*!
		 * Allocates for a matrix of `rows` x `width` and the given split factor.
		 * @param rows Number of points.
		 * @param width Padded row width.
		 * @param split Maximum number of partials per `i`.
		 * @param n_params Number of curve parameters.
		 */
		void
		allocate (int rows, int width, int split, int n_params)
		{
			n					= rows;
			W					= width;
			max_split = split;

			const std::size_t cells = (std::size_t)rows * (std::size_t)width;
			check (cudaMalloc (&theta, cells * sizeof (T)), "cudaMalloc(theta)");
			check (cudaMalloc (&len, cells * sizeof (T)), "cudaMalloc(len)");
			check (cudaMalloc (&next, cells * sizeof (int)), "cudaMalloc(next)");
			check (cudaMalloc (&row_len, (std::size_t)rows * sizeof (int)), "cudaMalloc(row_len)");
			check (cudaMalloc (&params, (std::size_t)n_params * sizeof (T)), "cudaMalloc(params)");

			const std::size_t parts = (std::size_t)width * (std::size_t)split;
			check (cudaMalloc (&part_len, parts * sizeof (T)), "cudaMalloc(part_len)");
			check (cudaMalloc (&part_idx, parts * sizeof (int)), "cudaMalloc(part_idx)");

			check (cudaMalloc (&out_angles, (std::size_t)rows * sizeof (T)), "cudaMalloc(out_angles)");
			check (cudaMalloc (&out_len, sizeof (T)), "cudaMalloc(out_len)");
		}
	};

	//! A pair of CUDA events, destroyed even when a solve throws.
	struct EventPair {
		cudaEvent_t start = nullptr;
		cudaEvent_t stop	= nullptr;

		EventPair()
		{
			check (cudaEventCreate (&start), "cudaEventCreate(start)");
			check (cudaEventCreate (&stop), "cudaEventCreate(stop)");
		}
		~EventPair()
		{
			cudaEventDestroy (start);
			cudaEventDestroy (stop);
		}
		EventPair (const EventPair&)						 = delete;
		EventPair& operator= (const EventPair&) = delete;
	};

	//////////////////////////////////////////////////////////////////////////////
	// Host driver
	//////////////////////////////////////////////////////////////////////////////

	/*!
	 * The solver, specialised on precision and curve family.
	 * @param points The points; angles are updated if `opts.save_angles`.
	 * @param fixedAngles One flag per point.
	 * @param params Curve parameters.
	 * @param opts Solver options.
	 * @return The best length, the angles, and timing information.
	 */
	template <typename T, class CurveT>
	Result
	solveDPImpl (
			std::vector<Configuration2>& points,
			const std::vector<bool>& fixedAngles,
			const std::vector<real_type>& params,
			const Options& opts)
	{
		const int n = (int)points.size();

		std::vector<Configuration2> comp_points;
		comp_points.reserve (n);
		for (const auto& p : points)
		{
			comp_points.push_back (Configuration2 (p.x(), p.y(), p.th()));
		}

		const double kmax = (double)params[0];
		
		const int hn			= opts.discr / 2;
		const int w_first = opts.discr + 10;
		const int w_ref		= 2 * hn + 1 + 10;
		const int W				= std::max (w_first, w_ref);
		const int split		= std::max (1, std::min (opts.max_split, (W + 31) / 32));

		DeviceWorkspace<T> ws;
		ws.allocate (n, W, split, CurveT::kNumParams);

		{
			std::vector<T> h_params (CurveT::kNumParams);
			for (int i = 0; i < CurveT::kNumParams; ++i) { h_params[i] = (T)params[i]; }
			check (
					cudaMemcpy (
							ws.params, h_params.data(), h_params.size() * sizeof (T),
							cudaMemcpyHostToDevice),
					"cudaMemcpy(params)");
		}

		std::vector<T> h_theta ((std::size_t)n * W);
		std::vector<int> h_row_len (n);
		std::vector<T> h_angles (n);
		T h_len = T (0);

		EventPair ev;

		Result result;
		double hrange = 2.0 * 3.14159265358979323846264338328;

		const int warps_y = std::max (1, opts.threads_y);

		for (int round = 0; round <= opts.nref; ++round)
		{
			std::vector<std::vector<double>> rows;
			if (round == 0)
			{
				rows = samplingAnglesFirst (opts.discr, fixedAngles, comp_points, kmax);
			}
			else
			{
				// Same schedule as the CPU: the window shrinks by discr/1.5 each time.
				hrange = hrange / opts.discr * 1.5;
				rows	 = samplingAnglesRefine (comp_points, fixedAngles, hrange, hn, kmax);
			}

			int w_used = 0;
			for (int i = 0; i < n; ++i)
			{
				h_row_len[i] = (int)rows[i].size();
				w_used			 = std::max (w_used, h_row_len[i]);
				for (int j = 0; j < h_row_len[i]; ++j)
				{
					h_theta[(std::size_t)i * W + j] = (T)rows[i][j];
				}
			}
			if (w_used > W)
			{
				throw std::runtime_error (
						"mpdp::gpu::solveDP: row wider than the reserved padding");
			}

			check (cudaEventRecord (ev.start), "cudaEventRecord");

			check (
					cudaMemcpy (
							ws.theta, h_theta.data(), (std::size_t)n * W * sizeof (T),
							cudaMemcpyHostToDevice),
					"cudaMemcpy(theta)");
			check (
					cudaMemcpy (
							ws.row_len, h_row_len.data(), (std::size_t)n * sizeof (int),
							cudaMemcpyHostToDevice),
					"cudaMemcpy(row_len)");

			{
				const int threads = 256;
				const int blocks	= std::min (4096, (n * W + threads - 1) / threads);
				dpInitKernel<T><<<blocks, threads>>> (ws.len, ws.next, ws.row_len, W, n);
				MPDP_CHECK_LAUNCH ("dpInitKernel");
			}

			for (int idx = n - 1; idx >= 1; --idx)
			{
				const int i_count = h_row_len[idx - 1];
				const int j_count = h_row_len[idx];
				result.curves += (unsigned long long)i_count * (unsigned long long)j_count;

				const dim3 block (32, warps_y);
				const int gx = (i_count + warps_y - 1) / warps_y;
				const int gy =
						std::max (1, std::min (opts.max_split, (j_count + 31) / 32));
				const dim3 grid (gx, gy);

				// With a single block in y the kernel can write straight into the
				// destination row; the combine pass exists only for wider splits.
				T* dst_len		= (gy == 1) ? ws.len + (std::size_t)(idx - 1) * W : ws.part_len;
				int* dst_idx	= (gy == 1) ? ws.next + (std::size_t)(idx - 1) * W : ws.part_idx;
				const int dst_stride = (gy == 1) ? 1 : gy;

				dpStageKernel<T, CurveT><<<grid, block>>> (
						(T)comp_points[idx - 1].x(), (T)comp_points[idx - 1].y(),
						(T)comp_points[idx].x(), (T)comp_points[idx].y(),
						ws.theta + (std::size_t)(idx - 1) * W, ws.theta + (std::size_t)idx * W,
						ws.len + (std::size_t)idx * W, ws.params, dst_len, dst_idx, i_count,
						j_count, dst_stride);
				MPDP_CHECK_LAUNCH ("dpStageKernel");

				if (gy > 1)
				{
					const int threads = 128;
					const int blocks	= (i_count + threads - 1) / threads;
					dpStageCombineKernel<T><<<blocks, threads>>> (
							ws.part_len, ws.part_idx, ws.len + (std::size_t)(idx - 1) * W,
							ws.next + (std::size_t)(idx - 1) * W, i_count, gy);
					MPDP_CHECK_LAUNCH ("dpStageCombineKernel");
				}
			}

			dpTraceKernel<T><<<1, 1>>> (
					ws.theta, ws.len, ws.next, ws.row_len, W, n, ws.out_angles, ws.out_len);
			MPDP_CHECK_LAUNCH ("dpTraceKernel");

			check (
					cudaMemcpy (
							h_angles.data(), ws.out_angles, (std::size_t)n * sizeof (T),
							cudaMemcpyDeviceToHost),
					"cudaMemcpy(out_angles)");
			check (
					cudaMemcpy (&h_len, ws.out_len, sizeof (T), cudaMemcpyDeviceToHost),
					"cudaMemcpy(out_len)");

			check (cudaEventRecord (ev.stop), "cudaEventRecord");
			check (cudaEventSynchronize (ev.stop), "cudaEventSynchronize");
			float ms = 0.f;
			check (cudaEventElapsedTime (&ms, ev.start, ev.stop), "cudaEventElapsedTime");
			result.device_ms += (double)ms;

			if (h_len < T (0))
			{
				throw std::runtime_error (
						"mpdp::gpu::solveDP: dynamic programming produced no feasible path");
			}

			for (int i = 0; i < n; ++i) { comp_points[i].th ((Angle)h_angles[i]); }
			result.length = (LEN_T)h_len;
		}

		result.angles.resize (n);
		for (int i = 0; i < n; ++i) { result.angles[i] = comp_points[i].th(); }

		if (opts.save_angles)
		{
			for (int i = 0; i < n; ++i) { points[i].th (result.angles[i]); }
		}

		return result;
	}

	/*!
	 * Chooses the curve family for a given precision.
	 * @param points The points.
	 * @param fixedAngles One flag per point.
	 * @param params Curve parameters.
	 * @param opts Solver options.
	 * @return The result of the specialised solver.
	 */
	template <typename T>
	Result
	dispatchCurve (
			std::vector<Configuration2>& points,
			const std::vector<bool>& fixedAngles,
			const std::vector<real_type>& params,
			const Options& opts)
	{
		switch (opts.curve)
		{
			case CurveKind::DUBINS:
				return solveDPImpl<T, Dubins<T>> (points, fixedAngles, params, opts);
			case CurveKind::REEDS_SHEPP:
				if constexpr (ReedsShepp<T>::kImplemented)
				{
					return solveDPImpl<T, ReedsShepp<T>> (points, fixedAngles, params, opts);
				}
				else
				{
					throw std::runtime_error (
							"mpdp::gpu::solveDP: the Reeds-Shepp curve is not implemented yet");
				}
			default:
				throw std::runtime_error ("mpdp::gpu::solveDP: unknown curve kind");
		}
	}

}	 // namespace

Result
solveDP (
		std::vector<Configuration2>& points,
		const std::vector<bool>& fixedAngles,
		const std::vector<real_type>& params,
		const Options& opts)
{
	if (points.size() < 2)
	{
		throw std::runtime_error ("mpdp::gpu::solveDP: at least two points are needed");
	}
	if (points.size() != fixedAngles.size())
	{
		throw std::runtime_error (
				"mpdp::gpu::solveDP: points and fixedAngles have different sizes");
	}
	if (params.empty())
	{
		throw std::runtime_error ("mpdp::gpu::solveDP: params[0] must hold the curvature");
	}
	if (opts.discr < 1)
	{
		throw std::runtime_error ("mpdp::gpu::solveDP: discr must be positive");
	}
	if (opts.nref < 0)
	{
		throw std::runtime_error ("mpdp::gpu::solveDP: nref must not be negative");
	}

	switch (opts.precision)
	{
		case Precision::FP32:
			return dispatchCurve<float> (points, fixedAngles, params, opts);
		case Precision::FP64:
		default:
			return dispatchCurve<double> (points, fixedAngles, params, opts);
	}
}

std::string
deviceName()
{
	int count = 0;
	if (cudaGetDeviceCount (&count) != cudaSuccess || count == 0)
	{
		return "<no CUDA device>";
	}
	int dev = 0;
	cudaGetDevice (&dev);
	cudaDeviceProp prop;
	if (cudaGetDeviceProperties (&prop, dev) != cudaSuccess) { return "<unknown device>"; }
	std::ostringstream oss;
	oss << prop.name << " (sm_" << prop.major << prop.minor << ", " << prop.multiProcessorCount
			<< " SMs)";
	return oss.str();
}

}	 // namespace gpu
}	 // namespace mpdp
