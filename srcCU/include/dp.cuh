/**
 * @file dp.cuh
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Public interface of the GPU dynamic programming solver.
 *
 * The GPU counterpart of `srcCC`'s `DP::solveDP`. Worth knowing:
 *
 *  - The angle sets are built exactly like the CPU solver builds them
 *    (uniform samples plus the per-pair "guess" angles, and the same
 *    refinement schedule), so the result is directly comparable with the CPU
 *    rather than merely similar. With `nref = 0` the two agree bit for bit.
 *  - Rows of the DP matrix have different lengths - the endpoints are fixed
 *    and the guess angles are not evenly distributed - so the matrix is stored
 *    padded to the longest row with an explicit per-row length. Kernels never
 *    read past a row's length, so padding costs memory but no work.
 *  - Every stage of the backward sweep is one grid of threads over the full
 *    (previous angle x next angle) product followed by a reduction.
 *  - The kernels are templated on the scalar type and on the curve family, so
 *    the same source serves fp32/fp64 and, once `rs.cuh` is filled in,
 *    Reeds-Shepp.
 */

#ifndef MPDP_DP_CUH
#define MPDP_DP_CUH

#include <string>
#include <vector>

#include <configuration.cuh>
#include <curve.cuh>
#include <dubins.cuh>
#include <rs.cuh>
#include <typedefs.hh>

namespace mpdp {
namespace gpu {

//! Floating point precision the kernels run in.
enum class Precision {
	FP64 = 0,	 ///< Double, directly comparable with the CPU solver.
	FP32 = 1	 ///< Single, much faster on consumer cards.
};

//! Knobs for `solveDP`.
struct Options {
	int discr				= 90;		///< Angle samples per point in the first round.
	int nref				= 4;		///< Number of refinement rounds after the first one.
	Precision precision = Precision::FP64;	 ///< Precision of the kernels.
	CurveKind curve = CurveKind::DUBINS;		 ///< Point-to-point curve family.
	int threads_y		= 4;	 ///< Warps per block; the x dimension is always a warp.
	//! Cap on the number of blocks used to split one row, bounding scratch memory.
	int max_split		= 256;
	bool save_angles = true;	///< Write the best angles back into `points`.
};

//! What `solveDP` gives back.
struct Result {
	LEN_T length = 0.0;					///< Length of the best path found.
	std::vector<Angle> angles;	///< One angle per point.
	//! Device time: kernels plus transfers, excluding host-side angle sampling.
	double device_ms = 0.0;
	//! Number of point-to-point curves evaluated, summed over rounds and stages.
	unsigned long long curves = 0;
};

/*!
 * Solves the multi-point problem by dynamic programming on the GPU.
 *
 * @param points The points the path must pass through. The angle of a point
 *        whose entry in `fixedAngles` is true is a constraint; the angle of a
 *        free point is only used as the centre of the first sampling window.
 * @param fixedAngles One flag per point, true where the angle is constrained.
 * @param params Curve parameters; `params[0]` is the maximum curvature.
 * @param opts Solver options.
 * @return The best length, the angles achieving it and timing information.
 * @throws std::runtime_error on a CUDA failure, an inconsistent input, or a
 *         curve family that is not implemented yet.
 */
Result
solveDP (
		std::vector<Configuration2>& points,
		const std::vector<bool>& fixedAngles,
		const std::vector<real_type>& params,
		const Options& opts);

/*!
 * Name of the device the solver will run on, for reporting.
 * @return The device name, or a message explaining why there is none.
 */
std::string
deviceName();

}	 // namespace gpu
}	 // namespace mpdp

#endif	// MPDP_DP_CUH
