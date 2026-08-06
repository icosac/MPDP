/**
 * @file bench_gpu.cu
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @brief Driver for the GPU solver (`mpdp::gpu::solveDP`).
 *
 * Links against `MPDPCU` only - see the note in bench_io.hh about why the CPU
 * and CUDA libraries never share a binary.
 */

#include <cuda_runtime.h>

#include <iostream>
#include <string>
#include <vector>

#include <bench_driver.hh>

#include <configuration.cuh>
#include <dp.cuh>
#include <dubins.cuh>

using namespace mpdp;
using namespace mpdp::gpu;

namespace {

/*!
 * Rebuilds the point list of an instance.
 * @param p The instance.
 * @return The `Configuration2` list, with free intermediate angles.
 */
std::vector<Configuration2>
makePoints (const bench::Problem& p)
{
	std::vector<Configuration2> points;
	points.reserve (p.size());
	for (std::size_t i = 0; i < p.size(); ++i)
	{
		double th = ANGLE::FREE;
		if (i == 0) { th = p.th_init; }
		else if (i + 1 == p.size()) { th = p.th_final; }
		points.push_back (Configuration2 (p.x[i], p.y[i], th));
	}
	return points;
}

/*!
 * Independent length check, always in double precision and always with the
 * Dubins solver, so an fp32 run is scored on the true length of the path it
 * chose rather than on its own rounded arithmetic.
 * @param p The instance.
 * @param angles One angle per point.
 * @return The total length, or NaN if the angles are unusable.
 */
double
pathLength (const bench::Problem& p, const std::vector<double>& angles)
{
	if (angles.size() != p.size()) { return std::nan (""); }
	const double params[1] = {p.kmax};
	double total					 = 0.0;
	for (std::size_t i = 0; i + 1 < p.size(); ++i)
	{
		total += mpdp::gpu::Dubins<double>::length (
				p.x[i], p.y[i], angles[i], p.x[i + 1], p.y[i + 1], angles[i + 1], params);
	}
	return total;
}

}	 // namespace

int
main (int argc, char** argv)
{
	bench::DriverOptions opts;
	if (!bench::parseCommonOptions (argc, argv, opts)) { return 1; }

	mpdp::gpu::Precision precision = mpdp::gpu::Precision::FP64;
	std::string precision_name		 = "fp64";
	int threads_y									 = 4;
	for (int i = 1; i < argc; ++i)
	{
		const std::string a = argv[i];
		if (a == "--precision" && i + 1 < argc)
		{
			precision_name = argv[++i];
			if (precision_name == "fp32") { precision = mpdp::gpu::Precision::FP32; }
			else if (precision_name == "fp64") { precision = mpdp::gpu::Precision::FP64; }
			else
			{
				std::cerr << "Error: --precision must be fp32 or fp64\n";
				return 1;
			}
		}
		else if (a == "--threads-y" && i + 1 < argc) { threads_y = std::stoi (argv[++i]); }
	}

	// Pay for context creation before anything is timed.
	if (cudaFree (0) != cudaSuccess)
	{
		std::cerr << "[gpu] no usable CUDA device\n";
		return 1;
	}
	std::cout << "[gpu] device: " << mpdp::gpu::deviceName()
						<< ", precision: " << precision_name << std::endl;

	auto solve = [&] (const bench::Problem& p) {
		std::vector<Configuration2> points = makePoints (p);
		std::vector<bool> fixedAngles (points.size(), false);
		fixedAngles.front() = true;
		fixedAngles.back()	= true;
		std::vector<real_type> params = {p.kmax};

		mpdp::gpu::Options o;
		o.discr			= p.discr;
		o.nref			= p.nref;
		o.precision = precision;
		o.curve			= mpdp::gpu::CurveKind::DUBINS;
		o.threads_y = threads_y;

		mpdp::gpu::Result res = mpdp::gpu::solveDP (points, fixedAngles, params, o);

		bench::SolveOutput out;
		out.length				 = res.length;
		out.angles				 = res.angles;
		out.time_solver_ms = res.device_ms;
		return out;
	};

	try
	{
		return bench::runDriver (opts, "gpu_" + precision_name, solve, pathLength);
	}
	catch (const std::exception& e)
	{
		std::cerr << "[gpu] fatal: " << e.what() << std::endl;
		return 1;
	}
}
