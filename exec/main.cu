/**
 * @file main.cu
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Small driver for the GPU multi-point Dubins solver.
 *
 * For the CPU/GPU comparison harness see `examples/MPMDBenchmark/`.
 */

#include <cuda_runtime.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <configuration.cuh>
#include <dp.cuh>
#include <dubins.cuh>
#include <timeperf.hh>

namespace {

//! Kaya's fourth example (see examples/MPMD/MPMD.hh), whose optimum is known.
const std::vector<Configuration2> kExample = {
		Configuration2 (0.5, 1.2, 5.0 * M_PI / 6.0), Configuration2 (0.0, 0.5, ANGLE::FREE),
		Configuration2 (0.5, 0.5, ANGLE::FREE),			Configuration2 (1.0, 0.5, ANGLE::FREE),
		Configuration2 (1.5, 0.5, ANGLE::FREE),			Configuration2 (2.0, 0.5, ANGLE::FREE),
		Configuration2 (2.0, 0.0, ANGLE::FREE),			Configuration2 (1.5, 0.0, ANGLE::FREE),
		Configuration2 (1.0, 0.0, ANGLE::FREE),			Configuration2 (0.5, 0.0, ANGLE::FREE),
		Configuration2 (0.0, 0.0, ANGLE::FREE),			Configuration2 (0.0, -0.5, 0)};

//! Published length of the example above, for `kmax = 3`.
constexpr LEN_T kExampleLength = 7.46756219733842652175326293218;

//! Prints usage.
void
usage (const char* argv0)
{
	std::cout << "Usage: " << argv0 << " [options]\n"
		<< "  --discr <n>       Angle samples per point (default 90)\n"
		<< "  --nref <n>        Refinement rounds (default 4)\n"
		<< "  --kmax <v>        Maximum curvature (default 3)\n"
		<< "  --precision <p>   fp32 | fp64 (default fp64)\n";
}

}	 // namespace

int
main (int argc, char** argv)
{
	int discr = 90;
	int nref = 4;
	double kmax = 3.0;
	std::string precision_name = "fp64";

	for (int i = 1; i < argc; ++i)
	{
		const std::string a = argv[i];
		const bool has_next	= (i + 1 < argc);
		if (a == "--discr" && has_next) { discr = std::stoi (argv[++i]); }
		else if (a == "--nref" && has_next) { nref = std::stoi (argv[++i]); }
		else if (a == "--kmax" && has_next) { kmax = std::stod (argv[++i]); }
		else if (a == "--precision" && has_next) { precision_name = argv[++i]; }
		else if (a == "-h" || a == "--help") { usage (argv[0]); return 0; }
		else
		{
			std::cerr << "Unknown argument: " << a << "\n";
			usage (argv[0]);
			return 1;
		}
	}

	mpdp::gpu::Options opts;
	opts.discr = discr;
	opts.nref	 = nref;
	if (precision_name == "fp32") { opts.precision = mpdp::gpu::Precision::FP32; }
	else if (precision_name == "fp64") { opts.precision = mpdp::gpu::Precision::FP64; }
	else
	{
		std::cerr << "Error: --precision must be fp32 or fp64\n";
		return 1;
	}

	// Pay for the CUDA context before timing anything.
	if (cudaFree (0) != cudaSuccess)
	{
		std::cerr << "No usable CUDA device\n";
		return 1;
	}
	std::cout << "Device: " << mpdp::gpu::deviceName() << "\n"
						<< "discr=" << discr << " nref=" << nref << " kmax=" << kmax
						<< " precision=" << precision_name << "\n\n";

	std::vector<Configuration2> points = kExample;
	std::vector<bool> fixedAngles (points.size(), false);
	fixedAngles.front() = true;
	fixedAngles.back()	= true;
	std::vector<real_type> params = {kmax};

	try
	{
		TimePerf tp;
		tp.start();
		mpdp::gpu::Result res = mpdp::gpu::solveDP (points, fixedAngles, params, opts);
		const double ms				= tp.getTime();

		// The DP accumulates in the kernels' precision; recompute in double from
		// the returned angles to get the length the path actually measures.
		const double dparams[1] = {kmax};
		double exact						= 0.0;
		for (std::size_t i = 0; i + 1 < points.size(); ++i)
		{
			exact += mpdp::gpu::Dubins<double>::length (
					points[i].x(), points[i].y(), res.angles[i], points[i + 1].x(),
					points[i + 1].y(), res.angles[i + 1], dparams);
		}

		std::cout << std::setprecision (17);
		std::cout << "length (DP)        : " << res.length << "\n"
							<< "length (recomputed): " << exact << "\n"
							<< "reference          : " << kExampleLength << "\n"
							<< "error              : " << (exact - kExampleLength) << "\n\n";
		std::cout << std::setprecision (6);
		std::cout << "angles             : ";
		for (const auto a : res.angles) { std::cout << a << " "; }
		std::cout << "\n\n";
		std::cout << "curves evaluated   : " << res.curves << "\n"
							<< "device time        : " << res.device_ms << " ms\n"
							<< "wall time          : " << ms << " ms" << std::endl;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Solver failed: " << e.what() << std::endl;
		return 1;
	}

	return 0;
}
