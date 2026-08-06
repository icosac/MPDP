/**
 * @file bench_cpu.cc
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @brief Reference driver: the current CPU multi-point Dubins solver (srcCC).
 *
 * Links against `MPDPCC` only. It must never be linked together with the CUDA
 * library: both trees define `Dubins`, `Curve` and `Configuration2` with the
 * same mangled names.
 */

#include <iostream>
#include <string>
#include <vector>

#include <bench_driver.hh>

#include <configuration.hh>
#include <dp.hh>
#include <dubins.hh>

using namespace mpdp;
using namespace mpdp::cpu;

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
		points.emplace_back (p.x[i], p.y[i], th);
	}
	return points;
}

/*!
 * The reference length of a solution: the sum of the point-to-point Dubins
 * lengths for the given angles. Every driver reports this alongside the length
 * its own DP claims, which is what catches broken DP bookkeeping.
 * @param p The instance.
 * @param angles One angle per point.
 * @return The total length, or NaN if the angles are unusable.
 */
double
pathLength (const bench::Problem& p, const std::vector<double>& angles)
{
	if (angles.size() != p.size()) { return std::nan (""); }
	std::vector<real_type> params = {p.kmax};
	double total									= 0.0;
	for (std::size_t i = 0; i + 1 < p.size(); ++i)
	{
		Dubins d (p.x[i], p.y[i], angles[i], p.x[i + 1], p.y[i + 1], angles[i + 1], params);
		total += d.l();
	}
	return total;
}

}	 // namespace

int
main (int argc, char** argv)
{
	bench::DriverOptions opts;
	if (!bench::parseCommonOptions (argc, argv, opts)) { return 1; }

#ifdef _OPENMP
	std::cout << "[cpu] built with OpenMP" << std::endl;
#else
	std::cout << "[cpu] built without OpenMP (single-threaded DP)" << std::endl;
#endif

	auto solve = [] (const bench::Problem& p) {
		std::vector<Configuration2> points = makePoints (p);
		std::vector<bool> fixedAngles (points.size(), false);
		fixedAngles.front() = true;
		fixedAngles.back()	= true;
		std::vector<real_type> params = {p.kmax};

		auto ret = DP().solveDP<Dubins> (
				points, fixedAngles, params, p.discr, p.nref, /*saveAngles=*/true);

		bench::SolveOutput o;
		o.length = ret.first;
		o.angles.assign (ret.second.begin(), ret.second.end());
		return o;
	};

	try
	{
		return bench::runDriver (opts, "cpu", solve, pathLength);
	}
	catch (const std::exception& e)
	{
		std::cerr << "[cpu] fatal: " << e.what() << std::endl;
		return 1;
	}
}
