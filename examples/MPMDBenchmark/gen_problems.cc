/**
 * @file gen_problems.cc
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @brief Generates the shared problem set consumed by every solver driver.
 *
 * The whole point of a separate generator is that the CPU, legacy-GPU and
 * new-GPU drivers all read the *same* file, so any difference in the reported
 * lengths is a difference between the solvers and never between their inputs.
 *
 * Two profiles are produced:
 *  - `budget`: the discretisation of every configuration is derived from a
 *    total curve budget (~1e6 by default), so the whole sweep costs about that
 *    many point-to-point Dubins evaluations.
 *  - `stress`: realistic discretisations (90..360) at growing N. This is where
 *    a GPU actually has enough work per kernel to pay off. Bounded by the
 *    per-solver wall-clock cap rather than by the curve budget.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <bench_io.hh>

namespace {

//! Upper bound for the discretisation the `budget` profile may pick.
constexpr int kMaxBudgetDiscr = 512;

//! One point of the (N, discr, nref) sweep.
struct ConfigSpec {
	const char* tag;
	int n_points;
	int discr;	///< Used as-is by `stress`, recomputed from the budget otherwise.
	int nref;
	double kmax;
};

//! Realistic discretisations, growing N. Used by the `stress` profile.
const std::vector<ConfigSpec> kStressSweep = {
		{"n3", 3, 360, 4, 3.0},			{"n5", 5, 360, 4, 3.0},
		{"n10", 10, 360, 4, 3.0},		{"n25", 25, 180, 4, 3.0},
		{"n50", 50, 180, 4, 3.0},		{"n100", 100, 90, 4, 3.0},
		{"n200", 200, 90, 4, 3.0},
};

//! Same N ladder; `discr` is overwritten from the curve budget.
const std::vector<ConfigSpec> kBudgetSweep = {
		{"n3", 3, 0, 4, 3.0},				{"n5", 5, 0, 4, 3.0},
		{"n10", 10, 0, 4, 3.0},			{"n25", 25, 0, 4, 3.0},
		{"n50", 50, 0, 4, 3.0},			{"n100", 100, 0, 4, 3.0},
		{"n200", 200, 0, 4, 3.0},
};

/*!
 * Draws one random instance.
 *
 * The bounding box grows with sqrt(N) so that the point density - and with it
 * the difficulty relative to the turning radius - stays comparable across the
 * whole N ladder.
 *
 * @param id The instance id.
 * @param cfg The configuration the instance belongs to.
 * @param gen The random engine.
 * @return The generated instance.
 */
bench::Problem
makeProblem (int id, const ConfigSpec& cfg, std::mt19937& gen)
{
	const double side = 10.0 * std::sqrt (cfg.n_points / 10.0 > 1.0 ? cfg.n_points / 10.0 : 1.0);
	std::uniform_real_distribution<double> coord (0.0, side);
	std::uniform_real_distribution<double> angle (0.0, 2.0 * M_PI);

	bench::Problem p;
	p.id			 = id;
	p.tag			 = cfg.tag;
	p.kmax		 = cfg.kmax;
	p.discr		 = cfg.discr;
	p.nref		 = cfg.nref;
	p.th_init	 = angle (gen);
	p.th_final = angle (gen);
	for (int i = 0; i < cfg.n_points; ++i)
	{
		p.x.push_back (coord (gen));
		p.y.push_back (coord (gen));
	}
	return p;
}

//! Prints usage.
void
usage (const char* argv0)
{
	std::cout
			<< "Usage: " << argv0 << " [options]\n"
			<< "  --out <path>       Output problem file (required)\n"
			<< "  --profile <name>   budget | stress          (default: budget)\n"
			<< "  --budget <n>       Target total curves for the budget profile"
				 " (default: 1000000)\n"
			<< "  --instances <n>    Random instances per configuration (default: 1)\n"
			<< "  --seed <n>         RNG seed (default: 13)\n";
}

}	 // namespace

int
main (int argc, char** argv)
{
	std::string out_path;
	std::string profile = "budget";
	std::uint64_t budget = 1000000;
	int instances				 = 1;
	unsigned seed				 = 13;

	for (int i = 1; i < argc; ++i)
	{
		const std::string a = argv[i];
		auto next						= [&] (const char* what) -> std::string {
			 if (i + 1 >= argc)
			 {
				 throw std::runtime_error (std::string ("missing value for ") + what);
			 }
			 return argv[++i];
		};
		if (a == "--out") { out_path = next ("--out"); }
		else if (a == "--profile") { profile = next ("--profile"); }
		else if (a == "--budget") { budget = std::stoull (next ("--budget")); }
		else if (a == "--instances") { instances = std::stoi (next ("--instances")); }
		else if (a == "--seed") { seed = (unsigned)std::stoul (next ("--seed")); }
		else if (a == "-h" || a == "--help") { usage (argv[0]); return 0; }
		else
		{
			std::cerr << "Unknown argument: " << a << "\n";
			usage (argv[0]);
			return 1;
		}
	}

	if (out_path.empty())
	{
		std::cerr << "Error: --out is required\n";
		usage (argv[0]);
		return 1;
	}
	if (profile != "budget" && profile != "stress")
	{
		std::cerr << "Error: --profile must be 'budget' or 'stress'\n";
		return 1;
	}
	if (instances < 1)
	{
		std::cerr << "Error: --instances must be >= 1\n";
		return 1;
	}

	std::vector<ConfigSpec> sweep =
			(profile == "stress" ? kStressSweep : kBudgetSweep);

	if (profile == "budget")
	{
		// Split the budget evenly over the sweep and invert
		//     curves ~= (N-1) * discr^2 * (nref+1)
		// to get the discretisation each configuration can afford.
		const std::uint64_t per_instance =
				budget / (static_cast<std::uint64_t> (sweep.size()) * instances);
		for (auto& cfg : sweep)
		{
			// Invert  rounds * (N-3) * discr^2 ~= per_instance. With N == 3 there is
			// no quadratic stage at all (both endpoints are fixed), so the cost grows
			// only linearly with discr and the budget cannot be spent on it; cap the
			// discretisation instead of blowing it up to five digits.
			const double rounds = static_cast<double> (cfg.nref + 1);
			if (cfg.n_points <= 3) { cfg.discr = kMaxBudgetDiscr; }
			else
			{
				const double denom = static_cast<double> (cfg.n_points - 3) * rounds;
				int discr = static_cast<int> (std::lround (std::sqrt (per_instance / denom)));
				cfg.discr = std::min (kMaxBudgetDiscr, std::max (8, discr));
			}
		}
	}

	std::mt19937 gen (seed);
	std::vector<bench::Problem> problems;
	int id = 0;
	for (const auto& cfg : sweep)
	{
		for (int r = 0; r < instances; ++r)
		{
			problems.push_back (makeProblem (id++, cfg, gen));
		}
	}

	std::uint64_t projected = 0;
	for (const auto& p : problems) { projected += p.projectedCurves(); }

	bench::saveProblems (out_path, problems);

	std::cout << "profile          : " << profile << "\n"
						<< "instances        : " << problems.size() << "\n"
						<< "projected curves : " << projected << "\n"
						<< "written to       : " << out_path << "\n\n";
	std::cout << "  tag       N   discr  nref      curves\n";
	for (const auto& p : problems)
	{
		std::printf (
				"  %-6s %4zu  %5d  %4d  %10llu\n", p.tag.c_str(), p.size(), p.discr, p.nref,
				(unsigned long long)p.projectedCurves());
	}
	return 0;
}
