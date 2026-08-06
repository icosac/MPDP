/**
 * @file bench_driver.hh
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @brief Common main-loop for the solver drivers.
 *
 * Every driver reads the same problem file, walks it in the same order and
 * honours the same wall-clock cap, so the three result files line up row by
 * row. Instances not reached before the cap are written out as `skipped`
 * rather than dropped, which keeps the comparison explicit about what a slow
 * solver did not get to.
 */

#ifndef MPDP_UPDATE_TEST_BENCH_DRIVER_HH
#define MPDP_UPDATE_TEST_BENCH_DRIVER_HH

#include <chrono>
#include <cmath>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include <bench_io.hh>

namespace bench {

//! What a solver returns for one instance.
struct SolveOutput {
	double length = 0.0;					///< Length as the solver's DP accounts for it.
	std::vector<double> angles;		///< One angle per point.
	double time_solver_ms = 0.0;	///< Solver-internal time, 0 if not measured.
};

//! Options shared by every driver.
struct DriverOptions {
	std::string problems_path;
	std::string out_path;
	double max_seconds = 60.0;	///< Hard cap on the whole run, per solver.
	int reps					 = 1;			///< Timing repetitions; the fastest one is reported.
};

/*!
 * Parses the options every driver understands. Unknown arguments are left
 * alone so a driver can add its own.
 * @param argc Argument count.
 * @param argv Argument values.
 * @param opts Where to store the parsed options.
 * @return `true` on success.
 */
inline bool
parseCommonOptions (int argc, char** argv, DriverOptions& opts)
{
	for (int i = 1; i < argc; ++i)
	{
		const std::string a = argv[i];
		const bool has_next = (i + 1 < argc);
		if (a == "--problems" && has_next) { opts.problems_path = argv[++i]; }
		else if (a == "--out" && has_next) { opts.out_path = argv[++i]; }
		else if (a == "--max-seconds" && has_next)
		{
			opts.max_seconds = std::stod (argv[++i]);
		}
		else if (a == "--reps" && has_next) { opts.reps = std::stoi (argv[++i]); }
	}
	if (opts.problems_path.empty() || opts.out_path.empty())
	{
		std::cerr << "Error: --problems <path> and --out <path> are required\n";
		return false;
	}
	if (opts.reps < 1) { opts.reps = 1; }
	return true;
}

/*!
 * Runs `solve` over the whole problem set and writes the result file.
 *
 * @tparam SolveFn `SolveOutput(const Problem&)`.
 * @tparam LengthFn `double(const Problem&, const std::vector<double>&)`, used
 *         to independently recompute the length of the returned angles.
 * @param opts The parsed common options.
 * @param solver_name Name recorded in every row.
 * @param solve The solver under test.
 * @param pathLength The independent length check.
 * @return 0 on success.
 */
template <class SolveFn, class LengthFn>
int
runDriver (
		const DriverOptions& opts,
		const std::string& solver_name,
		SolveFn solve,
		LengthFn pathLength)
{
	const auto problems = loadProblems (opts.problems_path);
	auto out						= openResults (opts.out_path);

	const auto t_run_start = std::chrono::steady_clock::now();
	auto elapsedSeconds		 = [&t_run_start] () {
		 return std::chrono::duration<double> (
							std::chrono::steady_clock::now() - t_run_start)
				 .count();
	};

	std::size_t solved = 0, skipped = 0, failed = 0;
	bool out_of_time = false;

	for (const auto& p : problems)
	{
		RunResult r;
		r.id		 = p.id;
		r.solver = solver_name;

		if (out_of_time || elapsedSeconds() >= opts.max_seconds)
		{
			out_of_time = true;
			r.status		= "skipped";
			appendResult (out, r);
			++skipped;
			continue;
		}

		try
		{
			SolveOutput best;
			double best_ms = 0.0;
			for (int rep = 0; rep < opts.reps; ++rep)
			{
				const auto t0 = std::chrono::steady_clock::now();
				SolveOutput o = solve (p);
				const double ms =
						std::chrono::duration<double, std::milli> (
								std::chrono::steady_clock::now() - t0)
								.count();
				if (rep == 0 || ms < best_ms)
				{
					best_ms = ms;
					best		= std::move (o);
				}
				// Never let the repetitions push us past the cap.
				if (elapsedSeconds() >= opts.max_seconds) { break; }
			}

			r.status			 = "ok";
			r.length			 = best.length;
			r.angles			 = best.angles;
			r.time_ms			 = best_ms;
			r.time_solver_ms = best.time_solver_ms;
			r.length_check = pathLength (p, best.angles);
			++solved;
		}
		catch (const std::exception& e)
		{
			r.status = std::string ("error:") + e.what();
			++failed;
		}

		appendResult (out, r);

		std::cout << "[" << solver_name << "] id=" << r.id << " tag=" << p.tag
							<< " N=" << p.size() << " discr=" << p.discr << " status=" << r.status;
		if (r.status == "ok")
		{
			std::cout << " len=" << r.length << " check=" << r.length_check
								<< " t=" << r.time_ms << "ms";
		}
		std::cout << std::endl;
	}

	std::cout << "[" << solver_name << "] done: " << solved << " solved, " << skipped
						<< " skipped (wall-clock cap " << opts.max_seconds << "s), " << failed
						<< " failed, total " << elapsedSeconds() << "s" << std::endl;
	return failed == 0 ? 0 : 2;
}

}	 // namespace bench

#endif	// MPDP_UPDATE_TEST_BENCH_DRIVER_HH
