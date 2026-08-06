/**
 * @file bench_io.hh
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @brief Dependency-free problem/result I/O shared by every solver driver.
 *
 * This header is deliberately free of any MPDP dependency: it only uses the
 * standard library, so the CPU driver (compiled against `srcCC`) and the GPU
 * driver (compiled against `srcCU`) can exchange problems and results as plain
 * text and be guaranteed to see byte-identical inputs.
 *
 * The two libraries could now be linked into one binary - they live in
 * `mpdp::cpu` and `mpdp::gpu` - but keeping the drivers as separate processes
 * is still worth it: it keeps CUDA context creation out of the CPU timings,
 * and a crash in one solver leaves the other's results intact.
 */

#ifndef MPDP_UPDATE_TEST_BENCH_IO_HH
#define MPDP_UPDATE_TEST_BENCH_IO_HH

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace bench {

/*!
 * A single multi-point Dubins instance.
 *
 * The problem is fully defined by the initial direction, the final direction
 * and the N intermediate points; `kmax`, `discr` and `nref` describe how the
 * DP is asked to solve it and are part of the instance so that every solver
 * receives byte-identical settings.
 */
struct Problem {
	int id				 = 0;			///< Unique instance id, matches across result files.
	std::string tag = "";		///< Name of the configuration this instance belongs to.
	double kmax		 = 3.0;		///< Maximum curvature.
	int discr			 = 90;		///< Number of angle samples per point.
	int nref			 = 4;			///< Number of refinement rounds.
	double th_init = 0.0;		///< Fixed direction at the first point.
	double th_final = 0.0;	///< Fixed direction at the last point.
	std::vector<double> x;	///< Abscissas of the N points.
	std::vector<double> y;	///< Ordinates of the N points.

	//! Number of points of the instance.
	std::size_t
	size() const
	{
		return x.size();
	}

	/*!
	 * Count of the point-to-point curves the DP has to evaluate for this
	 * instance, ignoring the handful of extra "guess" angles per row.
	 *
	 * The first and the last angle are fixed, so those two rows of the DP matrix
	 * hold a single cell each. A stage between two free rows costs discr^2, but
	 * the two stages touching an endpoint only cost discr:
	 *     rounds * ( 2*discr + (N-3)*discr^2 )
	 * @return The projected number of evaluated curves.
	 */
	std::uint64_t
	projectedCurves() const
	{
		const std::int64_t n = static_cast<std::int64_t> (x.size());
		if (n < 2) { return 0; }
		const std::uint64_t d			 = static_cast<std::uint64_t> (discr);
		const std::uint64_t rounds = static_cast<std::uint64_t> (nref + 1);
		if (n == 2) { return rounds; }
		const std::uint64_t inner = static_cast<std::uint64_t> (n > 3 ? n - 3 : 0);
		return rounds * (2 * d + inner * d * d);
	}
};

//! Outcome of one solver on one instance.
struct RunResult {
	int id					= 0;						///< Instance id.
	std::string solver	= "";				///< Solver name.
	std::string status	= "ok";			///< `ok`, `skipped`, or `error:<msg>`.
	//! Length as reported by the solver's own DP bookkeeping.
	double length = std::numeric_limits<double>::quiet_NaN();
	//! Length recomputed from the returned angles with the reference CPU Dubins.
	//! A mismatch with `length` means the solver's DP accounting is broken.
	double length_check = std::numeric_limits<double>::quiet_NaN();
	double time_ms	= 0.0;	 ///< Wall-clock solve time in milliseconds.
	double time_solver_ms = 0.0;	 ///< Solver-internal time (GPU: kernels + copies).
	std::vector<double> angles;		 ///< The N angles of the solution.
};

namespace detail {

	//! Split `s` on `sep`.
	inline std::vector<std::string>
	split (const std::string& s, char sep)
	{
		std::vector<std::string> out;
		std::string item;
		std::istringstream iss (s);
		while (std::getline (iss, item, sep)) { out.push_back (item); }
		return out;
	}

	//! Join a vector of doubles with ';' at full precision.
	inline std::string
	joinDoubles (const std::vector<double>& v)
	{
		std::ostringstream oss;
		oss << std::setprecision (17);
		for (std::size_t i = 0; i < v.size(); ++i)
		{
			if (i != 0) { oss << ';'; }
			oss << v[i];
		}
		return oss.str();
	}

	/*!
	 * Makes a free-form string safe to drop into a CSV field: an exception
	 * message carrying a comma or a newline would otherwise shift every column
	 * after it.
	 * @param s The string to sanitise.
	 * @return The string with separators replaced.
	 */
	inline std::string
	csvSafe (std::string s)
	{
		for (char& c : s)
		{
			if (c == ',') { c = ';'; }
			else if (c == '\n' || c == '\r') { c = ' '; }
		}
		return s;
	}

	//! Parse a ';'-separated list of doubles.
	inline std::vector<double>
	parseDoubles (const std::string& s)
	{
		std::vector<double> out;
		if (s.empty()) { return out; }
		for (const auto& tok : split (s, ';'))
		{
			if (!tok.empty()) { out.push_back (std::stod (tok)); }
		}
		return out;
	}

}	 // namespace detail

/*!
 * Writes the problem set to disk.
 * @param path Destination file.
 * @param problems The instances to write.
 */
inline void
saveProblems (const std::string& path, const std::vector<Problem>& problems)
{
	std::ofstream out (path);
	if (!out.is_open())
	{
		throw std::runtime_error ("bench::saveProblems: cannot open " + path);
	}
	out << std::setprecision (17);
	out << "id,tag,n_points,kmax,discr,nref,th_init,th_final,x,y\n";
	for (const auto& p : problems)
	{
		out << p.id << ',' << p.tag << ',' << p.x.size() << ',' << p.kmax << ','
				<< p.discr << ',' << p.nref << ',' << p.th_init << ',' << p.th_final << ','
				<< detail::joinDoubles (p.x) << ',' << detail::joinDoubles (p.y) << '\n';
	}
}

/*!
 * Reads back a problem set written by `saveProblems`.
 * @param path Source file.
 * @return The instances, in file order.
 */
inline std::vector<Problem>
loadProblems (const std::string& path)
{
	std::ifstream in (path);
	if (!in.is_open())
	{
		throw std::runtime_error ("bench::loadProblems: cannot open " + path);
	}

	std::vector<Problem> problems;
	std::string line;
	bool first = true;
	while (std::getline (in, line))
	{
		if (line.empty()) { continue; }
		if (first)
		{
			first = false;
			if (line.rfind ("id,", 0) == 0) { continue; }	 // header
		}
		const auto f = detail::split (line, ',');
		if (f.size() < 10)
		{
			throw std::runtime_error ("bench::loadProblems: malformed line: " + line);
		}
		Problem p;
		p.id			 = std::stoi (f[0]);
		p.tag			 = f[1];
		p.kmax		 = std::stod (f[3]);
		p.discr		 = std::stoi (f[4]);
		p.nref		 = std::stoi (f[5]);
		p.th_init	 = std::stod (f[6]);
		p.th_final = std::stod (f[7]);
		p.x				 = detail::parseDoubles (f[8]);
		p.y				 = detail::parseDoubles (f[9]);
		if (p.x.size() != p.y.size() || p.x.size() != (std::size_t)std::stoul (f[2]))
		{
			throw std::runtime_error ("bench::loadProblems: inconsistent point count");
		}
		problems.push_back (std::move (p));
	}
	return problems;
}

/*!
 * Opens a result file and writes its header.
 * @param path Destination file.
 * @return The open stream.
 */
inline std::ofstream
openResults (const std::string& path)
{
	std::ofstream out (path);
	if (!out.is_open())
	{
		throw std::runtime_error ("bench::openResults: cannot open " + path);
	}
	out << "id,solver,status,length,length_check,time_ms,time_solver_ms,angles\n";
	return out;
}

/*!
 * Appends one result row and flushes, so a crash still leaves usable data.
 * @param out The stream returned by `openResults`.
 * @param r The result to write.
 */
inline void
appendResult (std::ofstream& out, const RunResult& r)
{
	out << std::setprecision (17);
	out << r.id << ',' << r.solver << ',' << detail::csvSafe (r.status) << ',' << r.length
			<< ','
			<< r.length_check << ',' << r.time_ms << ',' << r.time_solver_ms << ','
			<< detail::joinDoubles (r.angles) << '\n';
	out.flush();
}

}	 // namespace bench

#endif	// MPDP_UPDATE_TEST_BENCH_IO_HH
