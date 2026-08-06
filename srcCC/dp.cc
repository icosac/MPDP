/**
 * @file dp.cc
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief This file contains the source code for some functions for the dynamic
 * programming algorithm.
 */

#ifndef CUDA_ON
#include <dp.hh>
#include <stdexcept>

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
	double TOL = 1e-8;

	double q	= std::hypot (x2 - x1, y2 - y1);
	double x3 = 0.5 * (x1 + x2);
	double y3 = 0.5 * (y1 + y2);

	double delta = r * r - q * q / 4.;

	XC.clear();
	YC.clear();

	if (delta < -TOL) return;

	if (delta < TOL)
	{
		XC.push_back (x3);
		YC.push_back (y3);
	}
	else
	{
		double deltaS = std::sqrt (delta);
		XC.push_back (x3 + deltaS * (y1 - y2) / q);
		YC.push_back (y3 + deltaS * (x2 - x1) / q);
		XC.push_back (x3 - deltaS * (y1 - y2) / q);
		YC.push_back (y3 - deltaS * (x2 - x1) / q);
	}
}

void
DP::guessInitialAngles (
		const uint i,
		std::vector<DP::Cell>& thPrev,
		std::vector<DP::Cell>& thCur,
		const std::vector<Configuration2>& points)
{
	thPrev.clear();
	thCur.clear();

	thPrev.reserve (5);
	thCur.reserve (5);

	// aligned on straight line
	double th =
			std::atan2 (points[i].y() - points[i - 1].y(), points[i].x() - points[i - 1].x());
	thPrev.emplace_back (th);
	thCur.emplace_back (DP::Cell (th));

	// aligned on circle
	std::vector<double> XC, YC;
	circles (
			points[i - 1].x(), points[i - 1].y(), points[i].x(), points[i].y(),
			1. / this->k_max_, XC, YC);
	for (uint j = 0; j < XC.size(); ++j)
	{
		th = std::atan2 (points[i - 1].y() - YC[j], points[i - 1].x() - XC[j]);
		thPrev.emplace_back (DP::Cell (th + m_pi / 2.));
		thPrev.emplace_back (DP::Cell (th - m_pi / 2.));
		th = std::atan2 (points[i].y() - YC[j], points[i].x() - XC[j]);
		thCur.emplace_back (DP::Cell (th + m_pi / 2.));
		thCur.emplace_back (DP::Cell (th - m_pi / 2.));
	}
}

void
DP::setSamplingAngles (
		int discr,
		const std::vector<bool>& fixedAngles,
		const std::vector<Configuration2>& points)
{
	MATRIX.clear();
	MATRIX.resize (points.size());

	Angle dtheta = 2 * m_pi / discr;
	for (uint i = 0; i < points.size(); ++i)
	{
		MATRIX.reserve (discr + 10);
		for (int j = 0; j < discr; ++j)
		{
			MATRIX[i].push_back (DP::Cell (points[i].th() + dtheta * j));
		}
		if (i > 0)
		{
			std::vector<DP::Cell> thPrev, thCur;
			guessInitialAngles (i, thPrev, thCur, points);
			MATRIX[i - 1].insert (MATRIX[i - 1].end(), thPrev.begin(), thPrev.end());
			MATRIX[i].insert (MATRIX[i].end(), thCur.begin(), thCur.end());
		}
	}
	for (uint i = 0; i < points.size(); ++i)
	{
		if (fixedAngles[i])
		{
			MATRIX[i].clear();
			MATRIX[i] = {DP::Cell (points[i].th())};
		}
	}
}

void
DP::setSamplingAngles (
		const std::vector<Configuration2>& points,
		const std::vector<bool>& fixedAngles,
		double hrange,
		int hn)
{
	MATRIX.clear();
	MATRIX.resize (points.size());

	Angle dtheta = hrange / hn;
	for (uint i = 0; i < points.size(); ++i)
	{
		MATRIX.reserve (
				2 * hn + 11);	 // up to 10 "special" + one for thref + hn on both side of thref
		MATRIX[i].push_back (DP::Cell (points[i].th()));
		for (int j = 1; j <= hn; ++j)
		{
			MATRIX[i].push_back (DP::Cell (dtheta * j + points[i].th()));
			MATRIX[i].push_back (DP::Cell (-dtheta * j + points[i].th()));
		}
		if (i >= 1)
		{
			std::vector<DP::Cell> thPrev, thCur;
			guessInitialAngles (i, thPrev, thCur, points);
			MATRIX[i - 1].insert (MATRIX[i - 1].end(), thPrev.begin(), thPrev.end());
			MATRIX[i].insert (MATRIX[i].end(), thCur.begin(), thCur.end());
		}
	}
	for (uint i = 0; i < points.size(); ++i)
	{
		if (fixedAngles[i])
		{
			MATRIX[i].clear();
			MATRIX[i] = {DP::Cell (points[i].th())};
		}
	}
}

std::pair<LEN_T, std::vector<Angle>>
DP::bestAngles (std::vector<Configuration2>* points)
{
	if (MATRIX.empty() || MATRIX[0].empty())
	{
		throw std::runtime_error ("DP::bestAngles: empty DP matrix");
	}

	int bestIdx = -1;
	LEN_T bestL = std::numeric_limits<LEN_T>::max();
	// Find best path overall
	for (uint i = 0; i < MATRIX[0].size(); ++i)
	{
		if (MATRIX[0][i].l() < bestL)
		{
			bestL	= MATRIX[0][i].l();
			bestIdx = (int)(i);
		}
	}

	if (bestIdx == -1 || bestL == std::numeric_limits<LEN_T>::max())
	{
		throw std::runtime_error ("DP::bestAngles: no feasible path found");
	}

	const size_t numPoints = points ? points->size() : MATRIX.size();
	if (numPoints == 0 || MATRIX.size() < numPoints)
	{
		throw std::runtime_error ("DP::bestAngles: inconsistent points size");
	}

	std::vector<Angle> vtheta;
	vtheta.push_back (MATRIX[0][bestIdx].th());
	for (size_t i = 0; i + 1 < numPoints; ++i)
	{
		int nIdx = MATRIX[i][bestIdx].next();
		if (nIdx < 0 || nIdx >= (int)MATRIX[i + 1].size())
		{
			throw std::runtime_error ("DP::bestAngles: broken DP path");
		}
		vtheta.push_back (MATRIX[i + 1][nIdx].th());
		bestIdx = nIdx;
	}

	if (points != nullptr)
	{
		for (size_t i = 0; i < numPoints; ++i) { (*points)[i].th (vtheta[i]); }
	}

	return std::pair<LEN_T, std::vector<Angle>> (bestL, vtheta);
}

std::vector<std::pair<size_t, size_t>>
DP::bestPathIndices() const
{
	std::vector<std::pair<size_t, size_t>> indices;
	if (MATRIX.empty()) { return indices; }

	int bestIdx = -1;
	LEN_T bestL = std::numeric_limits<LEN_T>::max();
	for (size_t col = 0; col < MATRIX[0].size(); ++col)
	{
		if (MATRIX[0][col].l() < bestL)
		{
			bestL	 = MATRIX[0][col].l();
			bestIdx = static_cast<int> (col);
		}
	}

	if (bestIdx < 0) { return indices; }

	for (size_t row = 0; row < MATRIX.size() && bestIdx >= 0; ++row)
	{
		indices.emplace_back (row, static_cast<size_t> (bestIdx));
		if (row + 1 >= MATRIX.size()) { break; }
		bestIdx = MATRIX[row][bestIdx].next();
	}

	return indices;
}

void
DP::exportVisualizationData (const std::string& json_path) const
{
	if (!this->has_solution_)
	{
		throw std::runtime_error (
			"DP::exportVisualizationData: solveDP must finish successfully before exporting");
	}

	if (this->last_points_.size() != MATRIX.size())
	{
		throw std::runtime_error (
			"DP::exportVisualizationData: internal cache is inconsistent with DP matrix");
	}

	std::ofstream out (json_path);
	if (!out.is_open())
	{
		throw std::runtime_error ("DP::exportVisualizationData: unable to open output path");
	}
	out << std::setprecision (16);

	const LEN_T invalid_length = std::numeric_limits<LEN_T>::max() / 4.0;
	auto write_number = [] (std::ostream& stream, double value) {
		if (std::isfinite (value)) { stream << value; }
		else { stream << "null"; }
	};
	auto write_length = [invalid_length] (std::ostream& stream, double value) {
		if (!std::isfinite (value) || value >= invalid_length) { stream << "null"; }
		else { stream << value; }
	};

	const auto best_path = bestPathIndices();
	out << "{\n";
	out << "  \"k_max\": ";
	write_number (out, this->k_max_);
	out << ",\n";
	out << "  \"points\": [\n";
	for (size_t i = 0; i < this->last_points_.size(); ++i)
	{
		const auto& pt = this->last_points_[i];
		out << "    [";
		write_number (out, pt.x());
		out << ", ";
		write_number (out, pt.y());
		out << "]";
		if (i + 1 < this->last_points_.size()) { out << ","; }
		out << "\n";
	}
	out << "  ],\n";
	out << "  \"best_angles\": [";
	for (size_t i = 0; i < this->last_best_angles_.size(); ++i)
	{
		write_number (out, this->last_best_angles_[i]);
		if (i + 1 < this->last_best_angles_.size()) { out << ", "; }
	}
	out << "],\n";
	out << "  \"matrix\": [\n";
	for (size_t row = 0; row < MATRIX.size(); ++row)
	{
		out << "    [\n";
		const auto& cells = MATRIX[row];
		for (size_t col = 0; col < cells.size(); ++col)
		{
			const auto& cell = cells[col];
			out << "      {\"theta\": ";
			write_number (out, cell.th());
			out << ", \"length\": ";
			write_length (out, cell.l());
			out << ", \"next_col\": ";
			if (cell.next() >= 0) { out << cell.next(); }
			else { out << "null"; }
			out << "}";
			if (col + 1 < cells.size()) { out << ","; }
			out << "\n";
		}
		out << "    ]";
		if (row + 1 < MATRIX.size()) { out << ","; }
		out << "\n";
	}
	out << "  ],\n";
	out << "  \"best_path\": [\n";
	for (size_t i = 0; i < best_path.size(); ++i)
	{
		out << "    [" << best_path[i].first << ", " << best_path[i].second << "]";
		if (i + 1 < best_path.size()) { out << ","; }
		out << "\n";
	}
	out << "  ]\n";
	out << "}\n";
}
#endif
