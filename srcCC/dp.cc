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
			points[i - 1].x(), points[i - 1].y(), points[i].x(), points[i].y(), 1. / Kmax, XC,
			YC);
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
	int bestIdx = -1;
	LEN_T bestL = std::numeric_limits<LEN_T>::max();
	// Find best path overall
	for (uint i = 0; i < MATRIX[0].size(); ++i)
	{
		if (MATRIX[0][i].l() < bestL)
		{
			bestL		= MATRIX[0][i].l();
			bestIdx = (int)(i);
		}
	}

	std::vector<Angle> vtheta;
	vtheta.push_back (MATRIX[0][bestIdx].th());
	for (int i = 0; i < points->size() - 1; ++i)
	{
		int nIdx = MATRIX[i][bestIdx].next();
		vtheta.push_back (MATRIX[i + 1][nIdx].th());
		bestIdx = nIdx;
	}

	if (points != nullptr)
	{
		for (int i = 0; i < points->size(); ++i) { (*points)[i].th (vtheta[i]); }
	}

	return std::pair<LEN_T, std::vector<Angle>> (bestL, vtheta);
}
#endif