/**
 * @file IOUtils.hh
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License Agero 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Header file containing function to read the points from a file.
 * TODO This should be really improved.
 */

#ifndef IOUTILS_HH
#define IOUTILS_HH

#ifndef CUDA_ON
#include <configuration.hh>
#else
#include <configuration.cuh>
#endif

// System includes
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <cstdlib>

// Library includes
#include <asyplot.hh>

std::vector<std::string>
splitString (std::string str)
{
	size_t pos1 = 0;
	std::vector<std::string> ret;
	while (pos1 != std::string::npos)
	{
		pos1 = str.find (" ");
		// std::cout << "[" << str << "] " << pos1 << " " << std::string::npos << std::endl;
		if (pos1 != std::string::npos)
		{
			ret.emplace_back (str.substr (0, pos1));
			str = str.substr (pos1 + 1, str.length());
		}
		else
		{
			std::string str1 = str.substr (0, str.length());
			if (str1.length() > 0) { ret.emplace_back (str1); }
		}
	}
	return ret;
}

/*!
 * Function to read from a file the `Configuration2`.
 * @param inputFile A `fstream` to the input file.
 * @param type The type of file: 0 or 1. Default to 0.
 * @param close Whether to close the file after reading or not. Default to true.
 * @return Returns a pair where the first element is vector of `Configuration2` and the
 * second element is a vector of booleans for the fixed angles.
 */
std::pair<std::vector<Configuration2>, std::vector<bool>>
readConfigurationsFromFile (std::fstream& inputFile, int type = 0, bool close = true)
{
	std::vector<Configuration2> points;
	std::vector<bool> fixedAngles;

	if (type == 0)
	{
		real_type x, y;
		unsigned short fixed;
		std::string th, appString;

		while (getline (inputFile, appString))
		{
			int i = 0;
			for (auto el : splitString (appString))
			{
				switch (i)
				{
					case 0:
						x = atof (el.c_str());
						break;
					case 1:
						y = atof (el.c_str());
						break;
					case 2: {
						th = el;
						if (th.find ("FREE") != std::string::npos)
						{
							points.emplace_back (Configuration2 (x, y, ANGLE::FREE));
							fixedAngles.emplace_back (false);
						}
						else { points.emplace_back (Configuration2 (x, y, atof (th.c_str()))); }
						break;
					}
					case 3: {
						if (el == "0") { fixedAngles.emplace_back (false); }
						else { fixedAngles.emplace_back (true); }
					}
				}
				i++;
			}
		}
	}
	if (type == 1)
	{
		uint counter = 0;
		size_t dim;
		real_type value;
		std::string th;
		int fixed;

		inputFile >> dim;
		points.resize (dim);

		while (inputFile >> value)
		{
			if ((int)(counter / dim) == 0) { points[counter % dim].x (value); }
			else if ((int)(counter / dim) == 1) { points[counter % dim].y (value); }
			else { points[counter % dim].th (value); }
			counter++;
			if ((int)(counter / dim) == 2) { break; }
		}
		counter = 0;
		while (inputFile >> th)
		{
			if (th.find ("FREE") != std::string::npos) { points[counter].th (ANGLE::FREE); }
			else { points[counter].th (atof (th.c_str())); }
			counter++;
			if ((int)(counter / dim) == 1) { break; }
		}
		while (inputFile >> fixed)
		{
			if (fixed == 0) { fixedAngles.emplace_back (false); }
			else { fixedAngles.emplace_back (true); }
		}
	}

	if (close) { inputFile.close(); }

	return std::pair<std::vector<Configuration2>, std::vector<bool>> (points, fixedAngles);
}

/*!
 * Function to read from a file the `Configuration2`.
 * @param inputFile A string containing the path to the input file.
 * @param type The type of file: 0 or 1. Default to 0.
 * @return Returns a pair where the first element is vector of `Configuration2` and the
 * second element is a vector of booleans for the fixed angles.
 */
std::pair<std::vector<Configuration2>, std::vector<bool>>
readConfigurationsFromFile (const char* filename, int type = 0)
{
	std::fstream inputFile;
	inputFile.open (filename, std::fstream::in);
	return readConfigurationsFromFile (inputFile, type, true);
}

/*!
 * Function to read from a file the points for the `Configuration2`. The first and final
 * angle are manually set
 * @param inputFile A `fstream` to the input file.
 * @param thi The initial angle to be set. Default is ANGLE::FREE.
 * @param thf The initial angle to be set. Default is ANGLE::FREE.
 * @param type The type of file: 0 or 1. Default to 0.
 * @param close Whether to close the file after reading or not. Default to true.
 * @return Returns a pair where the first element is vector of `Configuration2` and the
 * second element is a vector of booleans for the fixed angles.
 */
std::pair<std::vector<Configuration2>, std::vector<bool>>
readPointsFromFile (
		std::fstream& inputFile,
		Angle thi	 = ANGLE::FREE,
		Angle thf	 = ANGLE::FREE,
		int type	 = 0,
		bool close = true)
{
	std::vector<Configuration2> points;

	if (type == 0)
	{
		real_type x, y;
		while (inputFile >> x >> y)
		{
			points.emplace_back (Configuration2 (x, y, ANGLE::FREE));
		}
	}

	if (type == 1)
	{
		uint counter = 0;
		size_t dim;
		real_type value;

		inputFile >> dim;
		points.resize (dim);

		while (inputFile >> value)
		{
			if (counter < dim) { points[counter % dim].x (value); }
			else { points[counter % dim].y (value); }
			counter++;
		}
	}

	points.front().th (thi);
	points.back().th (thf);

	if (close) { inputFile.close(); }

	std::vector<bool> fixedAngles (points.size(), false);
	fixedAngles.front() = true;
	fixedAngles.back()	= true;

	return std::pair<std::vector<Configuration2>, std::vector<bool>> (points, fixedAngles);
}

/*!
 * Function to read from a file the points for the `Configuration2`. The first and final
 * angle are manually set
 * @param inputFile A string containing the path to the input file.
 * @param thi The initial angle to be set. Default is ANGLE::FREE.
 * @param thf The initial angle to be set. Default is ANGLE::FREE.
 * @param type The type of file: 0 or 1. Default to 0.
 * @return Returns a pair where the first element is vector of `Configuration2` and the
 * second element is a vector of booleans for the fixed angles.
 */
std::pair<std::vector<Configuration2>, std::vector<bool>>
readPointsFromFile (
		const char* filename, Angle thi = ANGLE::FREE, Angle thf = ANGLE::FREE, int type = 0)
{
	std::fstream inputFile;
	inputFile.open (filename, std::fstream::in);
	return readPointsFromFile (inputFile, thi, thf, type, true);
}

/*!
 * This function returns informations regarding the Dubins composing the path.
 * @param points A vector of configurations with or without the newly computed angles.
 * @param kmax The maximum curvature used to compute the best angles.
 * @param vtheta The best angles. Default is nullptr, instead if set, the function will
 * assign the angles to the point before starting.
 * @param len The length the MPMD should have. Deafult is 0.0, instead if set, it's used
 * to check whether the now computed Dubins are correct.
 */
void
getMPMDInfo (
		std::vector<Configuration2> points,
		K_T kmax,
		const std::vector<Angle>* vtheta = NULL,
		LEN_T len												 = 0.0)
{
	std::vector<Dubins> dubinss;
	LEN_T Len = 0.0;
	for (uint i = 0; i < points.size() - 1; i++)
	{
		if (vtheta != NULL)
		{
			if (i == 0) { points[i].th ((*vtheta)[i]); }
			points[i + 1].th ((*vtheta)[i + 1]);
		}
		dubinss.emplace_back (Dubins (points[i], points[i + 1], kmax));
		Len += dubinss.back().l();
	}
	if (len != 0 && !eq<double> (len, Len, 1e-13))
	{
		std::cout << "ERROR before: " << std::setprecision (20) << len << " after: " << Len
							<< " difference: ";
		PrintScientific2D (std::abs (len - Len));
		std::cout << std::endl;
	}
	else { std::cout << "Total length: " << Len << std::endl; }
	for (Dubins d : dubinss) { std::cout << d << std::endl; }
}

#endif	// IOUTILS_HH