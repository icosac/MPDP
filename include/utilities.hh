/*!
 * @file utilities.hh
 * @brief Utility functions
 *
 * @author Enrico Saccon <enricosaccon96[at]gmail[dot]com>
 * @copyright Copyright (c) 2020 Enrico Saccon. All rights reserved.
 *
 * @brief This file contains utility functions for printing numbers in scientific notation
 * and initializing Asymptote files.
 */
#ifndef UTILITIES_HH
#define UTILITIES_HH

#include <fstream>
#include <cstdlib>
#include <cmath>

inline void
PrintScientific1D (real_type d)
{
	if (d == 0)
	{
		printf ("%*d", 6, 0);
		return;
	}

	int exponent = (int)floor (log10 (std::abs (d)));	 // This will round down the exponent
	real_type base = d * pow (10, -1.0 * exponent);

	printf ("%1.1lfe%+01d", base, exponent);
}

inline void
PrintScientific2D (real_type d)
{
	if (d == 0)
	{
		printf ("%*d", 7, 0);
		return;
	}

	int exponent = (int)floor (log10 (std::abs (d)));	 // This will round down the exponent
	real_type base = d * pow (10, -1.0 * exponent);

	printf ("%1.1lfe%+02d", base, exponent);
}

template<typename IntType = uint64_t>
std::string PrintScientificLargeInt(IntType num){
	std::stringstream sstream;
	sstream << std::scientific << (double)num;
	return std::string(sstream.str());
}

inline void
initAsyFile (std::ofstream& file)
{
	file << "import graph;\n"
			 << "include \"clothoidLib.asylib\";\n"
			 << "size(14cm,7cm);\n"
			 << "\n\n\n"
			 << "path p;\n";
}

#endif	// UTILITIES_HH
