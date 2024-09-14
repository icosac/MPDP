/*!
 * \file utils.hh
 * \brief Utility functions
 * \author Enrico Saccon <enricosaccon96[at]gmail[dot]com>
 * \copyright Copyright (c) 2020 Enrico Saccon. All rights reserved.
 * \brief This file contains utility functions of different nature, plus it provides the includes
 * for the other utility functions.
 */
#ifndef UTILS_HH
#define UTILS_HH

// System imports
#include <sstream>

// Library imports to expose
#include <typedefs.hh>
#include <utilities.hh>
#include <math_utils.hh>

// #define DEBUG

#ifdef DEBUG
#define COUT(x) std::cout << #x << ": " << x << std::endl;
#define printCV(v, d)                                    \
	printf ("<");                                          \
	for (uint i = 0; i < d; i++) { printf ("%s ", v[i]); } \
	printf ("\n");

#define printV(v)                         \
	std::cout << "<";                       \
	for (auto a : v) std::cout << a << " "; \
	std::cout << ">" << endl;

#define printM(M, discr, size)                                  \
	for (int i = 0; i < discr; i++)                               \
	{                                                             \
		cout << "th" << i;                                          \
		for (int j = 0; j < size; j++) { cout << M[i][j] << "\t"; } \
		cout << endl;                                               \
	}

#define printVM(M, discr, size)             \
	for (int i = 0; i < discr; i++)           \
	{                                         \
		std::cout << "th" << i;                 \
		for (int j = 0; j < size; j++)          \
		{                                       \
			std::cout << std::setw (30);          \
			std::cout << M[i * size + j] << "\t"; \
		}                                       \
		std::cout << std::endl;                 \
	}

#define printCVM(M, discr, size)                                                     \
	for (int i = 0; i < discr; i++)                                                    \
	{                                                                                  \
		printf ("th%d", i);                                                              \
		for (int j = 0; j < size; j++) { printf ("\t%-5f", (double)(M[i * size + j])); } \
		printf ("\n");                                                                   \
	}
#else
#define COUT(x)
#define printCV(v, d)
#define printV(v)
#define printM(M, discr, size)
#define printVM(M, discr, size)
#define printCVM(M, discr, size)
#endif	// DEBUG

#ifndef ASSERT
#define ASSERT(COND, MSG)                                                   \
	if (!(COND))                                                              \
	{                                                                         \
		std::ostringstream ost;                                                 \
		ost << "On line: " << __LINE__ << " file: " << __FILE__ << MSG << '\n'; \
		throw std::runtime_error (ost.str());                                   \
	}
#endif	// ASSERT

#endif	// UTILS_HH
