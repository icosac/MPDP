/**
* @file math_utils.hh
* @author Enrico Saccon <enricosaccon96@gmail.com>
* @license This project is released under the GNU Public License Agero 3.0.
* @copyright Copyright 2020 Enrico Saccon. All rights reserved.
* @brief This file contains various mathematical utility functions used in the code.
*/

#ifndef MPDP_MATH_UTILS_HH
#define MPDP_MATH_UTILS_HH

// Library includes
#include <typedefs.hh>
#include <configuration.hh>

// System includes
#include <cmath>
#include <limits>

extern real_type const epsi;				 ///< A small value to compare real numbers.
extern real_type const m_pi;				 ///< The value of \f$\displaystyle\pi\f$.
extern real_type const m_pi_2;			 ///< The value of \f$\displaystyle\frac{\pi}{2}\f$
extern real_type const m_2pi;				 ///< The value of \f$\displaystyle 2\times \pi\f$
extern real_type const m_1_pi;			 ///< The value of \f$\displaystyle\frac{1}{\pi}\f$
extern real_type const m_1_sqrt_pi;	 ///< The value of \f$\frac{1}{\sqrt{\pi}}\f$

/**
 * @brief Function to compare two real numbers.
 * @param x The first number.
 * @param y The second number.
 * @param EPSI The small difference that is allowed between the two numbers.
 * @return True if the two numbers are within the bounded difference.
 */
template <class T = double>
inline bool
eq (const T x, const T y, const T EPSI = std::numeric_limits<T>::epsilon()) {
	return ABS (x, y) <= (EPSI);
}

/*!
 * Function to standardize an angle between 0 and \f$2\pi\f$.
 * @param ang The angle to be standardized.
 * @return The standardized angle.
 */
Angle
mod2pi (Angle ang);

/**
 * @brief Function to compute the sinc of a value \f$\displaystyle\frac{sin(x)}{x}\f$.
 * @param x The value to compute the sinc function.
 * @return The value of the sinc function.
 */
double
sinc (double x);

double
f (double ell, double k, double th);

double
g (double ell, double k, double th);

/**
 * @brief Starting from a configuration, the function computes the arrival configuration
 * on a circle.
 * @param s The length of the circle line.
 * @param kur The curvature of the circle line.
 * @param c The initial configuration.
 * @return The final configuration.
 */
Configuration2
circleLine (double s, double kur, Configuration2 c);

#endif	// MPDP_MATH_UTILS_HH
