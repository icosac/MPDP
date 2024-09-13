/**
* @file math_utils.cc
* @author Enrico Saccon <enricosaccon96@gmail.com>
* @license This project is released under the GNU Public License Agero 3.0.
* @copyright Copyright 2020 Enrico Saccon. All rights reserved.
* @brief This file contains the source code for some of the mathematical utility functions.
*/

#include <math_utils.hh>

const real_type epsi				= std::numeric_limits<real_type>::epsilon();
const real_type m_pi				= 3.14159265358979323846264338328;
const real_type m_pi_2			= 1.57079632679489661923132169164;
const real_type m_2pi				= 6.28318530717958647692528676656;
const real_type m_1_pi			= 0.318309886183790671537767526745;
const real_type m_1_sqrt_pi = 0.564189583547756286948079451561;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Angle
mod2pi (Angle ang)
{
	while (ang < 0) ang += m_2pi;
	while (ang >= 2 * m_pi) ang -= m_2pi;
	return ang;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double
sinc (double x)
{
	if (std::abs (x) < 0.002)
	{
		double xs = x * x;
		return 1 - xs / 6. * (1 - xs / 20.0);
	}
	else { return std::sin (x) / x; }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double
f (double ell, double k, double th)
{
	double tmp = k * ell * 0.5;
	return ell * sinc (tmp) * std::cos (th + tmp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double
g (double ell, double k, double th)
{
	double tmp = k * ell * 0.5;
	return ell * sinc (tmp) * std::sin (th + tmp);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Configuration2
circleLine (double s, double kur, Configuration2 c)
{
	double xEnd, yEnd, thetaEnd;
	xEnd		 = c.x() - f (s, kur, mod2pi (c.th() + m_pi));
	yEnd		 = c.y() - g (s, kur, mod2pi (c.th() + m_pi));
	thetaEnd = mod2pi (c.th() + kur * s);

	Configuration2 res (xEnd, yEnd, thetaEnd, kur);

	return res;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////