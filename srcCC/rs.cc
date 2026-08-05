/**
 * @file rs.cc
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief This file contains the source code for some of the functions to compute the P2P
 * RS.
 */

#ifndef CUDA_ON

#include <rs.hh>

#define RADCURVMUL2 (2 * RADCURV)
#define RADCURVMUL4 (4 * RADCURV)
#define SQRADCURV (RADCURV * RADCURV)
#define SQRADCURVMUL2 (4 * RADCURV * RADCURV)

const double EPS1 = 1.0e-14;
const double EPS3 = 1.0e-14;
const double EPS4 = 1.0e-14;

const double MPI	 = 3.1415926535897932385;
const double MPIMUL2 = 6.2831853071795864770;
const double MPIDIV2 = 1.5707963267948966192;

Configuration2
RS::circleLine (double s, double dir, double kur, double kmax, Configuration2 c, int seg)
{
	double sigmaDir, sign;
	if (dir > 0)
	{
		sigmaDir = 0;
		sign		 = 1;
	}
	else
	{
		sigmaDir = 1;
		sign		 = -1;
	}
	double xEnd, yEnd, thetaEnd;
	xEnd		 = c.x() + f (s, sign * kur * kmax, mod2pi (c.th() + sigmaDir * MPI));
	yEnd		 = c.y() + g (s, sign * kur * kmax, mod2pi (c.th() + sigmaDir * MPI));
	thetaEnd = mod2pi (c.th() + sign * kur * kmax * s);
	if (seg == 0)
	{
		X[seg]	= c.x();
		Y[seg]	= c.y();
		TH[seg] = c.th();
	}
	this->X[seg + 1]	= xEnd;
	this->Y[seg + 1]	= yEnd;
	this->TH[seg + 1] = thetaEnd;
	// this->L[seg] = s;
	this->L[seg] = sign * s;
	this->D[seg] = sign;
	this->K[seg] = kur * kmax;
	Configuration2 res (xEnd, yEnd, thetaEnd, kur * kmax);
	if (kur > 0)
	{
		if (dir > 0) { this->ManType[seg] = "Lp"; }
		else { this->ManType[seg] = "Ln"; }
	}
	else if (kur < 0)
	{
		if (dir > 0) { this->ManType[seg] = "Rp"; }
		else { this->ManType[seg] = "Rn"; }
	}
	else
	{
		if (dir > 0) { this->ManType[seg] = "Sp"; }
		else { this->ManType[seg] = "Sn"; }
	}

	return res;
}

void
RS::buildRS (int man)
{
	// if (man == -1) { man = this->_Nman; }
	std::vector<double> dir;
	std::vector<double> kk;
	std::vector<double> ll = {_L1, _L2, _L3};
	Configuration2 c, d;
	int seg				= 0;
	double len		= 1;
	double R			= -1;
	double S			= 0;
	double F			= 1;
	double B			= -1;
	double lambda = 1 / _kmax;

	switch (man)
	{
		// C | C | C
		case 1:
		case 4:
			// LpRnLp
			this->_Nseg = 3;
			dir					= {F, B, F};			// fwd-back-fwd
			kk					= {len, R, len};	// len R len
			ll					= {_L1, _L2, _L3};
			//      std::cout << "-" << _L1 << " " << _L2 << " " << _L3 << std::endl;
			//      std::cout << "-" << _L1/lambda << " " << _L2/lambda << " " << _L3/lambda <<
			//      std::endl;
			break;

			//    case 2:
			//      // Substituted by 3
			//      // LnRpLn
			//      this->_Nseg = 3;
			//      dir = { B,F,B }; // back-fwd-back
			//      kk = {len, R, len}; // len R len
			//      ll = { _L1, _L2, _L3 };
			//      break;

		case 2:
		case 3:
			// RpLnRp
			this->_Nseg = 3;
			dir					= {F, B, F};		// fwd-back-fwd
			kk					= {R, len, R};	// R len R
			ll					= {_L1, _L2, _L3};
			break;

			//    case 4:
			//      // Substituted by 1
			//      // RnLpRn
			//      this->_Nseg = 3;
			//      dir = { B,F,B }; // back-fwd-back
			//      kk = {R, len, R }; // R len R
			//      ll = { _L1, _L2, _L3 };
			//      break;

			// C | C C
		case 5:
			// LpRnLn
			this->_Nseg = 3;
			dir					= {F, B, B};			// fwd-back-back
			kk					= {len, R, len};	// len R len
			ll					= {_L1, _L2, _L3};
			break;

		case 6:
			// LnRpLp
			this->_Nseg = 3;
			dir					= {B, F, F};			// back-fwd-fwd
			kk					= {len, R, len};	// len R len
			ll					= {_L1, _L2, _L3};
			break;

		case 7:
			// RpLnRn
			this->_Nseg = 3;
			dir					= {F, B, B};		// fwd-back-back
			kk					= {R, len, R};	//  R len R
			ll					= {_L1, _L2, _L3};
			break;

		case 8:
			// RnLpRp
			this->_Nseg = 3;
			dir					= {B, F, F};		// back-fwd-fwd
			kk					= {R, len, R};	//  R len R
			ll					= {_L1, _L2, _L3};
			break;

			// C S C
		case 9:
			// LpSpLp
			this->_Nseg = 3;
			dir					= {F, F, F};			// fwd-fwd-fwd
			kk					= {len, S, len};	//  len S len
			ll					= {_L1, _L2 / lambda, _L3};
			break;

		case 10:
			// RpSpRp
			this->_Nseg = 3;
			dir					= {F, F, F};	// fwd-fwd-fwd
			kk					= {R, S, R};	//  R S R
			ll					= {_L1, _L2 / lambda, _L3};
			break;

		case 11:
			// LnSnLn
			this->_Nseg = 3;
			dir					= {B, B, B};			// back-back-back
			kk					= {len, S, len};	//  len S len
			ll					= {_L1, _L2 / lambda, _L3};
			break;

		case 12:
			// RnSnRn
			this->_Nseg = 3;
			dir					= {B, B, B};	// back-back-back
			kk					= {R, S, R};	//  R S R
			ll					= {_L1, _L2 / lambda, _L3};
			break;

		case 13:
			// LpSpRp
			this->_Nseg = 3;
			dir					= {F, F, F};		// fwd-fwd-fwd
			kk					= {len, S, R};	//  len S R
			ll					= {_L1, _L2 / lambda, _L3};
			break;

		case 14:
			// RpSpLp
			this->_Nseg = 3;
			dir					= {F, F, F};		// fwd-fwd-fwd
			kk					= {R, S, len};	//  R S len
			ll					= {_L1, _L2 / lambda, _L3};
			break;

		case 15:
			// LnSnRn
			this->_Nseg = 3;
			dir					= {B, B, B};		// back-back-back
			kk					= {len, S, R};	//  len S R
			ll					= {_L1, _L2 / lambda, _L3};
			break;

		case 16:
			// RnSnLn
			this->_Nseg = 3;
			dir					= {B, B, B};		// back-back-back
			kk					= {R, S, len};	//  R S len
			ll					= {_L1, _L2 / lambda, _L3};
			break;

			// C C= | C= C
		case 17:
			// LpRpLnRn
			this->_Nseg = 4;
			dir					= {F, F, B, B};			 // fwd-fwd-back-back
			kk					= {len, R, len, R};	 // len R len R
			ll					= {_L1, _L2, _L2, _L3};
			break;

		case 18:
			// RpLpRnLn
			this->_Nseg = 4;
			dir					= {F, F, B, B};			 // fwd-fwd-back-back
			kk					= {R, len, R, len};	 // R len R len
			ll					= {_L1, _L2, _L2, _L3};
			break;

		case 19:
			// LnRnLpRp
			this->_Nseg = 4;
			dir					= {B, B, F, F};			 // back-back-fwd-fwd
			kk					= {len, R, len, R};	 // len R len R
			ll					= {_L1, _L2, _L2, _L3};
			break;

		case 20:
			// RnLnRpLp
			this->_Nseg = 4;
			dir					= {B, B, F, F};			 // back-back-fwd-fwd
			kk					= {R, len, R, len};	 // R len R len
			ll					= {_L1, _L2, _L2, _L3};
			break;

			// C | Cu Cu | C

		case 21:
			// LpRnLnRp
			this->_Nseg = 4;
			dir					= {F, B, B, F};			 // fwd-back-back-fwd
			kk					= {len, R, len, R};	 // len R len R
			ll					= {_L1, _L2, _L2, _L3};
			break;

		case 22:
			// RpLnRnLp
			this->_Nseg = 4;
			dir					= {F, B, B, F};			 // fwd-back-back-fwd
			kk					= {R, len, R, len};	 // R len R len
			ll					= {_L1, _L2, _L2, _L3};
			break;

		case 23:
			// LnRpLpRn
			this->_Nseg = 4;
			dir					= {B, F, F, B};			 // back-fwd-fwd-back
			kk					= {len, R, len, R};	 // len R len R
			ll					= {_L1, _L2, _L2, _L3};
			break;

		case 24:
			// RnLpRpLn
			this->_Nseg = 4;
			dir					= {B, F, F, B};			 // back-fwd-fwd-back
			kk					= {R, len, R, len};	 // R len R len
			ll					= {_L1, _L2, _L2, _L3};
			break;

			// C | C(MPIDIV2) S C
		case 25:
			// LpRnSnLn
			this->_Nseg = 4;
			dir					= {F, B, B, B};			 // fwd-back-back-back
			kk					= {len, R, S, len};	 // len R S len
			ll					= {_L1, MPIDIV2 * lambda, _L2 / lambda, _L3};
			break;

		case 26:
			// RpLnSnRn
			this->_Nseg = 4;
			dir					= {F, B, B, B};		 // fwd-back-back-back
			kk					= {R, len, S, R};	 // R len S R
			ll					= {_L1, MPIDIV2 * lambda, _L2 / lambda, _L3};
			break;

		case 27:
			// LnRpSpLp
			this->_Nseg = 4;
			dir					= {B, F, F, F};			 // back-fwd-fwd-fwd
			kk					= {len, R, S, len};	 // len R S len
			ll					= {_L1, MPIDIV2 * lambda, _L2 / lambda, _L3};
			break;

		case 28:
			// RnLpSpRp
			this->_Nseg = 4;
			dir					= {B, F, F, F};		 // back-fwd-fwd-fwd
			kk					= {R, len, S, R};	 // R len S R
			ll					= {_L1, MPIDIV2 * lambda, _L2 / lambda, _L3};
			break;

		case 29:
			// LpRnSnRn
			this->_Nseg = 4;
			dir					= {F, B, B, B};		 // fwd-back-back-back
			kk					= {len, R, S, R};	 // len R S R
			ll					= {_L1, MPIDIV2 * lambda, _L2 / lambda, _L3};
			break;

		case 30:
			// RpLnSnLn
			this->_Nseg = 4;
			dir					= {F, B, B, B};			 // fwd-back-back-back
			kk					= {R, len, S, len};	 // R len S len
			ll					= {_L1, MPIDIV2 * lambda, _L2 / lambda, _L3};
			break;

		case 31:
			// LnRpSpRp
			this->_Nseg = 4;
			dir					= {B, F, F, F};		 // back-fwd-fwd-fwd
			kk					= {len, R, S, R};	 // len R S R
			ll					= {_L1, MPIDIV2 * lambda, _L2 / lambda, _L3};
			break;

		case 32:
			// RnLpSpLp
			this->_Nseg = 4;
			dir					= {B, F, F, F};			 //  back-fwd-fwd-fwd
			kk					= {R, len, S, len};	 // R len S len
			ll					= {_L1, MPIDIV2 * lambda, _L2 / lambda, _L3};
			break;

			// C |  C(MPIDIV2) S  C(MPIDIV2) | C

		case 33:
			// LpRnSnLnRp
			this->_Nseg = 5;
			dir					= {F, B, B, B, F};			//  fwd-back-back-back-fwd
			kk					= {len, R, S, len, R};	// len R S len R
			ll					= {_L1, MPIDIV2 * lambda, _L2 / lambda, MPIDIV2 * lambda, _L3};
			break;

		case 34:
			// RpLnSnRnLp
			this->_Nseg = 5;
			dir					= {F, B, B, B, F};			//  fwd-back-back-back-fwd
			kk					= {R, len, S, R, len};	// R len S R len
			ll					= {_L1, MPIDIV2 * lambda, _L2 / lambda, MPIDIV2 * lambda, _L3};
			break;

		case 35:
			// LnRpSpLpRn
			this->_Nseg = 5;
			dir					= {B, F, F, F, B};			//  back-fwd-fwd-fwd-back
			kk					= {len, R, S, len, R};	// len R S len R
			ll					= {_L1, MPIDIV2 * lambda, _L2 / lambda, MPIDIV2 * lambda, _L3};
			break;

		case 36:
			// RnLpSpRpLn
			this->_Nseg = 5;
			dir					= {B, F, F, F, B};			//  back-fwd-fwd-fwd-back
			kk					= {R, len, S, R, len};	// R len S R len
			ll					= {_L1, MPIDIV2 * lambda, _L2 / lambda, MPIDIV2 * lambda, _L3};
			break;

			// C C | C
		case 37:
			// LpRpLn
			this->_Nseg = 3;
			dir					= {F, F, B};			//  fwd-fwd-back
			kk					= {len, R, len};	//  len R len
			ll					= {_L1, _L2, _L3};
			break;

		case 38:
			// RpLpRn
			this->_Nseg = 3;
			dir					= {F, F, B};		//  fwd-fwd-back
			kk					= {R, len, R};	// R len R
			ll					= {_L1, _L2, _L3};
			break;

		case 39:
			// LnRnLp
			this->_Nseg = 3;
			dir					= {B, B, F};			//  back-back-fwd
			kk					= {len, R, len};	// len R len
			ll					= {_L1, _L2, _L3};
			break;

		case 40:
			// RnLnRp
			this->_Nseg = 3;
			dir					= {B, B, F};		//  back-back-fwd
			kk					= {R, len, R};	// R len R
			ll					= {_L1, _L2, _L3};
			break;

			// C S C(MPIDIV2) | C

		case 41:
			// LpSpRpLn
			this->_Nseg = 4;
			dir					= {F, F, F, B};			 //  fwd-fwd-fwd-back
			kk					= {len, S, R, len};	 // len S R len
			ll					= {_L1, _L2 / lambda, MPIDIV2 * lambda, _L3};
			break;

		case 42:
			// RpSpLpRn
			this->_Nseg = 4;
			dir					= {F, F, F, B};		 //  fwd-fwd-fwd-back
			kk					= {R, S, len, R};	 // R S len R
			ll					= {_L1, _L2 / lambda, MPIDIV2 * lambda, _L3};
			break;

		case 43:
			// LnSnRnLp
			this->_Nseg = 4;
			dir					= {B, B, B, F};			 //  back-back-back-fwd
			kk					= {len, S, R, len};	 // len S R len
			ll					= {_L1, _L2 / lambda, MPIDIV2 * lambda, _L3};
			break;

		case 44:
			// RnSnLnRp
			this->_Nseg = 4;
			dir					= {B, B, B, F};		 //  back-back-back-fwd
			kk					= {R, S, len, R};	 // R S len R
			ll					= {_L1, _L2 / lambda, MPIDIV2 * lambda, _L3};
			break;

		case 45:
			// LpSpLpRn
			this->_Nseg = 4;
			dir					= {F, F, F, B};			 //  back-back-back-fwd
			kk					= {len, S, len, R};	 // len S len R
			ll					= {_L1, _L2 / lambda, MPIDIV2 * lambda, _L3};
			break;

		case 46:
			// RpSpRpLn
			this->_Nseg = 4;
			dir					= {F, F, F, B};		 //  fwd-fwd-fwd-back
			kk					= {R, S, R, len};	 // R S R len
			ll					= {_L1, _L2 / lambda, MPIDIV2 * lambda, _L3};
			break;

		case 47:
			// LnSnLnRp
			this->_Nseg = 4;
			dir					= {B, B, B, F};			 //  back-back-back-fwd
			kk					= {len, S, len, R};	 // len S len R
			ll					= {_L1, _L2 / lambda, MPIDIV2 * lambda, _L3};
			break;

		case 48:
			// RnSnRnLp
			this->_Nseg = 4;
			dir					= {B, B, B, F};		 //  back-back-back-fwd
			kk					= {R, S, R, len};	 // R S R len
			ll					= {_L1, _L2 / lambda, MPIDIV2 * lambda, _L3};
			break;

		default:
			break;
	}

	c.copy (*this->ci());

	for (seg = 0; seg < this->_Nseg; ++seg)
	{
		d = this->circleLine (ll[seg], dir[seg], kk[seg], _kmax, c, seg);
		c.copy (d);
	}

	//  std::cout << "diff (x,y,th): "
	//            << (std::abs(this->cf()->x()-d.x())) << "  "
	//            << (std::abs(this->cf()->y()-d.y())) << "  "
	//            << (mod2pi(this->cf()->th()-d.th())) << std::endl;

	this->X[seg]	= this->cf()->x();
	this->Y[seg]	= this->cf()->y();
	this->TH[seg] = this->cf()->th();
}

/***********************************************************/
static double
my_atan2 (double y, double x)
{
	double a;
	if ((x == 0.0) && (y == 0.0)) return 0.0;
	if (x == 0.0)
	{
		if (y > 0)
			return MPIDIV2;
		else
			return -MPIDIV2;
	}
	a = atan (y / x);
	if (a > 0.0)
		if (x > 0)
			return a;
		else
			return (a + MPI);
	else if (x > 0)
		return (a + MPIMUL2);
	else
		return (a + MPI);
}

/***********************************************************/
static double
c_c_c (
		double RADCURV,
		double x,
		double y,
		double phi,
		double rs,
		double rc,
		double* t,
		double* u,
		double* v)
{
	double a, b, u1, theta, alpha, length_rs;

	a = x - rs;
	b = y + rc;
	if ((fabs (a) < EPS3) && (fabs (b) < EPS3))
		return (std::numeric_limits<double>::infinity());
	u1 = sqrt (a * a + b * b);
	if (u1 > RADCURVMUL4) return (std::numeric_limits<double>::infinity());
	theta = my_atan2 (b, a);
	alpha = acos (u1 / RADCURVMUL4);
	*t		= mod2pi (MPIDIV2 + alpha + theta);
	*u		= mod2pi (MPI - 2 * alpha);
	*v		= mod2pi (phi - *t - *u);

	length_rs = RADCURV * (*t + *u + *v);
	//  std::cout << "*" << *t*RADCURV << " " << *u*RADCURV << " " << *v*RADCURV <<
	//  std::endl;

	return (length_rs);
}

/***********************************************************/
static double
c_cc (
		double RADCURV,
		double x,
		double y,
		double phi,
		double rs,
		double rc,
		double* t,
		double* u,
		double* v)
{
	double a, b, u1, theta, alpha, length_rs;

	a = x - rs;
	b = y + rc;
	if ((fabs (a) < EPS3) && (fabs (b) < EPS3))
		return (std::numeric_limits<double>::infinity());
	u1 = sqrt (a * a + b * b);
	if (u1 > RADCURVMUL4) return (std::numeric_limits<double>::infinity());
	theta = my_atan2 (b, a);
	alpha = acos (u1 / RADCURVMUL4);
	*t		= mod2pi (MPIDIV2 + alpha + theta);
	*u		= mod2pi (MPI - 2 * alpha);
	*v		= mod2pi (*t + *u - phi);

	length_rs = RADCURV * (*t + *u + *v);
	return (length_rs);
}

/***********************************************************/
static double
csca (
		double RADCURV,
		double x,
		double y,
		double phi,
		double rs,
		double rc,
		double* t,
		double* u,
		double* v)
{
	double a, b, length_rs;

	a	 = x - rs;
	b	 = y + rc;
	*t = mod2pi (my_atan2 (b, a));
	*u = sqrt (a * a + b * b);
	*v = mod2pi (phi - *t);

	length_rs = RADCURV * (*t + *v) + *u;
	return (length_rs);
}

/***********************************************************/
static double
cscb (
		double RADCURV,
		double x,
		double y,
		double phi,
		double rs,
		double rc,
		double* t,
		double* u,
		double* v)
{
	double a, b, u1, theta, alpha, length_rs;

	a	 = x + rs;
	b	 = y - rc;
	u1 = sqrt (a * a + b * b);
	if (u1 < RADCURVMUL2) return (std::numeric_limits<double>::infinity());
	theta = my_atan2 (b, a);
	*u		= sqrt (u1 * u1 - SQRADCURVMUL2);
	alpha = my_atan2 (RADCURVMUL2, *u);
	*t		= mod2pi (theta + alpha);
	*v		= mod2pi (*t - phi);

	length_rs = RADCURV * (*t + *v) + *u;
	return (length_rs);
}

/***********************************************************/
static double
ccu_cuc (
		double RADCURV,
		double x,
		double y,
		double phi,
		double rs,
		double rc,
		double* t,
		double* u,
		double* v)
{
	double a, b, u1, theta, alpha, length_rs;

	a = x + rs;
	b = y - rc;
	if ((fabs (a) < EPS3) && (fabs (b) < EPS3))
		return (std::numeric_limits<double>::infinity());
	u1 = sqrt (a * a + b * b);
	if (u1 > RADCURVMUL4) return (std::numeric_limits<double>::infinity());
	theta = my_atan2 (b, a);
	if (u1 > RADCURVMUL2)
	{
		alpha = acos ((u1 / 2 - RADCURV) / RADCURVMUL2);
		*t		= mod2pi (MPIDIV2 + theta - alpha);
		*u		= mod2pi (MPI - alpha);
		*v		= mod2pi (phi - *t + 2 * (*u));
	}
	else
	{
		alpha = acos ((u1 / 2 + RADCURV) / (RADCURVMUL2));
		*t		= mod2pi (MPIDIV2 + theta + alpha);
		*u		= mod2pi (alpha);
		*v		= mod2pi (phi - *t + 2 * (*u));
	}

	length_rs = RADCURV * (2 * (*u) + *t + *v);
	return (length_rs);
}

/***********************************************************/
static double
c_cucu_c (
		double RADCURV,
		double x,
		double y,
		double phi,
		double rs,
		double rc,
		double* t,
		double* u,
		double* v)
{
	double a, b, u1, theta, alpha, length_rs, va1, va2;

	a = x + rs;
	b = y - rc;
	if ((fabs (a) < EPS3) && (fabs (b) < EPS3))
		return (std::numeric_limits<double>::infinity());
	u1 = sqrt (a * a + b * b);
	if (u1 > 6 * RADCURV) return (std::numeric_limits<double>::infinity());
	theta = my_atan2 (b, a);
	va1		= (5 * SQRADCURV - u1 * u1 / 4) / SQRADCURVMUL2;
	if ((va1 < 0.0) || (va1 > 1.0)) return (std::numeric_limits<double>::infinity());
	*u		= acos (va1);
	va2		= sin (*u);
	alpha = asin (RADCURVMUL2 * va2 / u1);
	*t		= mod2pi (MPIDIV2 + theta + alpha);
	*v		= mod2pi (*t - phi);

	length_rs = RADCURV * (2 * (*u) + *t + *v);
	return (length_rs);
}

/***********************************************************/
static double
c_c2sca (
		double RADCURV,
		double x,
		double y,
		double phi,
		double rs,
		double rc,
		double* t,
		double* u,
		double* v)
{
	double a, b, u1, theta, alpha, length_rs;

	a	 = x - rs;
	b	 = y + rc;
	u1 = sqrt (a * a + b * b);
	if (u1 < RADCURVMUL2) return (std::numeric_limits<double>::infinity());
	theta = my_atan2 (b, a);
	*u		= sqrt (u1 * u1 - SQRADCURVMUL2) - RADCURVMUL2;
	if (*u < 0.0) return (std::numeric_limits<double>::infinity());
	alpha = my_atan2 (RADCURVMUL2, (*u + RADCURVMUL2));
	*t		= mod2pi (MPIDIV2 + theta + alpha);
	*v		= mod2pi (*t + MPIDIV2 - phi);

	length_rs = RADCURV * (*t + MPIDIV2 + *v) + *u;
	return (length_rs);
}

/***********************************************************/
static double
c_c2scb (
		double RADCURV,
		double x,
		double y,
		double phi,
		double rs,
		double rc,
		double* t,
		double* u,
		double* v)
{
	double a, b, u1, theta, length_rs;

	a	 = x + rs;
	b	 = y - rc;
	u1 = sqrt (a * a + b * b);
	if (u1 < RADCURVMUL2) return (std::numeric_limits<double>::infinity());
	theta = my_atan2 (b, a);
	*t		= mod2pi (MPIDIV2 + theta);
	*u		= u1 - RADCURVMUL2;
	*v		= mod2pi (phi - *t - MPIDIV2);

	length_rs = RADCURV * (*t + MPIDIV2 + *v) + *u;
	return (length_rs);
}

/***********************************************************/
static double
c_c2sc2_c (
		double RADCURV,
		double x,
		double y,
		double phi,
		double rs,
		double rc,
		double* t,
		double* u,
		double* v)
{
	double a, b, u1, theta, alpha, length_rs;

	a	 = x + rs;
	b	 = y - rc;
	u1 = sqrt (a * a + b * b);
	if (u1 < RADCURVMUL4) return (std::numeric_limits<double>::infinity());
	theta = my_atan2 (b, a);
	*u		= sqrt (u1 * u1 - SQRADCURVMUL2) - RADCURVMUL4;
	if (*u < 0.0) return (std::numeric_limits<double>::infinity());
	alpha = my_atan2 (RADCURVMUL2, (*u + RADCURVMUL4));
	*t		= mod2pi (MPIDIV2 + theta + alpha);
	*v		= mod2pi (*t - phi);

	length_rs = RADCURV * (*t + MPI + *v) + *u;
	return (length_rs);
}

/***********************************************************/
static double
cc_c (
		double RADCURV,
		double x,
		double y,
		double phi,
		double rs,
		double rc,
		double* t,
		double* u,
		double* v)
{
	double a, b, u1, theta, alpha, length_rs, va;

	a = x - rs;
	b = y + rc;
	if ((fabs (a) < EPS3) && (fabs (b) < EPS3))
		return (std::numeric_limits<double>::infinity());
	u1 = sqrt (a * a + b * b);
	if (u1 > RADCURVMUL4) return (std::numeric_limits<double>::infinity());
	theta = my_atan2 (b, a);
	*u		= acos ((8 * SQRADCURV - u1 * u1) / (8 * SQRADCURV));
	va		= sin (*u);
	if (fabs (va) < 0.001) va = 0.0;
	if ((fabs (va) < 0.001) && (fabs (u1) < 0.001))
		return (std::numeric_limits<double>::infinity());
	alpha = asin (RADCURVMUL2 * va / u1);
	*t		= mod2pi (MPIDIV2 - alpha + theta);
	*v		= mod2pi (*t - *u - phi);

	length_rs = RADCURV * (*t + *u + *v);
	return (length_rs);
}

/***********************************************************/
static double
csc2_ca (
		double RADCURV,
		double x,
		double y,
		double phi,
		double rs,
		double rc,
		double* t,
		double* u,
		double* v)
{
	double a, b, u1, theta, alpha, length_rs;

	a	 = x - rs;
	b	 = y + rc;
	u1 = sqrt (a * a + b * b);
	if (u1 < RADCURVMUL2) return (std::numeric_limits<double>::infinity());
	theta = my_atan2 (b, a);
	*u		= sqrt (u1 * u1 - SQRADCURVMUL2) - RADCURVMUL2;
	if (*u < 0.0) return (std::numeric_limits<double>::infinity());
	alpha = my_atan2 ((*u + RADCURVMUL2), RADCURVMUL2);
	*t		= mod2pi (MPIDIV2 + theta - alpha);
	*v		= mod2pi (*t - MPIDIV2 - phi);

	length_rs = RADCURV * (*t + MPIDIV2 + *v) + *u;
	return (length_rs);
}

/***********************************************************/
static double
csc2_cb (
		double RADCURV,
		double x,
		double y,
		double phi,
		double rs,
		double rc,
		double* t,
		double* u,
		double* v)
{
	double a, b, u1, theta, length_rs;

	a	 = x + rs;
	b	 = y - rc;
	u1 = sqrt (a * a + b * b);
	if (u1 < RADCURVMUL2) return (std::numeric_limits<double>::infinity());
	theta = my_atan2 (b, a);
	*t		= mod2pi (theta);
	*u		= u1 - RADCURVMUL2;
	*v		= mod2pi (-*t - MPIDIV2 + phi);

	length_rs = RADCURV * (*t + MPIDIV2 + *v) + *u;
	return (length_rs);
}

/***********************************************************/

double
RS::reeds_shepp (int Nman, std::vector<double>* debug)
{
	double x1 = this->ci()->x();
	double y1 = this->ci()->y();
	double t1 = this->ci()->th();
	double x2 = this->cf()->x();
	double y2 = this->cf()->y();
	double t2 = this->cf()->th();

	double RADCURV = 1.0 / this->getKmax();

	double x, y, phi;
	double t, u, v, tn, un, vn;
	int num;
	double var, vard, theta, alpha, dx, dy;
	double length = this->_L;
	double sphi, cphi;
	double ap, am, b1, b2;

	/* coordinate change */
	dx		= x2 - x1;
	dy		= y2 - y1;
	theta = my_atan2 (dy, dx);
	alpha = theta - t1;
	vard	= sqrt (dx * dx + dy * dy);
	x			= cos (alpha) * vard;
	y			= sin (alpha) * vard;
	phi		= t2 - t1;

	sphi = sin (phi);
	cphi = cos (phi);

	ap = RADCURV * sphi;
	am = -RADCURV * sphi;
	b1 = RADCURV * (cphi - 1);
	b2 = RADCURV * (cphi + 1);

	// There is no num = 0 in this code, so I need to add this
	SAVE_DEBUG (debug, 0.0)
	//  SAVE_DEBUG(debug, std::numeric_limits<double>::infinity())

	/*   C | C | C   */

	if (Nman == -1 || Nman == 1)
	{
		length = c_c_c (RADCURV, x, y, phi, ap, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, length)
		num = 1;
		t		= tn;
		u		= un;
		v		= vn;
	}

	if (Nman == -1 || Nman == 2)
	{
		var = c_c_c (RADCURV, -x, y, -phi, am, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length && false)
		{
			length = var;
			num		 = 2;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 3)
	{
		var = c_c_c (RADCURV, x, -y, -phi, am, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 3;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 4)
	{
		var = c_c_c (RADCURV, -x, -y, phi, ap, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length && false)
		{
			length = var;
			num		 = 4;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	/*   C | C C   */

	if (Nman == -1 || Nman == 5)
	{
		var = c_cc (RADCURV, x, y, phi, ap, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 5;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 6)
	{
		var = c_cc (RADCURV, -x, y, -phi, am, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 6;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 7)
	{
		var = c_cc (RADCURV, x, -y, -phi, am, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 7;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 8)
	{
		var = c_cc (RADCURV, -x, -y, phi, ap, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 8;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	/*   C S C   */

	if (Nman == -1 || Nman == 9)
	{
		var = csca (RADCURV, x, y, phi, ap, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 9;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 10)
	{
		var = csca (RADCURV, x, -y, -phi, am, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 10;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 11)
	{
		var = csca (RADCURV, -x, y, -phi, am, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 11;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 12)
	{
		var = csca (RADCURV, -x, -y, phi, ap, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 12;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 13)
	{
		var = cscb (RADCURV, x, y, phi, ap, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 13;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 14)
	{
		var = cscb (RADCURV, x, -y, -phi, am, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 14;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 15)
	{
		var = cscb (RADCURV, -x, y, -phi, am, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 15;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 16)
	{
		var = cscb (RADCURV, -x, -y, phi, ap, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 16;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	/*   C Cu | Cu C   */

	if (Nman == -1 || Nman == 17)
	{
		var = ccu_cuc (RADCURV, x, y, phi, ap, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 17;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 18)
	{
		var = ccu_cuc (RADCURV, x, -y, -phi, am, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 18;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 19)
	{
		var = ccu_cuc (RADCURV, -x, y, -phi, am, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 19;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 20)
	{
		var = ccu_cuc (RADCURV, -x, -y, phi, ap, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 20;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	/*   C | Cu Cu | C   */

	if (Nman == -1 || Nman == 21)
	{
		var = c_cucu_c (RADCURV, x, y, phi, ap, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 21;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 22)
	{
		var = c_cucu_c (RADCURV, x, -y, -phi, am, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 22;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 23)
	{
		var = c_cucu_c (RADCURV, -x, y, -phi, am, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 23;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 24)
	{
		var = c_cucu_c (RADCURV, -x, -y, phi, ap, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 24;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	/*   C | C2 S C   */

	if (Nman == -1 || Nman == 25)
	{
		var = c_c2sca (RADCURV, x, y, phi, ap, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 25;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 26)
	{
		var = c_c2sca (RADCURV, x, -y, -phi, am, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 26;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 27)
	{
		var = c_c2sca (RADCURV, -x, y, -phi, am, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 27;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 28)
	{
		var = c_c2sca (RADCURV, -x, -y, phi, ap, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 28;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 29)
	{
		var = c_c2scb (RADCURV, x, y, phi, ap, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 29;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 30)
	{
		var = c_c2scb (RADCURV, x, -y, -phi, am, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 30;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 31)
	{
		var = c_c2scb (RADCURV, -x, y, -phi, am, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 31;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 32)
	{
		var = c_c2scb (RADCURV, -x, -y, phi, ap, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 32;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	/*   C | C2 S C2 | C   */

	if (Nman == -1 || Nman == 33)
	{
		var = c_c2sc2_c (RADCURV, x, y, phi, ap, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 33;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 34)
	{
		var = c_c2sc2_c (RADCURV, x, -y, -phi, am, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 34;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 35)
	{
		var = c_c2sc2_c (RADCURV, -x, y, -phi, am, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 35;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 36)
	{
		var = c_c2sc2_c (RADCURV, -x, -y, phi, ap, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 36;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	/*   C C | C   */

	if (Nman == -1 || Nman == 37)
	{
		var = cc_c (RADCURV, x, y, phi, ap, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 37;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 38)
	{
		var = cc_c (RADCURV, x, -y, -phi, am, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 38;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 39)
	{
		var = cc_c (RADCURV, -x, y, -phi, am, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 39;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 40)
	{
		var = cc_c (RADCURV, -x, -y, phi, ap, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 40;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	/*   C S C2 | C   */

	if (Nman == -1 || Nman == 41)
	{
		var = csc2_ca (RADCURV, x, y, phi, ap, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 41;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 42)
	{
		var = csc2_ca (RADCURV, x, -y, -phi, am, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 42;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 43)
	{
		var = csc2_ca (RADCURV, -x, y, -phi, am, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 43;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 44)
	{
		var = csc2_ca (RADCURV, -x, -y, phi, ap, b1, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 44;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 45)
	{
		var = csc2_cb (RADCURV, x, y, phi, ap, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 45;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 46)
	{
		var = csc2_cb (RADCURV, x, -y, -phi, am, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 46;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 47)
	{
		var = csc2_cb (RADCURV, -x, y, -phi, am, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 47;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	if (Nman == -1 || Nman == 48)
	{
		var = csc2_cb (RADCURV, -x, -y, phi, ap, b2, &tn, &un, &vn);
		SAVE_DEBUG (debug, var)
		if (var < length)
		{
			length = var;
			num		 = 48;
			t			 = tn;
			u			 = un;
			v			 = vn;
		}
	}

	this->_L1		= t * RADCURV;
	this->_L2		= u * RADCURV;
	this->_L3		= v * RADCURV;
	this->_Nman = num;
	this->_L		= length;


	return length;
}

std::vector<RSSegment>
RS::getSegmentsData()
{
	std::vector<RSSegment> ret (this->getNseg());
	for (size_t i = 0; i < this->getNseg(); i++)
	{
		ret[i] = RSSegment (
				this->X[i], this->Y[i], this->TH[i], this->TH[i + 1], this->L[i], this->K[i],
				(int)this->D[i]);
	}
	return ret;
}

std::vector<std::vector<double>>
RS::split_wise()
{
	return {};
}

static void
drawAngleArrow (
		std::ofstream& file,
		const Configuration2& c,
		double length,
		const std::string& pen)
{
	const double x1 = c.x() + length * std::cos (c.th());
	const double y1 = c.y() + length * std::sin (c.th());
	file << "draw((" << c.x() << "," << c.y() << ")--(" << x1 << "," << y1 << "), "
			 << pen << ", Arrow);" << std::endl;
}

[[nodiscard]]
static
Configuration2
circleLine_helper(double s, double dir, double kur, double kmax, Configuration2 c)
{
	double sigmaDir, sign;
	if (dir > 0)
	{
		sigmaDir = 0;
		sign	 = 1;
	}
	else
	{
		sigmaDir = 1;
		sign	 = -1;
	}
	double xEnd, yEnd, thetaEnd;
	xEnd     = c.x() + f (s, sign * kur * kmax, mod2pi (c.th() + sigmaDir * MPI));
	yEnd	 = c.y() + g (s, sign * kur * kmax, mod2pi (c.th() + sigmaDir * MPI));
	thetaEnd = mod2pi (c.th() + sign * kur * kmax * s);
	return Configuration2(xEnd, yEnd, thetaEnd, kur * kmax);
}

// #ifdef MPDP_DRAW
void
RS::draw (
		std::ofstream& file,
		size_t width,
		size_t height,
		bool solve,
		bool close,
		bool init,
		bool axes,
		std::pair<std::string, std::string> arrows_pen,
		std::pair<std::string, std::string> points_pen, 
		std::vector<std::string> segments_pen
){
	if (solve) { this->solve(); }

	if (init) { initAsyFile (file); }

	Configuration2 c0 (this->ci()->x(), this->ci()->y(), this->ci()->th());
	double xmin = this->ci()->x();
	double xmax = this->ci()->x();
	double ymin = this->ci()->y();
	double ymax = this->ci()->y();
	for (size_t i = 0; i <= this->getNseg(); i++)
	{
		xmin = std::min (xmin, this->X[i]);
		xmax = std::max (xmax, this->X[i]);
		ymin = std::min (ymin, this->Y[i]);
		ymax = std::max (ymax, this->Y[i]);
	}
	xmin = std::min (xmin, this->cf()->x());
	xmax = std::max (xmax, this->cf()->x());
	ymin = std::min (ymin, this->cf()->y());
	ymax = std::max (ymax, this->cf()->y());
	const double arrowLength = std::max (0.15, 0.12 * std::hypot (xmax - xmin, ymax - ymin));

	for (size_t i = 0; i < this->getNseg(); i++)
	{
		const double dir = this->D[i];
		const double drawK = this->K[i];
		const double drawL = this->L[i];

		const double x0 = std::abs(c0.x()) < 1e-12? 0.0 : c0.x();
		const double y0 = std::abs(c0.y()) < 1e-12? 0.0 : c0.y();
		const double theta0 = std::abs(c0.th()) < 1e-12? 0.0 : c0.th();
		const double k0 = std::abs(drawK) < 1e-12? 0.0 : drawK;
		const double l0 = std::abs(drawL) < 1e-12? 0.0 : drawL;

		file << "p = clothoidPoints((" << x0 << "," << y0 << "), " << theta0 << ","
				 << k0 << ", 0, " << l0 << ");" << std::endl;
		file << "draw(p," << segments_pen[0] << ");" << std::endl;
		file << "dot((" << x0 << "," << y0 << "), " << (i == 0 ? points_pen.first : points_pen.second) << ");" << std::endl;

		double kur = this->K[i] > 0 ? 1 : -1;
		kur = std::abs(this->K[i]) < 1e-12? 0.0 : kur;

		c0 = circleLine_helper(std::abs(drawL), dir, kur, std::abs(k0), c0);
	}
	// Final point
	drawAngleArrow (file, *this->ci(), arrowLength, arrows_pen.first);
	file << "dot((" << this->cf()->x() << "," << this->cf()->y() << "), " << points_pen.first << ");" << std::endl;
	drawAngleArrow (file, *this->cf(), arrowLength, arrows_pen.second);

	if (axes)
	{
		file << "xaxis(\"$x$\", BottomTop(), Ticks(Label(\"$%.2f$\")));" << std::endl;
		file << "yaxis(\"$y$\", LeftRight(), Ticks(Label(\"$%.2f$\")));" << std::endl;
	}

	if (close) { file.close(); }
}
// #endif	// MPDP_DRAW

#endif	// CUDA_ON





// Warning: Optimal maneuver not found in top-4 predictions. Expected: 27, Optimal length: 2.31905, Best predicted length: 2.31911
// thi: 0.884363, thf: -0.697575, kmax: 0.740491
// Top-4 predictions (maneuver, logit): (6, 10.1079) (3, 9.19569) (31, 9.03442) (10, 7.57637) 
// Optimal maneuver: 27, Optimal length: 2.31905
// Best predicted length: 2.31911 with maneuver 6
// All class probabilities (maneuver, probability, descending): (6, 5.3619294356518121e-01) (3, 2.1535634862777936e-01) (31, 1.8328250732422996e-01) (10, 4.2648081138540947e-02) (27, 1.8755843284294332e-02) (14, 3.7642757053288306e-03) (39, 3.5141501191166978e-10) (47, 1.8660359066963506e-12) (32, 1.3163931344810862e-12) (35, 4.5916714090304227e-14) (23, 1.8433295139750520e-15) (20, 1.4207454671329626e-16) (11, 1.4762957453124174e-17) (43, 9.0321375007597795e-20) (15, 1.6110840332767411e-22) (1, 1.9811119141286992e-33) (21, 1.6954529747793814e-33) (40, 2.1095869520311942e-44) (16, 4.4316906214644634e-47) (33, 1.2320765279369168e-57) (12, 1.1908462820562216e-61) (48, 1.6652289529724813e-80) (44, 4.6343843798379190e-102) 

// Warning: Optimal maneuver not found in top-4 predictions. Expected: 23, Optimal length: 3.5327, Best predicted length: 3.53299
// thi: 2.08025, thf: 1.06108, kmax: 0.926986
// Top-4 predictions (maneuver, logit): (20, 14.4419) (40, 13.8888) (21, 13.6641) (6, 13.3097) 
// Optimal maneuver: 23, Optimal length: 3.5327
// Best predicted length: 3.53299 with maneuver 21
// All class probabilities (maneuver, probability, descending): (20, 4.0503010357640551e-01) (40, 2.3295883447864194e-01) (21, 1.8607083891045367e-01) (6, 1.3054739717178004e-01) (23, 4.5392744111073614e-02) (33, 4.8330713645591912e-08) (35, 2.9144561209584800e-08) (3, 4.0659546634363764e-09) (44, 1.6502811249313744e-10) (27, 3.8428044501071808e-11) (47, 5.5049961551523132e-12) (32, 1.4544018384861334e-12) (1, 5.7111670774278195e-17) (11, 2.9350493328271763e-21) (39, 5.8425494640601974e-27) (43, 1.9317192444108328e-32) (15, 4.7035867882246467e-33) (48, 1.8022686632630279e-35) (14, 2.0990115359255723e-42) (31, 5.6758975053713648e-50) (16, 1.1330568318182306e-51) (12, 3.2340544055154197e-57) (10, 3.0716746655285624e-87) 

// Warning: Optimal maneuver not found in top-4 predictions. Expected: 31, Optimal length: 2.85252, Best predicted length: 2.85252
// thi: 3.12478, thf: -0.00322971, kmax: 1.33908
// Top-4 predictions (maneuver, logit): (27, -4.38877) (43, -5.38799) (48, -5.75033) (44, -6.8964) 
// Optimal maneuver: 31, Optimal length: 2.85252
// Best predicted length: 2.85252 with maneuver 48
// All class probabilities (maneuver, probability, descending): (27, 5.2354576776202943e-01) (43, 1.9275304380188080e-01) (48, 1.3416473355580852e-01) (44, 4.2648483676134108e-02) (31, 3.7415520895385315e-02) (6, 3.1417957782083943e-02) (40, 1.2575958841334592e-02) (39, 1.2410599638793238e-02) (32, 1.0403855739534024e-02) (47, 2.2296096410607301e-03) (1, 4.1230456817357151e-04) (3, 2.2146354213597680e-05) (11, 1.7022084841881735e-08) (14, 7.2136106662348854e-10) (20, 1.2135949982814340e-13) (35, 3.4131755981763137e-16) (33, 2.0150857788549851e-16) (23, 1.5818158302110591e-16) (15, 1.0820057766601036e-16) (21, 1.6704075989013887e-21) (10, 1.3933881811904772e-24) (16, 1.0585933590789597e-35) (12, 1.6788322850744358e-60) 

// Warning: Optimal maneuver not found in top-4 predictions. Expected: 48, Optimal length: 3.05283, Best predicted length: 3.05996
// thi: 1.57809, thf: -1.56737, kmax: 1.08065
// Top-4 predictions (maneuver, logit): (3, -8.99394) (1, -10.4754) (31, -11.714) (47, -13.6928) 
// Optimal maneuver: 48, Optimal length: 3.05283
// Best predicted length: 3.05996 with maneuver 31
// All class probabilities (maneuver, probability, descending): (3, 7.6002253298981193e-01) (1, 1.7276517225221674e-01) (31, 5.0063133151920590e-02) (47, 6.9202154605260318e-03) (32, 4.2578422174165478e-03) (48, 3.1043563353754123e-03) (10, 1.6363660943711231e-03) (11, 7.1110164040218929e-04) (12, 4.6050523154518618e-04) (39, 4.1187620526450723e-05) (15, 1.6148296636015453e-05) (43, 1.2726535414656979e-06) (35, 9.3631075548228437e-08) (21, 3.0482675480524684e-08) (6, 2.5449901698195846e-08) (33, 1.5105684090975357e-08) (14, 1.3861060458301268e-09) (20, 2.3993243349938460e-13) (23, 2.3158189498395564e-14) (27, 2.6423972377440189e-15) (40, 1.2746202756207469e-15) (16, 4.0802611184607036e-16) (44, 1.3528012208550621e-36) 

// Warning: Optimal maneuver not found in top-4 predictions. Expected: 14, Optimal length: 2.27049, Best predicted length: 9.88002
// thi: 0.909844, thf: -0.604671, kmax: 0.825669
// Top-4 predictions (maneuver, logit): (10, 10.6682) (31, 10.177) (27, 9.78652) (6, 9.63331) 
// Optimal maneuver: 14, Optimal length: 2.27049
// Best predicted length: 9.88002 with maneuver 27
// All class probabilities (maneuver, probability, descending): (10, 3.8559422111579317e-01) (31, 2.3594066653695986e-01) (27, 1.5967574124155606e-01) (6, 1.3699350939104610e-01) (14, 8.1500274170944839e-02) (3, 2.9558752229446217e-04) (39, 1.9276515050154631e-11) (47, 1.1756686412290057e-12) (32, 9.4930296535904507e-13) (35, 3.8723023638850585e-15) (23, 7.3148398335300008e-17) (11, 1.2476368836317330e-17) (20, 1.0886768621652694e-18) (43, 6.2920802205908075e-20) (15, 2.8352579500042700e-23) (21, 2.5215433617824167e-34) (1, 3.6623320698286147e-36) (40, 6.4541235104224583e-45) (16, 2.6246848392625383e-47) (33, 1.6555418497136959e-58) (12, 1.3454607105021156e-62) (48, 3.9934655360299138e-81) (44, 1.2605003182865163e-100) 

// Warning: Optimal maneuver not found in top-4 predictions. Expected: 27, Optimal length: 2.27131, Best predicted length: 2.27132
// thi: 0.831105, thf: -0.746854, kmax: 0.719206
// Top-4 predictions (maneuver, logit): (10, 9.78096) (3, 9.69231) (6, 9.18104) (31, 7.52883) 
// Optimal maneuver: 27, Optimal length: 2.27131
// Best predicted length: 2.27132 with maneuver 6
// All class probabilities (maneuver, probability, descending): (10, 3.8756084307540678e-01) (3, 3.5468049816348912e-01) (6, 2.1271343883233590e-01) (31, 4.0761521166371738e-02) (14, 2.7529039603049884e-03) (27, 1.5307941047010532e-03) (39, 6.9475674414717758e-10) (47, 1.4540201952094990e-12) (32, 1.0282580983730188e-12) (35, 1.4782640678996131e-13) (23, 3.4122289034286873e-15) (20, 1.5666147817470403e-16) (11, 2.0991499627812628e-17) (43, 3.0453071855597101e-19) (15, 2.1531020770790116e-21) (1, 8.4870891111570852e-34) (21, 3.1200824730707383e-34) (40, 2.3498987457422572e-45) (16, 6.8426364451522154e-47) (33, 6.3154184293087292e-58) (12, 1.8317084472514705e-61) (48, 2.1224712757452167e-82) (44, 2.3893013554769062e-104) 





// Warning: Optimal maneuver not found in top-5 predictions. Expected: 48, Optimal length: 3.05283, Best predicted length: 3.05996
// thi: 1.57809, thf: -1.56737, kmax: 1.08065
// Top-5 predictions (maneuver, logit): (3, -8.99394) (1, -10.4754) (31, -11.714) (47, -13.6928) (32, -14.1785) 
// Optimal maneuver: 48, Optimal length: 3.05283
// Best predicted length: 3.05996 with maneuver 31
// All class probabilities (maneuver, probability, descending): (3, 7.6002253298981193e-01) (1, 1.7276517225221674e-01) (31, 5.0063133151920590e-02) (47, 6.9202154605260318e-03) (32, 4.2578422174165478e-03) (48, 3.1043563353754123e-03) (10, 1.6363660943711231e-03) (11, 7.1110164040218929e-04) (12, 4.6050523154518618e-04) (39, 4.1187620526450723e-05) (15, 1.6148296636015453e-05) (43, 1.2726535414656979e-06) (35, 9.3631075548228437e-08) (21, 3.0482675480524684e-08) (6, 2.5449901698195846e-08) (33, 1.5105684090975357e-08) (14, 1.3861060458301268e-09) (20, 2.3993243349938460e-13) (23, 2.3158189498395564e-14) (27, 2.6423972377440189e-15) (40, 1.2746202756207469e-15) (16, 4.0802611184607036e-16) (44, 1.3528012208550621e-36) 

// Warning: Optimal maneuver not found in top-5 predictions. Expected: 27, Optimal length: 2.27131, Best predicted length: 2.27132
// thi: 0.831105, thf: -0.746854, kmax: 0.719206
// Top-5 predictions (maneuver, logit): (10, 9.78096) (3, 9.69231) (6, 9.18104) (31, 7.52883) (14, 4.83375) 
// Optimal maneuver: 27, Optimal length: 2.27131
// Best predicted length: 2.27132 with maneuver 6
// All class probabilities (maneuver, probability, descending): (10, 3.8756084307540678e-01) (3, 3.5468049816348912e-01) (6, 2.1271343883233590e-01) (31, 4.0761521166371738e-02) (14, 2.7529039603049884e-03) (27, 1.5307941047010532e-03) (39, 6.9475674414717758e-10) (47, 1.4540201952094990e-12) (32, 1.0282580983730188e-12) (35, 1.4782640678996131e-13) (23, 3.4122289034286873e-15) (20, 1.5666147817470403e-16) (11, 2.0991499627812628e-17) (43, 3.0453071855597101e-19) (15, 2.1531020770790116e-21) (1, 8.4870891111570852e-34) (21, 3.1200824730707383e-34) (40, 2.3498987457422572e-45) (16, 6.8426364451522154e-47) (33, 6.3154184293087292e-58) (12, 1.8317084472514705e-61) (48, 2.1224712757452167e-82) (44, 2.3893013554769062e-104) 