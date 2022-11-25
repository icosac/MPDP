#ifndef CUDA_ON

#include <rs.hh>

#define RADCURVMUL2 (2 * RADCURV)
#define RADCURVMUL4 (4 * RADCURV)
#define SQRADCURV (RADCURV * RADCURV)
#define SQRADCURVMUL2 (4 * RADCURV * RADCURV)

const double EPS1 = 1.0e-14;
const double EPS3 = 1.0e-14;
const double EPS4 = 1.0e-14;

const double MPI = 3.1415926535897932385;
const double MPIMUL2 = 6.2831853071795864770;
const double MPIDIV2 = 1.5707963267948966192;

Configuration2 RS::circleLine(double s, double dir, double kur, double kmax, Configuration2 c, int seg) {
  double sigmaDir, sign;
  if (dir > 0) {
    sigmaDir = 0;
    sign = 1;
  }
  else {
    sigmaDir = 1;
    sign = -1;
  }
  double xEnd, yEnd, thetaEnd;
  xEnd = c.x() + f(s, sign*kur*kmax, mod2pi(c.th() + sigmaDir * M_PI));
  yEnd = c.y() + g(s, sign*kur*kmax, mod2pi(c.th() + sigmaDir * M_PI));
  thetaEnd = mod2pi(c.th() + sign* kur * kmax * s);
  if (seg == 0) {
    X[seg] = c.x();
    Y[seg] = c.y();
    TH[seg] = c.th();

  }
  this->X[seg+1] = xEnd;
  this->Y[seg+1] = yEnd;
  this->TH[seg+1] = thetaEnd;
  // this->L[seg] = s;
  this->L[seg] = sign*s;
  this->D[seg] = sign;
  this->K[seg] =  kur * kmax;
  Configuration2 res(xEnd,yEnd,thetaEnd, kur * kmax);
  if (kur > 0) {
    if (dir > 0) {
      this->ManType[seg] = "Lp";
    }
    else {
      this->ManType[seg] = "Ln";
    }
  }
  else if (kur < 0) {
    if (dir > 0) {
      this->ManType[seg] = "Rp";
    }
    else {
      this->ManType[seg] = "Rn";
    }
  }
  else {
    if (dir > 0) {
      this->ManType[seg] = "Sp";
    }
    else {
      this->ManType[seg] = "Sn";
    }

  }

  return res;
}

void RS::buildRS(int man) {
  // if (man == -1) { man = this->_Nman; }
  std::vector<double> dir;
  std::vector<double> kk;
  std::vector<double> ll = { _L1, _L2, _L3 };
  Configuration2 c, d;
  int seg = 0;
  double len = 1;
  double R = -1;
  double S = 0;
  double F = 1;
  double B = -1;
  double lambda = 1/_kmax;

  switch (man)
  {
    // C | C | C
    case 1:
      // LpRnLp
      this->_Nseg = 3;
      dir = { F,B,F }; // fwd-back-fwd
      kk = {len, R , len}; // len R len
      ll = { _L1, _L2, _L3 };
//      std::cout << "-" << _L1 << " " << _L2 << " " << _L3 << std::endl;
//      std::cout << "-" << _L1/lambda << " " << _L2/lambda << " " << _L3/lambda << std::endl;
      break;


    case 2:
      // LnRpLn
      this->_Nseg = 3;
      dir = { B,F,B }; // back-fwd-back
      kk = {len, R, len}; // len R len
      ll = { _L1, _L2, _L3 };
      break;

    case 3:
      // RpLnRp
      this->_Nseg = 3;
      dir = { F, B, F }; // fwd-back-fwd
      kk = {R, len, R }; // R len R
      ll = { _L1, _L2, _L3 };
      break;

    case 4:
      // RnLpRn
      this->_Nseg = 3;
      dir = { B,F,B }; // back-fwd-back
      kk = {R, len, R }; // R len R
      ll = { _L1, _L2, _L3 };
      break;

      // C | C C
    case 5:
      // LpRnLn
      this->_Nseg = 3;
      dir = { F,B,B}; // fwd-back-back
      kk = {len, R, len }; // len R len
      ll = { _L1, _L2, _L3 };
      break;

    case 6:
      // LnRpLp
      this->_Nseg = 3;
      dir = { B,F,F }; // back-fwd-fwd
      kk = {len, R, len }; // len R len
      ll = { _L1, _L2, _L3 };
      break;

    case 7:
      // RpLnRn
      this->_Nseg = 3;
      dir = { F,B,B }; // fwd-back-back
      kk = {R, len, R }; //  R len R
      ll = { _L1, _L2, _L3 };
      break;

    case 8:
      // RnLpRp
      this->_Nseg = 3;
      dir = {B,F,F }; // back-fwd-fwd
      kk = {R, len, R }; //  R len R
      ll = { _L1, _L2, _L3 };
      break;

      // C S C
    case 9:
      // LpSpLp
      this->_Nseg = 3;
      dir = {F,F,F }; // fwd-fwd-fwd
      kk = {len, S, len }; //  len S len
      ll = { _L1, _L2/lambda, _L3 };
      break;

    case 10:
      // RpSpRp
      this->_Nseg = 3;
      dir = { F,F,F }; // fwd-fwd-fwd
      kk = { R,S,R }; //  R S R
      ll = { _L1, _L2/lambda, _L3 };
      break;

    case 11:
      // LnSnLn
      this->_Nseg = 3;
      dir = { B,B,B }; // back-back-back
      kk = {len, S, len }; //  len S len
      ll = { _L1, _L2/lambda, _L3 };
      break;

    case 12:
      // RnSnRn
      this->_Nseg = 3;
      dir = { B,B,B }; // back-back-back
      kk = { R,S,R }; //  R S R
      ll = { _L1, _L2/lambda, _L3 };
      break;

    case 13:
      // LpSpRp
      this->_Nseg = 3;
      dir = { F,F,F }; // fwd-fwd-fwd
      kk = {len, S, R }; //  len S R
      ll = { _L1, _L2/lambda, _L3 };
      break;

    case 14:
      // RpSpLp
      this->_Nseg = 3;
      dir = { F,F,F }; // fwd-fwd-fwd
      kk = {R, S, len }; //  R S len
      ll = { _L1, _L2/lambda, _L3 };
      break;

    case 15:
      // LnSnRn
      this->_Nseg = 3;
      dir = { B,B,B }; // back-back-back
      kk = {len, S, R }; //  len S R
      ll = { _L1, _L2/lambda, _L3 };
      break;

    case 16:
      // RnSnLn
      this->_Nseg = 3;
      dir = {B,B,B }; // back-back-back
      kk = {R, S, len }; //  R S len
      ll = { _L1, _L2/lambda, _L3 };
      break;

      // C C= | C= C


    case 17:
      // LpRpLnRn
      this->_Nseg = 4;
      dir = { F,F,B,B }; // fwd-fwd-back-back
      kk = {len, R, len, R }; // len R len R
      ll = { _L1, _L2, _L2, _L3 };
      break;

    case 18:
      // RpLpRnLn
      this->_Nseg = 4;
      dir = { F,F,B,B}; // fwd-fwd-back-back
      kk = {R, len, R, len }; // R len R len
      ll = { _L1, _L2, _L2, _L3 };
      break;

    case 19:
      // LnRnLpRp
      this->_Nseg = 4;
      dir = { B,B,F,F }; // back-back-fwd-fwd
      kk = {len, R, len, R }; // len R len R
      ll = { _L1, _L2, _L2, _L3 };
      break;

    case 20:
      // RnLnRpLp
      this->_Nseg = 4;
      dir = { B,B,F,F }; // back-back-fwd-fwd
      kk = {R, len, R, len }; // R len R len
      ll = { _L1, _L2, _L2, _L3 };
      break;

      // C | Cu Cu | C

    case 21:
      // LpRnLnRp
      this->_Nseg = 4;
      dir = { F,B,B,F }; // fwd-back-back-fwd
      kk = {len, R, len, R }; // len R len R
      ll = { _L1, _L2, _L2, _L3 };
      break;

    case 22:
      // RpLnRnLp
      this->_Nseg = 4;
      dir = {F,B,B,F }; // fwd-back-back-fwd
      kk = {R, len, R, len }; // R len R len
      ll = { _L1, _L2, _L2, _L3 };
      break;

    case 23:
      // LnRpLpRn
      this->_Nseg = 4;
      dir = {B,F,F,B }; // back-fwd-fwd-back
      kk = {len, R, len, R }; // len R len R
      ll = { _L1, _L2, _L2, _L3 };
      break;

    case 24:
      // RnLpRpLn
      this->_Nseg = 4;
      dir = { B,F,F,B }; // back-fwd-fwd-back
      kk = {R, len, R, len }; // R len R len
      ll = { _L1, _L2, _L2, _L3 };
      break;

      // C | C(M_PI_2) S C
    case 25:
      // LpRnSnLn
      this->_Nseg = 4;
      dir = { F,B,B,B}; // fwd-back-back-back
      kk = {len, R, S, len }; // len R S len
      ll = { _L1, M_PI_2 * lambda, _L2/lambda, _L3 };
      break;

    case 26:
      // RpLnSnRn
      this->_Nseg = 4;
      dir = {F,B,B,B }; // fwd-back-back-back
      kk = {R, len, S, R }; // R len S R
      ll = { _L1, M_PI_2 * lambda, _L2/lambda, _L3 };
      break;

    case 27:
      // LnRpSpLp
      this->_Nseg = 4;
      dir = { B,F,F,F }; // back-fwd-fwd-fwd
      kk = {len, R, S, len }; // len R S len
      ll = { _L1, M_PI_2 * lambda, _L2/lambda, _L3 };
      break;

    case 28:
      // RnLpSpRp
      this->_Nseg = 4;
      dir = { B,F,F,F }; // back-fwd-fwd-fwd
      kk = {R, len, S, R }; // R len S R
      ll = { _L1, M_PI_2 * lambda, _L2/lambda, _L3 };
      break;

    case 29:
      // LpRnSnRn
      this->_Nseg = 4;
      dir = { F,B,B,B }; // fwd-back-back-back
      kk = {len, R, S, R }; // len R S R
      ll = { _L1, M_PI_2 * lambda, _L2/lambda, _L3 };
      break;

    case 30:
      // RpLnSnLn
      this->_Nseg = 4;
      dir = {F,B,B,B}; // fwd-back-back-back
      kk = {R, len, S, len }; // R len S len
      ll = { _L1, M_PI_2 * lambda, _L2/lambda, _L3 };
      break;

    case 31:
      // LnRpSpRp
      this->_Nseg = 4;
      dir = { B,F,F,F }; // back-fwd-fwd-fwd
      kk = {len, R, S, R }; // len R S R
      ll = { _L1, M_PI_2 * lambda, _L2/lambda, _L3 };
      break;

    case 32:
      // RnLpSpLp
      this->_Nseg = 4;
      dir = { B,F,F,F}; //  back-fwd-fwd-fwd
      kk = {R, len, S, len }; // R len S len
      ll = { _L1, M_PI_2 * lambda, _L2/lambda, _L3 };
      break;


      // C |  C(M_PI_2) S  C(M_PI_2) | C

    case 33:
      // LpRnSnLnRp
      this->_Nseg = 5;
      dir = { F,B,B,B,F }; //  fwd-back-back-back-fwd
      kk = {len, R, S, len, R }; // len R S len R
      ll = { _L1, M_PI_2 * lambda, _L2 / lambda, M_PI_2 * lambda, _L3 };
      break;

    case 34:
      // RpLnSnRnLp
      this->_Nseg = 5;
      dir = { F,B,B,B,F }; //  fwd-back-back-back-fwd
      kk = {R, len, S, R, len }; // R len S R len
      ll = { _L1, M_PI_2 * lambda, _L2 / lambda, M_PI_2 * lambda, _L3 };
      break;

    case 35:
      // LnRpSpLpRn
      this->_Nseg = 5;
      dir = {B,F,F,F,B }; //  back-fwd-fwd-fwd-back
      kk = {len, R, S, len, R }; // len R S len R
      ll = { _L1, M_PI_2 * lambda, _L2 / lambda, M_PI_2 * lambda, _L3 };
      break;

    case 36:
      // RnLpSpRpLn
      this->_Nseg = 5;
      dir = { B,F,F,F,B }; //  back-fwd-fwd-fwd-back
      kk = {R, len, S, R, len }; // R len S R len
      ll = { _L1, M_PI_2 * lambda, _L2 / lambda, M_PI_2 * lambda, _L3 };
      break;


      // C C | C
    case 37:
      // LpRpLn
      this->_Nseg = 3;
      dir = { F,F,B }; //  fwd-fwd-back
      kk = {len, R, len }; //  len R len
      ll = { _L1, _L2, _L3 };
      break;

    case 38:
      // RpLpRn
      this->_Nseg = 3;
      dir = { F,F,B }; //  fwd-fwd-back
      kk = {R, len, R }; // R len R
      ll = { _L1, _L2, _L3 };
      break;

    case 39:
      // LnRnLp
      this->_Nseg = 3;
      dir = { B,B,F }; //  back-back-fwd
      kk = {len, R, len }; // len R len
      ll = { _L1, _L2, _L3 };
      break;

    case 40:
      // RnLnRp
      this->_Nseg = 3;
      dir = { B,B,F }; //  back-back-fwd
      kk = {R, len, R }; // R len R
      ll = { _L1, _L2, _L3 };
      break;


      // C S C(M_PI_2) | C

    case 41:
      // LpSpRpLn
      this->_Nseg = 4;
      dir = {F,F,F,B }; //  fwd-fwd-fwd-back
      kk = {len, S, R, len }; // len S R len
      ll = { _L1, _L2 / lambda, M_PI_2 * lambda, _L3 };
      break;

    case 42:
      // RpSpLpRn
      this->_Nseg = 4;
      dir = { F,F,F,B }; //  fwd-fwd-fwd-back
      kk = {R, S, len, R }; // R S len R
      ll = { _L1, _L2/lambda, M_PI_2 * lambda, _L3 };
      break;

    case 43:
      // LnSnRnLp
      this->_Nseg = 4;
      dir = { B,B,B,F }; //  back-back-back-fwd
      kk = {len, S, R, len }; // len S R len
      ll = { _L1, _L2/lambda, M_PI_2 * lambda, _L3 };
      break;

    case 44:
      // RnSnLnRp
      this->_Nseg = 4;
      dir = { B,B,B,F }; //  back-back-back-fwd
      kk = {R, S, len, R }; // R S len R
      ll = { _L1, _L2/lambda, M_PI_2 * lambda, _L3 };
      break;

    case 45:
      // LpSpLpRn
      this->_Nseg = 4;
      dir = { F,F,F,B }; //  back-back-back-fwd
      kk = {len, S, len, R }; // len S len R
      ll = { _L1, _L2/lambda, M_PI_2 * lambda, _L3 };
      break;

    case 46:
      // RpSpRpLn
      this->_Nseg = 4;
      dir = { F,F,F,B }; //  fwd-fwd-fwd-back
      kk = {R, S, R, len }; // R S R len
      ll = { _L1, _L2/lambda, M_PI_2 * lambda, _L3 };
      break;

    case 47:
      // LnSnLnRp
      this->_Nseg = 4;
      dir = { B,B,B,F }; //  back-back-back-fwd
      kk = {len, S, len, R }; // len S len R
      ll = { _L1, _L2/lambda, M_PI_2 * lambda, _L3 };
      break;

    case 48:
      // RnSnRnLp
      this->_Nseg = 4;
      dir = { B,B,B,F }; //  back-back-back-fwd
      kk = {R, S, R, len }; // R S R len
      ll = { _L1, _L2/lambda, M_PI_2 * lambda, _L3 };
      break;

    default:
      break;
  }

  c.copy(*this->ci());

  for (seg = 0; seg < this->_Nseg; ++seg) {
    d = this->circleLine(ll[seg], dir[seg], kk[seg], _kmax, c, seg);
    c.copy(d);
  }

//  std::cout << "diff (x,y,th): "
//            << (std::abs(this->cf()->x()-d.x())) << "  "
//            << (std::abs(this->cf()->y()-d.y())) << "  "
//            << (mod2pi(this->cf()->th()-d.th())) << std::endl;

  this->X[seg] = this->cf()->x();
  this->Y[seg] = this->cf()->y();
  this->TH[seg] = this->cf()->th();

}

/***********************************************************/
static double my_atan2(double y, double x) {
  double a;
  if ((x == 0.0) && (y == 0.0))
    return 0.0;
  if (x == 0.0) {
    if (y > 0)
      return MPIDIV2;
    else
      return -MPIDIV2;
  }
  a = atan(y / x);
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
static double c_c_c(double RADCURV,
                    double x,
                    double y,
                    double phi,
                    double rs,
                    double rc,
                    double* t,
                    double* u,
                    double* v) {
  double a, b, u1, theta, alpha, length_rs;

  a = x - rs;
  b = y + rc;
  if ((fabs(a) < EPS3) && (fabs(b) < EPS3))
    return (std::numeric_limits<double>::infinity());
  u1 = sqrt(a * a + b * b);
  if (u1 > RADCURVMUL4)
    return (std::numeric_limits<double>::infinity());
  theta = my_atan2(b, a);
  alpha = acos(u1 / RADCURVMUL4);
  *t = mod2pi(MPIDIV2 + alpha + theta);
  *u = mod2pi(MPI - 2 * alpha);
  *v = mod2pi(phi - *t - *u);

  length_rs = RADCURV * (*t + *u + *v);
//  std::cout << "*" << *t*RADCURV << " " << *u*RADCURV << " " << *v*RADCURV << std::endl;

  return (length_rs);
}

/***********************************************************/
static double c_cc(double RADCURV,
                   double x,
                   double y,
                   double phi,
                   double rs,
                   double rc,
                   double* t,
                   double* u,
                   double* v) {
  double a, b, u1, theta, alpha, length_rs;

  a = x - rs;
  b = y + rc;
  if ((fabs(a) < EPS3) && (fabs(b) < EPS3))
    return (std::numeric_limits<double>::infinity());
  u1 = sqrt(a * a + b * b);
  if (u1 > RADCURVMUL4)
    return (std::numeric_limits<double>::infinity());
  theta = my_atan2(b, a);
  alpha = acos(u1 / RADCURVMUL4);
  *t = mod2pi(MPIDIV2 + alpha + theta);
  *u = mod2pi(MPI - 2 * alpha);
  *v = mod2pi(*t + *u - phi);

  length_rs = RADCURV * (*t + *u + *v);
  return (length_rs);
}

/***********************************************************/
static double csca(double RADCURV,
                   double x,
                   double y,
                   double phi,
                   double rs,
                   double rc,
                   double* t,
                   double* u,
                   double* v) {
  double a, b, length_rs;

  a = x - rs;
  b = y + rc;
  *t = mod2pi(my_atan2(b, a));
  *u = sqrt(a * a + b * b);
  *v = mod2pi(phi - *t);

  length_rs = RADCURV * (*t + *v) + *u;
  return (length_rs);
}

/***********************************************************/
static double cscb(double RADCURV,
                   double x,
                   double y,
                   double phi,
                   double rs,
                   double rc,
                   double* t,
                   double* u,
                   double* v) {
  double a, b, u1, theta, alpha, length_rs;

  a = x + rs;
  b = y - rc;
  u1 = sqrt(a * a + b * b);
  if (u1 < RADCURVMUL2)
    return (std::numeric_limits<double>::infinity());
  theta = my_atan2(b, a);
  *u = sqrt(u1 * u1 - SQRADCURVMUL2);
  alpha = my_atan2(RADCURVMUL2, *u);
  *t = mod2pi(theta + alpha);
  *v = mod2pi(*t - phi);

  length_rs = RADCURV * (*t + *v) + *u;
  return (length_rs);
}

/***********************************************************/
static double ccu_cuc(double RADCURV,
                      double x,
                      double y,
                      double phi,
                      double rs,
                      double rc,
                      double* t,
                      double* u,
                      double* v) {
  double a, b, u1, theta, alpha, length_rs;

  a = x + rs;
  b = y - rc;
  if ((fabs(a) < EPS3) && (fabs(b) < EPS3))
    return (std::numeric_limits<double>::infinity());
  u1 = sqrt(a * a + b * b);
  if (u1 > RADCURVMUL4)
    return (std::numeric_limits<double>::infinity());
  theta = my_atan2(b, a);
  if (u1 > RADCURVMUL2) {
    alpha = acos((u1 / 2 - RADCURV) / RADCURVMUL2);
    *t = mod2pi(MPIDIV2 + theta - alpha);
    *u = mod2pi(MPI - alpha);
    *v = mod2pi(phi - *t + 2 * (*u));
  } else {
    alpha = acos((u1 / 2 + RADCURV) / (RADCURVMUL2));
    *t = mod2pi(MPIDIV2 + theta + alpha);
    *u = mod2pi(alpha);
    *v = mod2pi(phi - *t + 2 * (*u));
  }

  length_rs = RADCURV * (2 * (*u) + *t + *v);
  return (length_rs);
}

/***********************************************************/
static double c_cucu_c(double RADCURV,
                       double x,
                       double y,
                       double phi,
                       double rs,
                       double rc,
                       double* t,
                       double* u,
                       double* v) {
  double a, b, u1, theta, alpha, length_rs, va1, va2;

  a = x + rs;
  b = y - rc;
  if ((fabs(a) < EPS3) && (fabs(b) < EPS3))
    return (std::numeric_limits<double>::infinity());
  u1 = sqrt(a * a + b * b);
  if (u1 > 6 * RADCURV)
    return (std::numeric_limits<double>::infinity());
  theta = my_atan2(b, a);
  va1 = (5 * SQRADCURV - u1 * u1 / 4) / SQRADCURVMUL2;
  if ((va1 < 0.0) || (va1 > 1.0))
    return (std::numeric_limits<double>::infinity());
  *u = acos(va1);
  va2 = sin(*u);
  alpha = asin(RADCURVMUL2 * va2 / u1);
  *t = mod2pi(MPIDIV2 + theta + alpha);
  *v = mod2pi(*t - phi);

  length_rs = RADCURV * (2 * (*u) + *t + *v);
  return (length_rs);
}

/***********************************************************/
static double c_c2sca(double RADCURV,
                      double x,
                      double y,
                      double phi,
                      double rs,
                      double rc,
                      double* t,
                      double* u,
                      double* v) {
  double a, b, u1, theta, alpha, length_rs;

  a = x - rs;
  b = y + rc;
  u1 = sqrt(a * a + b * b);
  if (u1 < RADCURVMUL2)
    return (std::numeric_limits<double>::infinity());
  theta = my_atan2(b, a);
  *u = sqrt(u1 * u1 - SQRADCURVMUL2) - RADCURVMUL2;
  if (*u < 0.0)
    return (std::numeric_limits<double>::infinity());
  alpha = my_atan2(RADCURVMUL2, (*u + RADCURVMUL2));
  *t = mod2pi(MPIDIV2 + theta + alpha);
  *v = mod2pi(*t + MPIDIV2 - phi);

  length_rs = RADCURV * (*t + MPIDIV2 + *v) + *u;
  return (length_rs);
}

/***********************************************************/
static double c_c2scb(double RADCURV,
                      double x,
                      double y,
                      double phi,
                      double rs,
                      double rc,
                      double* t,
                      double* u,
                      double* v) {
  double a, b, u1, theta, length_rs;

  a = x + rs;
  b = y - rc;
  u1 = sqrt(a * a + b * b);
  if (u1 < RADCURVMUL2)
    return (std::numeric_limits<double>::infinity());
  theta = my_atan2(b, a);
  *t = mod2pi(MPIDIV2 + theta);
  *u = u1 - RADCURVMUL2;
  *v = mod2pi(phi - *t - MPIDIV2);

  length_rs = RADCURV * (*t + MPIDIV2 + *v) + *u;
  return (length_rs);
}

/***********************************************************/
static double c_c2sc2_c(double RADCURV,
                        double x,
                        double y,
                        double phi,
                        double rs,
                        double rc,
                        double* t,
                        double* u,
                        double* v) {
  double a, b, u1, theta, alpha, length_rs;

  a = x + rs;
  b = y - rc;
  u1 = sqrt(a * a + b * b);
  if (u1 < RADCURVMUL4)
    return (std::numeric_limits<double>::infinity());
  theta = my_atan2(b, a);
  *u = sqrt(u1 * u1 - SQRADCURVMUL2) - RADCURVMUL4;
  if (*u < 0.0)
    return (std::numeric_limits<double>::infinity());
  alpha = my_atan2(RADCURVMUL2, (*u + RADCURVMUL4));
  *t = mod2pi(MPIDIV2 + theta + alpha);
  *v = mod2pi(*t - phi);

  length_rs = RADCURV * (*t + MPI + *v) + *u;
  return (length_rs);
}

/***********************************************************/
static double cc_c(double RADCURV,
                   double x,
                   double y,
                   double phi,
                   double rs,
                   double rc,
                   double* t,
                   double* u,
                   double* v) {
  double a, b, u1, theta, alpha, length_rs, va;

  a = x - rs;
  b = y + rc;
  if ((fabs(a) < EPS3) && (fabs(b) < EPS3))
    return (std::numeric_limits<double>::infinity());
  u1 = sqrt(a * a + b * b);
  if (u1 > RADCURVMUL4)
    return (std::numeric_limits<double>::infinity());
  theta = my_atan2(b, a);
  *u = acos((8 * SQRADCURV - u1 * u1) / (8 * SQRADCURV));
  va = sin(*u);
  if (fabs(va) < 0.001)
    va = 0.0;
  if ((fabs(va) < 0.001) && (fabs(u1) < 0.001))
    return (std::numeric_limits<double>::infinity());
  alpha = asin(RADCURVMUL2 * va / u1);
  *t = mod2pi(MPIDIV2 - alpha + theta);
  *v = mod2pi(*t - *u - phi);

  length_rs = RADCURV * (*t + *u + *v);
  return (length_rs);
}

/***********************************************************/
static double csc2_ca(double RADCURV,
                      double x,
                      double y,
                      double phi,
                      double rs,
                      double rc,
                      double* t,
                      double* u,
                      double* v) {
  double a, b, u1, theta, alpha, length_rs;

  a = x - rs;
  b = y + rc;
  u1 = sqrt(a * a + b * b);
  if (u1 < RADCURVMUL2)
    return (std::numeric_limits<double>::infinity());
  theta = my_atan2(b, a);
  *u = sqrt(u1 * u1 - SQRADCURVMUL2) - RADCURVMUL2;
  if (*u < 0.0)
    return (std::numeric_limits<double>::infinity());
  alpha = my_atan2((*u + RADCURVMUL2), RADCURVMUL2);
  *t = mod2pi(MPIDIV2 + theta - alpha);
  *v = mod2pi(*t - MPIDIV2 - phi);

  length_rs = RADCURV * (*t + MPIDIV2 + *v) + *u;
  return (length_rs);
}

/***********************************************************/
static double csc2_cb(double RADCURV,
                      double x,
                      double y,
                      double phi,
                      double rs,
                      double rc,
                      double* t,
                      double* u,
                      double* v) {
  double a, b, u1, theta, length_rs;

  a = x + rs;
  b = y - rc;
  u1 = sqrt(a * a + b * b);
  if (u1 < RADCURVMUL2)
    return (std::numeric_limits<double>::infinity());
  theta = my_atan2(b, a);
  *t = mod2pi(theta);
  *u = u1 - RADCURVMUL2;
  *v = mod2pi(-*t - MPIDIV2 + phi);

  length_rs = RADCURV * (*t + MPIDIV2 + *v) + *u;
  return (length_rs);
}

/***********************************************************/

double RS::reeds_shepp(int Nman, std::vector<double>* debug) {
  double x1 = this->ci()->x();
  double y1 = this->ci()->y();
  double t1 = this->ci()->th();
  double x2 = this->cf()->x();
  double y2 = this->cf()->y();
  double t2 = this->cf()->th();

  double RADCURV = 1.0/this->getKmax();

  double x, y, phi;
  double t, u, v, tn, un, vn;
  int num;
  double var, vard, theta, alpha, dx, dy;
  double length = this->_L;
  double sphi, cphi;
  double ap, am, b1, b2;

  /* coordinate change */
  dx = x2 - x1;
  dy = y2 - y1;
  theta = my_atan2(dy, dx);
  alpha = theta - t1;
  vard = sqrt(dx * dx + dy * dy);
  x = cos(alpha) * vard;
  y = sin(alpha) * vard;
  phi = t2 - t1;

  sphi = sin(phi);
  cphi = cos(phi);

  ap = RADCURV * sphi;
  am = -RADCURV * sphi;
  b1 = RADCURV * (cphi - 1);
  b2 = RADCURV * (cphi + 1);

  // There is no num = 0 in this code, so I need to add this
  SAVE_DEBUG(debug, 0.0)
//  SAVE_DEBUG(debug, std::numeric_limits<double>::infinity())

  /*   C | C | C   */

  if (Nman == -1 || Nman == 1)
  {
    length = c_c_c(RADCURV, x, y, phi, ap, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, length)
    num = 1;
    t = tn;
    u = un;
    v = vn;
  }

  if (Nman == -1 || Nman == 2) {
    var = c_c_c(RADCURV, -x, y, -phi, am, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 2;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 3) {
    var = c_c_c(RADCURV, x, -y, -phi, am, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 3;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 4) {
    var = c_c_c(RADCURV, -x, -y, phi, ap, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 4;
      t = tn;
      u = un;
      v = vn;
    }
  }


  /*   C | C C   */

  if (Nman == -1 || Nman == 5) {
    var = c_cc(RADCURV, x, y, phi, ap, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 5;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 6) {
    var = c_cc(RADCURV, -x, y, -phi, am, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 6;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 7) {
    var = c_cc(RADCURV, x, -y, -phi, am, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 7;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 8) {
    var = c_cc(RADCURV, -x, -y, phi, ap, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 8;
      t = tn;
      u = un;
      v = vn;
    }
  }


  /*   C S C   */

  if (Nman == -1 || Nman == 9) {
    var = csca(RADCURV, x, y, phi, ap, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 9;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 10) {
    var = csca(RADCURV, x, -y, -phi, am, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 10;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 11) {
    var = csca(RADCURV, -x, y, -phi, am, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 11;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 12) {
    var = csca(RADCURV, -x, -y, phi, ap, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 12;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 13) {
    var = cscb(RADCURV, x, y, phi, ap, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 13;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 14) {
    var = cscb(RADCURV, x, -y, -phi, am, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 14;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 15) {
    var = cscb(RADCURV, -x, y, -phi, am, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 15;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 16) {
    var = cscb(RADCURV, -x, -y, phi, ap, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 16;
      t = tn;
      u = un;
      v = vn;
    }
  }


  /*   C Cu | Cu C   */

  if (Nman == -1 || Nman == 17) {
    var = ccu_cuc(RADCURV, x, y, phi, ap, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 17;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 18) {
    var = ccu_cuc(RADCURV, x, -y, -phi, am, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 18;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 19) {
    var = ccu_cuc(RADCURV, -x, y, -phi, am, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 19;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 20) {
    var = ccu_cuc(RADCURV, -x, -y, phi, ap, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 20;
      t = tn;
      u = un;
      v = vn;
    }
  }


  /*   C | Cu Cu | C   */

  if (Nman == -1 || Nman == 21) {
    var = c_cucu_c(RADCURV, x, y, phi, ap, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 21;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 22) {
    var = c_cucu_c(RADCURV, x, -y, -phi, am, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 22;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 23) {
    var = c_cucu_c(RADCURV, -x, y, -phi, am, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 23;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 24) {
    var = c_cucu_c(RADCURV, -x, -y, phi, ap, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 24;
      t = tn;
      u = un;
      v = vn;
    }
  }


  /*   C | C2 S C   */

  if (Nman == -1 || Nman == 25) {
    var = c_c2sca(RADCURV, x, y, phi, ap, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 25;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 26) {
    var = c_c2sca(RADCURV, x, -y, -phi, am, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 26;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 27) {
    var = c_c2sca(RADCURV, -x, y, -phi, am, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 27;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 28) {
    var = c_c2sca(RADCURV, -x, -y, phi, ap, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 28;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 29) {
    var = c_c2scb(RADCURV, x, y, phi, ap, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 29;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 30) {
    var = c_c2scb(RADCURV, x, -y, -phi, am, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 30;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 31) {
    var = c_c2scb(RADCURV, -x, y, -phi, am, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 31;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 32) {
    var = c_c2scb(RADCURV, -x, -y, phi, ap, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 32;
      t = tn;
      u = un;
      v = vn;
    }
  }


  /*   C | C2 S C2 | C   */

  if (Nman == -1 || Nman == 33) {
    var = c_c2sc2_c(RADCURV, x, y, phi, ap, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 33;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 34) {
    var = c_c2sc2_c(RADCURV, x, -y, -phi, am, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 34;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 35) {
    var = c_c2sc2_c(RADCURV, -x, y, -phi, am, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 35;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 36) {
    var = c_c2sc2_c(RADCURV, -x, -y, phi, ap, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 36;
      t = tn;
      u = un;
      v = vn;
    }
  }


  /*   C C | C   */

  if (Nman == -1 || Nman == 37) {
    var = cc_c(RADCURV, x, y, phi, ap, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 37;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 38) {
    var = cc_c(RADCURV, x, -y, -phi, am, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 38;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 39) {
    var = cc_c(RADCURV, -x, y, -phi, am, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 39;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 40) {
    var = cc_c(RADCURV, -x, -y, phi, ap, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 40;
      t = tn;
      u = un;
      v = vn;
    }
  }


  /*   C S C2 | C   */

  if (Nman == -1 || Nman == 41) {
    var = csc2_ca(RADCURV, x, y, phi, ap, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 41;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 42) {
    var = csc2_ca(RADCURV, x, -y, -phi, am, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 42;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 43) {
    var = csc2_ca(RADCURV, -x, y, -phi, am, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 43;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 44) {
    var = csc2_ca(RADCURV, -x, -y, phi, ap, b1, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 44;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 45) {
    var = csc2_cb(RADCURV, x, y, phi, ap, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 45;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 46) {
    var = csc2_cb(RADCURV, x, -y, -phi, am, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 46;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 47) {
    var = csc2_cb(RADCURV, -x, y, -phi, am, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 47;
      t = tn;
      u = un;
      v = vn;
    }
  }


  if (Nman == -1 || Nman == 48) {
    var = csc2_cb(RADCURV, -x, -y, phi, ap, b2, &tn, &un, &vn);
    SAVE_DEBUG(debug, var)
    if (var < length) {
      length = var;
      num = 48;
      t = tn;
      u = un;
      v = vn;
    }
  }


  this->_L1 = t * RADCURV;
  this->_L2 = u * RADCURV;
  this->_L3 = v * RADCURV;
  this->_Nman = num;
  this->_L = length;
  return length;
}

void RS::draw() {

}

#endif // CUDA_ON