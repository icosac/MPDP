#ifndef CLOTHOID_HH
#define CLOTHOID_HH

#include <cmath>
#include <iostream>
using namespace std;

#include <curve.cuh>
#include <utils.cuh>

#define A_THRESOLD   0.01  //TODO cosa sono questi?
#define A_SERIE_SIZE 3     //TODO cosa sono questi?

template<class T1>
class ClothoidG1 : public Curve<T1>{
private:
  LEN_T _l;                          ///< Length
  real_type tol = 0.000000000001;   ///< Tolerance
  bool compute_deriv=false;         ///< Whether to compute derivatives or not
  K_T _dk;                          ///< Curvature derivative
  real_type L_D[2];                 ///< Derivateves of the length
  real_type k_D[2];                 ///< Curvature derivatives
  real_type dk_D[2];                ///< Sharpness derivatives

public:
  ClothoidG1() : Curve<T1>(CURVE_TYPE::CLOTHOID), _l(0) {}
  ClothoidG1(Configuration2<T1> ci, Configuration2<T1> cf, LEN_T l=0) : Curve<T1>(ci, cf, CURVE_TYPE::CLOTHOID), _l(l) {
    this->buildG1();
  }

  LEN_T l() const { return this->_l; }                  ///< Returns the length of the curve.
  K_T dk () const { return this->_dk; }                 ///< Returns the curvature derivative.
//  real_type k () { return ; } //  TODO how to compute these?

  K_T dk(real_type dk) {this->_dk=dk; return this->dk(); }   ///< Sets and returns the curvature derivative.
  LEN_T l(LEN_T l) { this->_l=l; return this->l(); }         ///< Sets and returns the length of the curve.

  int buildG1(){
    static real_type const CF[] = { 2.989696028701907,  0.716228953608281,
                                   -0.458969738821509, -0.502821153340377,
                                    0.261062141752652, -0.045854475238709 } ;
    // Translation (0,0)
    real_type dx   = this->cf()->x() - this->ci()->x() ;
    real_type dy   = this->cf()->y() - this->ci()->y() ;
    real_type r    = hypot( dx, dy ) ;
    real_type phi  = atan2( dy, dx ) ;
    real_type phi0 = this->ci()->th() - phi ;
    real_type phi1 = this->cf()->th() - phi ;

    phi0 -= m_2pi*round(phi0/m_2pi) ;
    phi1 -= m_2pi*round(phi1/m_2pi) ;

    if      ( phi0 >  m_pi ) phi0 -= m_2pi ;
    else if ( phi0 < -m_pi ) phi0 += m_2pi ;
    if      ( phi1 >  m_pi ) phi1 -= m_2pi ;
    else if ( phi1 < -m_pi ) phi1 += m_2pi ;

    real_type delta = phi1 - phi0 ;

    // Initial point
    real_type X  = phi0*m_1_pi ;
    real_type Y  = phi1*m_1_pi ;
    real_type xy = X*Y ;
    Y *= Y ; X *= X ;
    real_type A = (phi0+phi1)*(CF[0]+xy*(CF[1]+xy*CF[2])+(CF[3]+xy*CF[4])*(X+Y)+CF[5]*(X*X+Y*Y)) ;

    // Newton
    real_type g=0, dg, intC[3], intS[3] ;
    int_type  niter = 0 ;
    do {
      GeneralizedFresnelCS( 3, 2*A, delta-A, phi0, intC, intS ) ;
      g   = intS[0] ;
      dg  = intC[2] - intC[1] ;
      A  -= g / dg ;
    } while ( ++niter <= 10 && std::abs(g) > tol ) ;

    ASSERT( std::abs(g) <= this->tol, "Newton do not converge, g = " << g << " niter = " << niter ) ;
    GeneralizedFresnelCS( 2*A, delta-A, phi0, intC[0], intS[0] ) ;
    LEN_T L=this->l(r/intC[0]);

    ASSERT( L > 0, "Negative length L = " << L ) ;
    K_T kappa0=this->ci()->k((delta-A)/L);
    K_T dk=this->dk(2*A/L/L);

    if ( compute_deriv ) {
      real_type alpha = intC[0]*intC[1] + intS[0]*intS[1] ;
      real_type beta  = intC[0]*intC[2] + intS[0]*intS[2] ;
      real_type gamma = intC[0]*intC[0] + intS[0]*intS[0] ;
      real_type tx    = intC[1]-intC[2] ;
      real_type ty    = intS[1]-intS[2] ;
      real_type txy   = L*(intC[1]*intS[2]-intC[2]*intS[1]) ;
      real_type omega = L*(intS[0]*tx-intC[0]*ty) - txy ;

      delta = intC[0]*tx + intS[0]*ty ;

      L_D[0] = omega/delta ;
      L_D[1] = txy/delta ;

      delta *= L ;
      k_D[0] = (beta-gamma-kappa0*omega)/delta ;
      k_D[1] = -(beta+kappa0*txy)/delta ;

      delta  *= L/2 ;
      dk_D[0] = (gamma-alpha-dk*omega*L)/delta ;
      dk_D[1] = (alpha-dk*txy*L)/delta ;
    }
    return niter ;
  }

  void
  GeneralizedFresnelCS( real_type   a,
                        real_type   b,
                        real_type   c,
                        real_type & intC,
                        real_type & intS ) {
    real_type xx, yy ;
    if ( abs(a) < A_THRESOLD ) evalXYaSmall( a, b, A_SERIE_SIZE, xx, yy ) ;
    else                       evalXYaLarge( a, b, xx, yy ) ;

    real_type cosc = cos(c) ;
    real_type sinc = sin(c) ;

    intC = xx * cosc - yy * sinc ;
    intS = xx * sinc + yy * cosc ;
  }

  void
  GeneralizedFresnelCS( int_type  nk,
                        real_type a,
                        real_type b,
                        real_type c,
                        real_type intC[],
                        real_type intS[] ) {
    ASSERT( nk > 0 && nk < 4, "nk = " << nk << " must be in 1..3" ) ;

    if ( abs(a) < A_THRESOLD ) evalXYaSmall( nk, a, b, A_SERIE_SIZE, intC, intS ) ;
    else                       evalXYaLarge( nk, a, b, intC, intS ) ;

    real_type cosc = cos(c) ;
    real_type sinc = sin(c) ;

    for ( int_type k = 0 ; k < nk ; ++k ) {
      real_type xx = intC[k] ;
      real_type yy = intS[k] ;
      intC[k] = xx * cosc - yy * sinc ;
      intS[k] = xx * sinc + yy * cosc ;
    }
  }

  void FresnelCS(real_type x, real_type &C, real_type &S){
    // C = int_0^x cos(pi/2*t^2) dt
    // S = int_0^x sin(pi/2*t^2) dt
    // Computing Fresnel integrals via modified trapezium rules
    // Mohammad Alazah � Simon N. Chandler-Wilde � Scott La Porte
    // Numer. Math. (2014) 128:635�661
    // DOI 10.1007/s00211-014-0627-z

    const int N = 12;
    double h = std::sqrt(m_pi / (N + 0.5));
    double AN = m_pi / h;
    double rootpi = std::sqrt(m_pi);
    double t[N], t2[N], t4[N], et2[N];
    double x2pi2 = 0.5 * m_pi * x * x;
    double x4 = x2pi2 * x2pi2;

    for (int i = N; i > 0; --i) {
      t[i - 1] = h * (i - 0.5);
      t2[i - 1] = t[i - 1] * t[i - 1];
      t4[i - 1] = t2[i - 1] * t2[i - 1];
      et2[i - 1] = std::exp(-t2[i - 1]);
    }

    double a = et2[N - 1] / (x4 + t4[N - 1]);
    double b = t2[N - 1] * a;

    double term;
    for (int i = 2; i <= N; ++i) {
      term = et2[N - i] / (x4 + t4[N - i]);
      a = a + term;
      b = b + t2[N - i] * term;
    }
    a = a * x2pi2;
    double mx = rootpi * AN * x;
    double Mx = rootpi / AN * x;
    double Chalf = 0.5 * mx / abs(mx);
    double Shalf = Chalf;

    if (abs(mx) < 39) {
      double mxs = mx;
      double shx = std::sinh(mxs);
      double sx = std::sin(mxs);
      double den = 0.5 / (std::cos(mxs) + std::cosh(mxs));
      Chalf = (shx + sx) * den;
      double ssdiff = shx - sx;
      if (abs(mxs) < 1) {
        double mxs3 = mxs * mxs * mxs;
        double mxs4 = mxs3 * mxs;
        ssdiff = mxs3 * (1 / 3. + mxs4 * (1 / 2520. + mxs4 * ((1 / 19958400.) + (0.001 / 653837184.) * mxs4)));
      }
      Shalf = ssdiff * den;
    }
    double cx2 = std::cos(x2pi2);
    double sx2 = std::sin(x2pi2);
    C = Chalf + Mx * (a * sx2 - b * cx2);
    S = Shalf - Mx * (a * cx2 + b * sx2);
  }

  void FresnelCS(int_type  nk,
                 real_type t,
                 real_type C[],
                 real_type S[] ) {
    FresnelCS(t,C[0],S[0]) ;
    if ( nk > 1 ) {
      real_type tt = m_pi_2*(t*t) ;
      real_type ss = sin(tt) ;
      real_type cc = cos(tt) ;
      C[1] = ss*m_1_pi ;
      S[1] = (1-cc)*m_1_pi ;
      if ( nk > 2 ) {
        C[2] = (t*ss-S[0])*m_1_pi ;
        S[2] = (C[0]-t*cc)*m_1_pi ;
      }
    }
  }

  real_type LommelReduced(real_type mu, real_type nu, real_type b) {
    real_type tmp = 1/((mu+nu+1)*(mu-nu+1)) ;
    real_type res = tmp ;
    for ( int_type n = 1 ; n <= 100 ; ++n ) {
      tmp *= (-b/(2*n+mu-nu+1)) * (b/(2*n+mu+nu+1)) ;
      res += tmp ;
      if ( abs(tmp) < abs(res) * 1e-50 ) break ;
    }
    return res ;
  }

  void evalXYaLarge(real_type   a,
                    real_type   b,
                    real_type & X,
                    real_type & Y ) {
    real_type s = a > 0 ? +1 : -1;
    real_type absa = abs(a);
    real_type z = m_1_sqrt_pi * sqrt(absa);
    real_type ell = s * b * m_1_sqrt_pi / sqrt(absa);
    real_type g = -0.5 * s * (b * b) / absa;
    real_type cg = cos(g) / z;
    real_type sg = sin(g) / z;

    real_type Cl, Sl, Cz, Sz;
    this->FresnelCS(ell, Cl, Sl);
    this->FresnelCS(ell + z, Cz, Sz);

    real_type dC0 = Cz - Cl;
    real_type dS0 = Sz - Sl;

    X = cg * dC0 - s * sg * dS0;
    Y = sg * dC0 + s * cg * dS0;
  }

  void evalXYaLarge(int_type  nk,
                    real_type a,
                    real_type b,
                    real_type X[],
                    real_type Y[] ) {

    ASSERT( nk < 4 && nk > 0,
                  "In evalXYaLarge first argument nk must be in 1..3, nk " << nk ) ;

    real_type s    = a > 0 ? +1 : -1 ;
    real_type absa = abs(a) ;
    real_type z    = m_1_sqrt_pi*sqrt(absa) ;
    real_type ell  = s*b*m_1_sqrt_pi/sqrt(absa) ;
    real_type g    = -0.5*s*(b*b)/absa ;
    real_type cg   = cos(g)/z ;
    real_type sg   = sin(g)/z ;

    real_type Cl[3], Sl[3], Cz[3], Sz[3] ;

    this->FresnelCS( nk, ell,   Cl, Sl ) ;
    this->FresnelCS( nk, ell+z, Cz, Sz ) ;

    real_type dC0 = Cz[0] - Cl[0] ;
    real_type dS0 = Sz[0] - Sl[0] ;
    X[0] = cg * dC0 - s * sg * dS0 ;
    Y[0] = sg * dC0 + s * cg * dS0 ;
    if ( nk > 1 ) {
      cg /= z ;
      sg /= z ;
      real_type dC1 = Cz[1] - Cl[1] ;
      real_type dS1 = Sz[1] - Sl[1] ;
      real_type DC  = dC1-ell*dC0 ;
      real_type DS  = dS1-ell*dS0 ;
      X[1] = cg * DC - s * sg * DS ;
      Y[1] = sg * DC + s * cg * DS ;
      if ( nk > 2 ) {
        real_type dC2 = Cz[2] - Cl[2] ;
        real_type dS2 = Sz[2] - Sl[2] ;
        DC   = dC2+ell*(ell*dC0-2*dC1) ;
        DS   = dS2+ell*(ell*dS0-2*dS1) ;
        cg   = cg/z ;
        sg   = sg/z ;
        X[2] = cg * DC - s * sg * DS ;
        Y[2] = sg * DC + s * cg * DS ;
      }
    }
  }

  void  evalXYaZero(int_type nk,
                   real_type b,
                   real_type X[],
                   real_type Y[] ) {
    real_type sb = sin(b) ;
    real_type cb = cos(b) ;
    real_type b2 = b*b ;
    if ( abs(b) < 1e-3 ) {
      X[0] = 1-(b2/6)*(1-(b2/20)*(1-(b2/42))) ;
      Y[0] = (b/2)*(1-(b2/12)*(1-(b2/30))) ;
    } else {
      X[0] = sb/b ;
      Y[0] = (1-cb)/b ;
    }
    // use recurrence in the stable part
    int_type m = int_type(floor(2*b)) ;
    if ( m >= nk ) m = nk-1 ;
    if ( m < 1   ) m = 1 ;
    for ( int_type k = 1 ; k < m ; ++k ) {
      X[k] = (sb-k*Y[k-1])/b ;
      Y[k] = (k*X[k-1]-cb)/b ;
    }
    //  use Lommel for the unstable part
    if ( m < nk ) {
      real_type A   = b*sb ;
      real_type D   = sb-b*cb ;
      real_type B   = b*D ;
      real_type C   = -b2*sb ;
      real_type rLa = LommelReduced(m+0.5,1.5,b) ;
      real_type rLd = LommelReduced(m+0.5,0.5,b) ;
      for ( int_type k = m ; k < nk ; ++k ) {
        real_type rLb = LommelReduced(k+1.5,0.5,b) ;
        real_type rLc = LommelReduced(k+1.5,1.5,b) ;
        X[k] = ( k*A*rLa + B*rLb + cb ) / (1+k) ;
        Y[k] = ( C*rLc + sb ) / (2+k) + D*rLd ;
	      rLa  = rLc ;
  	    rLd  = rLb ;
      }
    }
  }

  void  evalXYaSmall( real_type   a,
                      real_type   b,
                      int_type    p,
                      real_type & X,
                      real_type & Y ) {
    ASSERT( p < 11 && p > 0,
                  "In evalXYaSmall p = " << p << " must be in 1..10" ) ;

    real_type X0[43], Y0[43] ;

    int_type nkk = 4*p + 3 ; // max 43
    this->evalXYaZero(nkk, b, X0, Y0) ;

    X = X0[0]-(a/2)*Y0[2] ;
    Y = Y0[0]+(a/2)*X0[2] ;

    real_type t  = 1 ;
    real_type aa = -a*a/4 ; // controllare!
    for ( int_type n=1 ; n <= p ; ++n ) {
      t *= aa/(2*n*(2*n-1)) ;
      real_type bf = a/(4*n+2) ;
      int_type  jj = 4*n ;
      X += t*(X0[jj]-bf*Y0[jj+2]) ;
      Y += t*(Y0[jj]+bf*X0[jj+2]) ;
    }
  }

  void
  evalXYaSmall( int_type  nk,
                real_type a,
                real_type b,
                int_type  p,
                real_type X[],
                real_type Y[] ) {
    int_type  nkk = nk + 4*p + 2 ; // max 45
    real_type X0[45], Y0[45] ;

    ASSERT( nkk < 46,
                  "In evalXYaSmall (nk,p) = (" << nk << "," << p << ")\n" <<
                  "nk + 4*p + 2 = " << nkk  << " must be less than 46\n" ) ;

    this->evalXYaZero(nkk, b, X0, Y0 ) ;

    for ( int_type j=0 ; j < nk ; ++j ) {
      X[j] = X0[j]-(a/2)*Y0[j+2] ;
      Y[j] = Y0[j]+(a/2)*X0[j+2] ;
    }

    real_type t  = 1 ;
    real_type aa = -a*a/4 ; // controllare!
    for ( int_type n=1 ; n <= p ; ++n ) {
      t *= aa/(2*n*(2*n-1)) ;
      real_type bf = a/(4*n+2) ;
      for ( int_type j = 0 ; j < nk ; ++j ) {
        int_type jj = 4*n+j ;
        X[j] += t*(X0[jj]-bf*Y0[jj+2]) ;
        Y[j] += t*(Y0[jj]+bf*X0[jj+2]) ;
      }
    }
  }
};


#endif //CLOTHOID_HH
