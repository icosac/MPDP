#ifndef DUBINS_HH
#define DUBINS_HH

#include <curve.cuh>
#include <utils.cuh>
#include <settings.hh>
#include <constants.cuh>
#include <trig.cuh>


#define DUBINS_DEFAULT_KMAX 0.01

#include <cmath>
#include <limits>

template<class T1>
class Dubins : public Curve<T1> {
public:
  enum D_TYPE {INVALID, LSL, RSR, LSR, RSL, RLR, LRL};
private:
  D_TYPE _Dtype;
  K_T _kmax=0.0, _k1=0.0, _k2=0.0, _k3=0.0;
  LEN_T _s1=0.0, _s2=0.0, _s3=0.0;

  BOTH void solve(){
    real_type lambda;
    K_T sKmax;
    Angle phi, sth0, sth1;
    scaleToStandard(phi, lambda, sth0, sth1, sKmax);
    computeBest(sth0, sth1, lambda, sKmax);
  }

public:
  Dubins() :
          Curve<T1>(CURVE_TYPE::DUBINS),
          _Dtype (D_TYPE::INVALID),
          _kmax(0) {}

  BOTH Dubins(Configuration2<T1> ci, Configuration2<T1> cf, real_type* params) :
          Curve<T1>(ci, cf, CURVE_TYPE::DUBINS, params),
          _Dtype(D_TYPE::INVALID) 
  {  
    if (params!=NULL) { this->_kmax=params[0]; }
    else              { this->_kmax=DUBINS_DEFAULT_KMAX; }
    solve();
  }

  BOTH Dubins(Configuration2<T1> ci, Configuration2<T1> cf, real_type kmax) :
          Curve<T1>(ci, cf, CURVE_TYPE::DUBINS),
          _Dtype(D_TYPE::INVALID),
          _kmax(kmax) 
  {
    solve();
  }

  BOTH K_T kmax() const { return this->_kmax; }
  BOTH K_T k1() const { return this->_k1; }
  BOTH K_T k2() const { return this->_k2; }
  BOTH K_T k3() const { return this->_k3; }
  BOTH LEN_T s1() const { return this->_s1; }
  BOTH LEN_T s2() const { return this->_s2; }
  BOTH LEN_T s3() const { return this->_s3; }
  BOTH LEN_T l() const { return (this->s1()+this->s2()+this->s3()); }
  BOTH D_TYPE type() const { return this->_Dtype; }

  BOTH K_T kmax(K_T kmax) { this->_kmax = kmax; return this->kmax(); }
  BOTH K_T k1(K_T k1) { this->_k1 = k1; return this->k1(); }
  BOTH K_T k2(K_T k2) { this->_k2 = k2; return this->k2(); }
  BOTH K_T k3(K_T k3) { this->_k3 = k3; return this->k3(); }
  BOTH LEN_T s1(LEN_T s1) { this->_s1 = s1; return this->s1(); }
  BOTH LEN_T s2(LEN_T s2) { this->_s2 = s2; return this->s2(); }
  BOTH LEN_T s3(LEN_T s3) { this->_s3 = s3; return this->s3(); }
  BOTH D_TYPE type(D_TYPE type) { this->_Dtype = type; return this->type(); }

  BOTH void scaleToStandard(Angle& phi, real_type& lambda, Angle& sth0, Angle& sth1, K_T& sKmax);
  BOTH void computeBest( Angle th0, Angle th1, real_type lambda, K_T& sKmax);

  // std::stringstream to_string (std::string str="") const {
  //   std::stringstream out;
  //   out << "Section1:\n\tx0: " << x0 << endl << "\ty0: " << y0 << endl << "\tth0: " << th0 << k
  //   return out;
  // }

  // /*! This function overload the << operator so to print with `std::cout` the most essential info about the `Configuration2`.
  // 		\param[in] out The out stream.
  // 		\param[in] data The configuration to print.
  // 		\returns An output stream to be printed.
  // */
  // friend std::ostream& operator<<(std::ostream &out, const Dubins& data) {
  //   out << data.to_string().str();
  //   return out;
  // }
  
};

#include <dubins.cut>

#endif //DUBINS_HH
