#ifndef DUBINS_CUH
#define DUBINS_CUH

#include <curve.cuh>
#include <utils.cuh>
#include <constants.cuh>

#define DUBINS_DEFAULT_KMAX 0.01

#include <cmath>
#include <limits>

class Dubins : public Curve {
public:
  enum D_TYPE {INVALID, LSL, RSR, LSR, RSL, RLR, LRL}; ///<The possible types of Dubins.
private:
  D_TYPE _type; ///<The possible types of Dubins.
  K_T _kmax=0.0, _k1=0.0, _k2=0.0, _k3=0.0; ///<The maximum curvature and the curvature for each part of the Dubins.
  LEN_T _s1=0.0, _s2=0.0, _s3=0.0; ///<The lengths of each part of the Dubins.

  /*!
   * Function to standardize the components. Credit to Marco Frego & Paolo Bevilacqua
   */
  BOTH void scaleToStandard(Angle& phi, real_type& lambda, Angle& sth0, Angle& sth1, K_T& sKmax);

  /*!
   * Given the standardized version, compute the best word. Credit to Marco Frego & Paolo Bevilacqua.
   * @param th0 The initial standardized angle.
   * @param th1 The final standardized angle.
   * @param lambda A multiplier.
   * @param sKmax The standardized curvature.
   */
  BOTH void computeBest( Angle th0, Angle th1, real_type lambda, K_T& sKmax);

  /*!
   * Function to solve the Dubins curve.
   */
  BOTH void solve(){
    real_type lambda;
    K_T sKmax;
    Angle phi, sth0, sth1;
    scaleToStandard(phi, lambda, sth0, sth1, sKmax);
    computeBest(sth0, sth1, lambda, sKmax);
  }

public:
  /*!
   * Void constructor to initialize a Dubins object.
   */
  Dubins() :
          Curve(CURVE_TYPE::DUBINS),
          _type (D_TYPE::INVALID),
          _kmax(0) {}

  /*!
   * Constructor to initialize a Dubins object with an initial and a final `Configuration2` and additional possible parameters.
   * @param ci The initial `Configuration2`.
   * @param cf The final `Configuration2`
   * @param params Additional parameters to pass. Default is `nullptr`, in such case DUBINS_DEFAULT_KMAX==0.01 is used.
   */
  BOTH Dubins(Configuration2 ci, Configuration2 cf, real_type* params=nullptr) :
    Curve(ci, cf, CURVE_TYPE::DUBINS, params),
    _type(D_TYPE::INVALID)
  {  
    if (params==nullptr) { this->_kmax=DUBINS_DEFAULT_KMAX; }
    else                 { this->_kmax=params[0]; }
    solve();
  }

  /*!
   * Constructor to initialize a Dubins object with an initial and a final `Configuration2` and additional possible parameters.
   * @param ci The initial `Configuration2`.
   * @param cf The final `Configuration2`
   * @param kmax The curvature of the Dubins parts.
   */
  BOTH Dubins(Configuration2 ci, Configuration2 cf, real_type kmax) :
    Curve(ci, cf, CURVE_TYPE::DUBINS),
    _type(D_TYPE::INVALID),
    _kmax(kmax)
  {
    solve();
  }

  BOTH K_T kmax() const { return this->_kmax; }                                ///<Returns the maximum curvature.
  BOTH K_T k1() const { return this->_k1; }                                    ///<Returns the curvature of the first part of the Dubins.
  BOTH K_T k2() const { return this->_k2; }                                    ///<Returns the curvature of the middle part of the Dubins.
  BOTH K_T k3() const { return this->_k3; }                                    ///<Returns the curvature of the final part of the Dubins.
  BOTH LEN_T s1() const { return this->_s1; }                                  ///<Returns the length of the first part of the Dubins.
  BOTH LEN_T s2() const { return this->_s2; }                                  ///<Returns the length of the first part of the Dubins.
  BOTH LEN_T s3() const { return this->_s3; }                                  ///<Returns the length of the first part of the Dubins.
  BOTH LEN_T l() const override { return (this->s1()+this->s2()+this->s3()); } ///<Returns the length of the Dubins.
  BOTH D_TYPE type() const { return this->_type; }                             ///<Returns the word of the Dubins.

  BOTH K_T kmax(K_T kmax) { this->_kmax = kmax; return this->kmax(); }         ///<Sets the maximum curvature and returns the new set value.
  BOTH K_T k1(K_T k1) { this->_k1 = k1; return this->k1(); }                   ///<Sets the curvature of the first part of the Dubins and returns the new set value.
  BOTH K_T k2(K_T k2) { this->_k2 = k2; return this->k2(); }                   ///<Sets the curvature of the middle part of the Dubins and returns the new set value.
  BOTH K_T k3(K_T k3) { this->_k3 = k3; return this->k3(); }                   ///<Sets the curvature of the final part of the Dubins and returns the new set value.
  BOTH LEN_T s1(LEN_T s1) { this->_s1 = s1; return this->s1(); }               ///<Sets the length of the first part of the Dubins and returns the new set value.
  BOTH LEN_T s2(LEN_T s2) { this->_s2 = s2; return this->s2(); }               ///<Sets the length of the middle part of the Dubins and returns the new set value.
  BOTH LEN_T s3(LEN_T s3) { this->_s3 = s3; return this->s3(); }               ///<Sets the length of the final part of the Dubins and returns the new set value.
  BOTH D_TYPE type(D_TYPE type) { this->_type = type; return this->type(); }   ///<Sets the word of the Dubins.

  /*!
   * Function to print the stringy word of the Dubins.
   * @return A string containing the word of the computed Dubins.
   */
  string type_to_string() {
    string ret="INVALID";
    switch(this->_type){
      case D_TYPE::LSL: { ret="LSL"; break; }
      case D_TYPE::LSR: { ret="LSR"; break; }
      case D_TYPE::RSR: { ret="RSR"; break; }
      case D_TYPE::RSL: { ret="RSL"; break; }
      case D_TYPE::LRL: { ret="LRL"; break; }
      case D_TYPE::RLR: { ret="RLR"; break; }
    }
    return ret;
  }

  /*!
   * Function to print the most essential info about `Dubins`.
   * @param str An additional string to add at the beginning.
   * @return A `std::stringstream` object containing the data of `Dubins`.
   */
  std::stringstream to_string (const std::string& str="") {
    std::stringstream out;
    out << "c0: " << this->ci()->to_string().str() << "\tc1: " << this->cf()->to_string().str() << "\tk: " << this->kmax() << "\tl: " << this->l();
    return out;
  }

  /*! This function overload the << operator so to print with `std::cout` the most essential info about the `Dubins`.
  		\param[in] out The out stream.
  		\param[in] data The Dubins to print.
  		\returns An output stream to be printed.
  */
  friend std::ostream& operator<<(std::ostream &out, Dubins& data) {
    out << data.to_string().str();
    return out;
  }
  
};

#endif //DUBINS_CUH
