#ifndef DUBINS_HH
#define DUBINS_HH

// Library includes
#include <curve.hh>
#include <utils.hh>

#define DUBINS_DEFAULT_KMAX 0.01

// System includes
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
  void scaleToStandard(Angle& phi, real_type& lambda, Angle& sth0, Angle& sth1, K_T& sKmax);

  /*!
   * Given the standardized version, compute the best word. Credit to Marco Frego & Paolo Bevilacqua.
   * @param th0 The initial standardized angle.
   * @param th1 The final standardized angle.
   * @param lambda A multiplier.
   * @param sKmax The standardized curvature.
   */
  void computeBest( Angle th0, Angle th1, real_type lambda, K_T& sKmax);

  /*!
   * Function to solve the Dubins curve.
   */
  void solve(){
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
  Dubins(Configuration2 ci, Configuration2 cf, std::vector<real_type> params = {}) :
    Curve(ci, cf, CURVE_TYPE::DUBINS, params),
    _type(D_TYPE::INVALID)
  {
    if (params.empty())   { this->_kmax = DUBINS_DEFAULT_KMAX; }
    else                  { this->_kmax = params[0]; }
    solve();
  }

  /*!
   * Constructor to initialize a Dubins object with an initial and a final `Configuration2` and additional possible parameters.
   * @param x0 The initial abscissa.
   * @param y0 The final ordinate.
   * @param th0 The initial angle.
   * @param x1 The initial abscissa.
   * @param y1 The final ordinate.
   * @param th1 The final angle.
   * @param params Additional parameters to pass. Default is `nullptr`, in such case DUBINS_DEFAULT_KMAX==0.01 is used.
   */
  Dubins(real_type x0, real_type y0, Angle th0, real_type x1, real_type y1, real_type th1, std::vector<real_type> params = {}) :
    Curve(Configuration2(x0, y0, th0), Configuration2(x1, y1, th1), CURVE_TYPE::DUBINS, params),
    _type(D_TYPE::INVALID)
  {
    if (params.empty())  { this->_kmax = DUBINS_DEFAULT_KMAX; }
    else                  { this->_kmax = params[0]; }
    solve();
  }

  /*!
   * Constructor to initialize a Dubins object with an initial and a final `Configuration2` and additional possible parameters.
   * @param ci The initial `Configuration2`.
   * @param cf The final `Configuration2`
   * @param kmax The curvature of the Dubins parts.
   */
  Dubins(Configuration2 ci, Configuration2 cf, real_type kmax) :
    Curve(ci, cf, CURVE_TYPE::DUBINS),
    _type(D_TYPE::INVALID),
    _kmax(kmax) 
  {
    solve();
  }

  K_T kmax() const { return this->_kmax; }                                 ///<Returns the maximum curvature.
  K_T k1() const { return this->_k1; }                                     ///<Returns the curvature of the first part of the Dubins.
  K_T k2() const { return this->_k2; }                                     ///<Returns the curvature of the middle part of the Dubins.
  K_T k3() const { return this->_k3; }                                     ///<Returns the curvature of the final part of the Dubins.
  LEN_T s1() const { return this->_s1; }                                   ///<Returns the length of the first part of the Dubins.
  LEN_T s2() const { return this->_s2; }                                   ///<Returns the length of the first part of the Dubins.
  LEN_T s3() const { return this->_s3; }                                   ///<Returns the length of the first part of the Dubins.
  LEN_T l() const override { return (this->s1()+this->s2()+this->s3()); }  ///<Returns the length of the Dubins.
  D_TYPE dtype() const { return this->_type; }                             ///<Returns the word of the Dubins.

  K_T kmax(K_T kmax) { this->_kmax = kmax; return this->kmax(); }          ///<Sets the maximum curvature and returns the new set value.
  K_T k1(K_T k1) { this->_k1 = k1; return this->k1(); }                    ///<Sets the curvature of the first part of the Dubins and returns the new set value.
  K_T k2(K_T k2) { this->_k2 = k2; return this->k2(); }                    ///<Sets the curvature of the middle part of the Dubins and returns the new set value.
  K_T k3(K_T k3) { this->_k3 = k3; return this->k3(); }                    ///<Sets the curvature of the final part of the Dubins and returns the new set value.
  LEN_T s1(LEN_T s1) { this->_s1 = s1; return this->s1(); }                ///<Sets the length of the first part of the Dubins and returns the new set value.
  LEN_T s2(LEN_T s2) { this->_s2 = s2; return this->s2(); }                ///<Sets the length of the middle part of the Dubins and returns the new set value.
  LEN_T s3(LEN_T s3) { this->_s3 = s3; return this->s3(); }                ///<Sets the length of the final part of the Dubins and returns the new set value.
  D_TYPE dtype(D_TYPE type) { this->_type = type; return this->dtype(); }  ///<Sets the word of the Dubins.

  /*!
   * Function to print the stringy word of the Dubins.
   * @return A string containing the word of the computed Dubins.
   */
  std::string type_to_string() {
    std::string ret="INVALID";
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
    out << (str.empty() ? "" : str+" ") << "c0: " << this->ci()->to_string().str() << "\tc1: " << this->cf()->to_string().str() << "\tk: " << this->kmax() << "\tl: " << this->l() << "\ttype: " << this->type_to_string();
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

  void draw() override;
};

#endif //DUBINS_HH
