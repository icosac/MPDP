#ifndef CURVE_HH
#define CURVE_HH

#include <configuration.cuh>

class Curve{
public: 
  enum CURVE_TYPE { INVALID, CLOTHOID, DUBINS, DUBINS_ARC }; ///< Possible types of CURVE
  
private:
  Configuration2 _ci; ///<Initial `Configuration`
  Configuration2 _cf; ///<Final `Configuration`
  CURVE_TYPE _type;   ///<Type of curve
  real_type* _params; ///<Parameters of curve

public:
  /*!
   * @brief Void constructor.
   */
  BOTH Curve() : _ci(Configuration2()), _cf(Configuration2()), _type(CURVE_TYPE::INVALID), _params(NULL) {}
  /*!
   * @brief Constructor to only set the type of the curve.
   */
  BOTH Curve(CURVE_TYPE type=CURVE_TYPE::INVALID) : _ci(), _cf(), _type(type), _params(NULL) {}

  /*!
   * @brief Constructor that takes two `Configuration2` and the type of the curve.
   * @param ci Initial configuration.
   * @param cf Final configuration.
   * @param type Type of the curve.
   * @param params The parameters of the curve, such as the curvature.
   */
  BOTH Curve(Configuration2 ci, Configuration2 cf, CURVE_TYPE type=CURVE_TYPE::INVALID, real_type* params=NULL) : _ci(ci), _cf(cf), _type(type), _params(params) {}

  BOTH Configuration2* ci() { return &(this->_ci); }   ///< Returns a pointer to the initial `Configuration2`.
  BOTH Configuration2* cf() { return &(this->_cf); }   ///< Returns a pointer to the final `Configuration2`.

  CURVE_TYPE type () const { return this->_type; }         ///< Returns type of curve.
  
  real_type* params () const { return this->_params; }     ///< Returns the parameters of the curve.

  BOTH virtual LEN_T l() const = 0;                        ///< Returns the length of the curve.
};

#endif //CURVE_HH
