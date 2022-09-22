#ifndef CURVE_HH
#define CURVE_HH

#include <configuration.hh>

/*!
 * Type-parameter for the `Configuration2`, i.e., specifies the coordinate's type.
 * @tparam T1
 */
template<class T1>
class Curve{
private:
  Configuration2<T1> _ci; ///<Initial `Configuration`
  Configuration2<T1> _cf; ///<Final `Configuration`
  CURVE_TYPE _type;       ///<Type of curve
  real_type* _params;     ///<Parameters of curve

public:
  /*!
   * @brief Void constructor.
   */
  Curve() : _ci(), _cf(), _type(CURVE_TYPE::INVALID), _params(NULL) {}
  /*!
   * @brief Constructor to only set the type of the curve.
   */
  Curve(CURVE_TYPE type=CURVE_TYPE::INVALID) : _ci(), _cf(), _type(type), _params(NULL) {}

  /*!
   * @brief Constructor that takes two `Configuration2` and the type of the curve.
   * @param ci Initial configuration.
   * @param cf Final configuration.
   * @param type Type of the curve.
   * @param params The parameters of the curve, such as the curvature.
   */
  Curve(Configuration2<T1> ci, Configuration2<T1> cf, CURVE_TYPE type=CURVE_TYPE::INVALID, real_type* params=NULL) : _ci(ci), _cf(cf), _type(type), _params(params) {}

  Configuration2<T1>* ci() { return &(this->_ci); }   ///< Returns a pointer to the initial `Configuration2`.
  Configuration2<T1>* cf() { return &(this->_cf); }   ///< Returns a pointer to the final `Configuration2`.

  CURVE_TYPE type () const { return this->_type; }    ///< Returns type of curve.
  
  real_type* params () const { return this->_params; }         ///< Returns the parameters of the curve.

  virtual LEN_T l() const = 0;                        ///< Returns the length of the curve.
};

#endif //CURVE_HH
