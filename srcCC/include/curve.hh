#ifndef CURVE_HH
#define CURVE_HH

// Library includes
#include <configuration.hh>

// System includes
#include <vector>

enum class CURVE_TYPE { INVALID, CLOTHOID, DUBINS, DUBINS_ARC, RS }; ///< Possible types of CURVE

class Curve{
private:
  Configuration2 _ci; ///<Initial `Configuration`
  Configuration2 _cf; ///<Final `Configuration`
  CURVE_TYPE _type;   ///<Type of curve
  std::vector<real_type> _params; ///<Parameters of curve

public:

  /*!
   * @brief Void constructor.
   */
  Curve() : _ci(), _cf(), _type(CURVE_TYPE::INVALID), _params({}) {}
  /*!
   * @brief Constructor to only set the type of the curve.
   */
  Curve(CURVE_TYPE type=CURVE_TYPE::INVALID) : _ci(), _cf(), _type(type), _params({}) {}

  /*!
   * @brief Constructor that takes two `Configuration2` and the type of the curve.
   * @param ci Initial configuration.
   * @param cf Final configuration.
   * @param type Type of the curve.
   * @param params The parameters of the curve, such as the curvature.
   */
  Curve(Configuration2 ci, Configuration2 cf, CURVE_TYPE type=CURVE_TYPE::INVALID, std::vector<real_type> params = {}) : _ci(ci), _cf(cf), _type(type), _params(params) {}

  Configuration2* ci() { return &(this->_ci); }   ///< Returns a pointer to the initial `Configuration2`.
  Configuration2* cf() { return &(this->_cf); }   ///< Returns a pointer to the final `Configuration2`.

  CURVE_TYPE type () const { return this->_type; }      ///< Returns type of curve.
  
  std::vector<real_type> params () const { return this->_params; }  ///< Returns the parameters of the curve.

  virtual LEN_T l() const = 0;                          ///< Returns the length of the curve.

  virtual void solve() = 0;                             ///< Solves the curve depending on the type.

  virtual void draw() = 0;                              ///< Draws the curve.
};

#endif //CURVE_HH
