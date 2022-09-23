#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include <iostream>
#include <sstream>
#include <ostream>
#include <string>


#include <typedefs.hh>

class Configuration2{
private:
  real_type _x, _y; ///<Coordinates
  Angle _th;        ///<Angle
  K_T _k;           ///<Curvature

public:
  /*!
   * Default constructor to return a `Configuration2` initialized to 0.
   */
  Configuration2() : _x(0), _y(0), _th(0), _k(ANGLE::FREE) {}
  /*!
   * Constructor to set the values of a `Configuration2`.
   * @param x The abscissa coordinate.
   * @param y The ordinate coordinate.
   * @param th The angle through the point.
   * @param k The curvature in that point. Default is 0.
   */
  Configuration2(real_type x, real_type y, Angle th, K_T k=0) : _x(x), _y(y), _th(th), _k(k) {}

  /*!
   * Returns the value of the abscissa.
   * @return The value of the abscissa in that point.
   */
  real_type x()  const { return this->_x;  }
  /*!
   * Returns the value of the ordinate.
   * @return The value of the ordinate in that point.
   */
  real_type y()  const { return this->_y;  }
  /*!
   * Returns the value of the angle.
   * @return The value of the angle in that point.
   */
  Angle th()  const { return this->_th; }
  /*!
   * Returns the value of the curvature.
   * @return The value of the curvature in that point.
   */
   K_T k()     const { return this->_k;  }

  /*!
   * Sets the value of the abscissa. Also returns the new value.
   * @param x The new value of the abscissa.
   * @return The new set value.
   */
  real_type x(const real_type x) { this->_x=((real_type)x); return this->x(); }
  /*!
   * Sets the value of the ordinate. Also returns the new value.
   * @param x The new value of the ordinate.
   * @return The new set value.
   */
  real_type y(const real_type y) { this->_y=((real_type)y); return this->y(); }
  /*!
   * Sets the value of the angle. Also returns the new value.
   * @param x The new value of the angle.
   * @return The new set value.
   */
  Angle th(const Angle th) { this->_th=th; return this->th(); }
  /*!
   * Sets the value of the curvature. Also returns the new value.
   * @param x The new value of the curvature.
   * @return The new set value.
   */
  K_T k(const real_type k) { this->_k=k;   return this->k();  }

  /*!
   * Function to deep copy one configuration onto another.
   * @param conf The `Configuration2` from which to copy.
   * @return The new copied `Configuration2`.
   */
  Configuration2& copy (const Configuration2& conf){
    this->x(conf.x());
    this->y(conf.y());
    this->th(conf.th());
    this->k(conf.k());
    return *this;
  }

  /*!
   * Override of the assignment (=) operator. It is overridden with the copy function.
   * @param conf The `Configuration2` to copy from.
   * @return The new copied `Configuration2`.
   */
  Configuration2& operator= (const Configuration2& conf){
    this->copy(conf);
    return *this;
  }

  /*!
   * Function to print the most essential info about `Configuration2`.
   * @param str An additional string to add at the beginning.
   * @return A `std::stringstream` object containing the data of `Configuration2`.
   */
  std::stringstream to_string (const std::string& str="") const {
    std::stringstream out;
    out << (str.empty() ? "" : str+" ") << "x: " << this->x() << "  y: " << this->y() << "  th: " << this->th();
    return out;
  }

  /*! This function overrides the << operator so to print with `std::cout` the most essential info about the `Configuration2`.
  		\param[in] out The out stream.
  		\param[in] data The configuration to print.
  		\returns An output stream to be printed.
  */
  friend std::ostream& operator<<(std::ostream &out, const Configuration2& data) {
    out << data.to_string().str();
    return out;
  }
};

#endif //CONFIGURATION_HH
