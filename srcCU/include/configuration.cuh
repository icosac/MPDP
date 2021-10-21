#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include <iostream>
#include <sstream>
#include <ostream>
#include <string>

using namespace std;

#include <typedefs.hh>

class Configuration2{
private:
  real_type _x, _y;     ///<Coordinates
  Angle _th;     ///<Angle
  K_T _k;  ///<Curvature

public:
  Configuration2() : _x(0), _y(0), _th(0), _k(0) {}
  BOTH Configuration2(real_type x, real_type y, Angle th, K_T k=0) : _x(x), _y(y), _th(th), _k(k) {}

  BOTH real_type x() const { return this->_x; }
  BOTH real_type y() const { return this->_y; }
  BOTH Angle th() const { return this->_th; }
  K_T k() const { return this->_k; }

  BOTH real_type x(const real_type x) { this->_x=((real_type)x); return this->x(); }
  BOTH real_type y(const real_type y) { this->_y=((real_type)y); return this->y(); }

  BOTH Angle th(const Angle th) { this->_th=th; return this->th(); }
  K_T k(const real_type k) { this->_k=k;   return this->k();  }

  Configuration2 copy (Configuration2 conf){
    this->x(conf.x());
    this->y(conf.y());
    this->th(conf.th());
    this->k(conf.k());
    return *this;
  }

  Configuration2 operator= (Configuration2 conf){
    return copy(conf);
  }

  std::stringstream to_string (std::string str="") const {
    std::stringstream out;
    out << (str!="" ? "" : str+" ") << "x: " << this->x() << "  y: " << this->y() << "  th: " << this->th() << "  k: " << this->k();
    return out;
  }

  /*! This function overload the << operator so to print with `std::cout` the most essential info about the `Configuration2`.
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
