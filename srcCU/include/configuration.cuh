#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include <iostream>
#include <sstream>
#include <ostream>
#include <string>

using namespace std;

#include <typedefs.hh>

template<class T1=float>
class Configuration2{
private:
  T1 _x, _y;     ///<Coordinates
  Angle _th;     ///<Angle
  K_T _k;  ///<Curvature

public:
  Configuration2() : _x(0), _y(0), _th(0), _k(0) {}
  BOTH Configuration2(T1 x, T1 y, Angle th, K_T k=0) : _x(x), _y(y), _th(th), _k(k) {}

  BOTH T1 x() const { return this->_x; }
  BOTH T1 y() const { return this->_y; }
  BOTH Angle th() const { return this->_th; }
  K_T k() const { return this->_k; }

  template<class T2>
  BOTH T1 x(const T2 x) { this->_x=((T1)x); return this->x(); }
  template<class T2>
  BOTH T1 y(const T2 y) { this->_y=((T1)y); return this->y(); }

  BOTH Angle th(const Angle th) { this->_th=th; return this->th(); }
  K_T k(const real_type k) { this->_k=k;   return this->k();  }

  template<class T2>
  Configuration2 copy (Configuration2<T2> conf){
    this->x(conf.x());
    this->y(conf.y());
    this->th(conf.th());
    this->k(conf.k());
    return *this;
  }

  template<class T2>
  Configuration2 operator= (Configuration2<T2> conf){
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
