#ifndef DP_HH
#define DP_HH

#include <settings.hh>
#include <utils.hh>
#include <typedefs.hh>
#include <configuration.hh>
#include <clothoidG1.hh>
#include <dubins.hh>

#include <AsyPlot.hh>
#include <DubinsUtils.hh>

#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>

namespace DP {
  namespace {
    class Cell { //TODO please change name
    private:
      Angle _th;
      LEN_T _l; //Length of the curve
      Cell* _next;
      int _i, _j;

    public:
      Cell() : _th(ANGLE::INVALID), _l(0), _next(NULL) {}

      Cell(Angle th, LEN_T l, Cell* next, int i=0, int j=0) : _th(th), _l(l),  _next(next), _i(i), _j(j) {}

      Angle th() const { return this->_th; }

      LEN_T l() const { return this->_l; }
      
      Cell* next() const { return this->_next; }

      int i() const { return this->_i; }
      int j() const { return this->_j; }

      Angle th(Angle th) {
        this->_th = th;
        return this->th();
      }

      LEN_T l(LEN_T l) {
        this->_l = l;
        return this->l();
      }

      Cell* next(Cell* next){
        this->_next = next;
        return this->next();
      }

      Cell copy(const Cell &d) {
        this->th(d.th());
        this->l(d.l());
        this->next(d.next());
        this->_i=d.i();
        this->_j=d.j();

        return *this;
      }

      Cell operator=(const Cell &d) {
        return copy(d);
      }

      std::stringstream to_string(bool pretty = false) const {
        std::stringstream out;
        if (pretty) {
          out << "th: " << this->th() << " l: " << this->l();
        } else {
          out << "<" << (Angle) (this->th() * 1.0) << ", " << (LEN_T) (this->l()) << " (" << this->_i << ", " << this->_j << ")" << ">";
        }
        return out;
      }

      friend std::ostream &operator<<(std::ostream &out, const Cell &data) {
        out << data.to_string().str();
        return out;
      }
    };

    //Finds the closest discriminant angle. 
    inline int closestDiscr(Angle th, int DISCR) {
      double app = th / (2 * m_pi) * DISCR;
      int a = (int) (app + 1.0e-5);
      int b = app;
      return (a == b ? b : a);
    }
    
  }//Anonymous namespace to hide information

  void solveDP    (std::vector<Configuration2<double> > points, int discr, int startFromBottom=0, int stopFromTop=0);
  void solveDPOld (std::vector<Configuration2<double> > points, int size, int discr, int startFromBottom=0, int stopFromTop=0);
} //namespace DP
#endif //DP_HH

















