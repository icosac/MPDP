#ifndef DP_HH
#define DP_HH

#ifndef CURVE
#define CURVE Dubins
#endif

#include <settings.hh>
#include <utils.cuh>
#include <typedefs.hh>
#include <configuration.cuh>
#include <dubins.cuh>
#include <constants.cuh>

#include <iostream>
#include <set>
#include <cmath>
#include <vector>
#include <sstream>
#include <algorithm>

namespace DP {
  namespace {
    class Cell { //TODO please change name
    private:
      Angle _th;
      LEN_T _l; //Length of the curve
      int _nextID;

    public:
      real_type* _results;
      Cell() : _th(ANGLE::FREE), _l(std::numeric_limits<LEN_T>::max()), _nextID(0) {}

      BOTH Cell(Angle th, LEN_T l, int nextID) : 
          _th(th), 
          _l(l),  
          _nextID(nextID)
        {/*printf("nextID: %d %u\n" nextID);*/}

      BOTH Angle th()                 const { return this->_th; }
      BOTH LEN_T l()                  const { return this->_l; }   
      BOTH int next()                 const { return this->_nextID; }

      BOTH Angle th(Angle th) {
        this->_th = th;
        return this->th();
      }

      BOTH LEN_T l(LEN_T l) {
        this->_l = l;
        return this->l();
      }

      BOTH int next(int nextID){
        //printf("nextID in class %d %u\n", nextID, nextID);
        this->_nextID = nextID;
        return this->next();
      }

      BOTH Cell copy(const Cell &d) {
        this->th(d.th());
        this->l(d.l());
        this->next(d.next());
        
        return *this;
      }

      BOTH Cell operator=(const Cell &d) {
        return copy(d);
      }

      std::stringstream to_string(bool pretty = false) const {
        std::stringstream out;
        out << std::setw(20) << std::setprecision(17);
        if (pretty) {
          out << "th: " << this->th() << " l: " << this->l();
        } else {
          out << "<" << (Angle)(this->th()*1.0) << ", " << (LEN_T)(this->l()) << ">";
        }
        return out;
      }

      friend std::ostream &operator<<(std::ostream &out, const Cell &data) {
        out << data.to_string().str();
        return out;
      }
      
    };
  } //Anonymous namespace to hide information

  std::vector<real_type> solveDP(std::vector<Configuration2>& points, int discr, const std::vector<bool> fixedAngles, std::vector<real_type> params, short type=2, bool guessInitialAnglesVal=false, uint nIter=1, uint threads=128, Angle _fullAngle=m_2pi);

} //namespace DP


#endif //DP_HH

















