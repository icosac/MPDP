#ifndef DP_HH
#define DP_HH

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
      Cell() : _th(ANGLE::INVALID), _l(std::numeric_limits<LEN_T>::max()), _nextID(0) {}

      BOTH Cell(Angle th, LEN_T l, int nextID) 
        : 
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

  uint guessInitialAngles(std::vector<std::set<Angle> >& moreAngles, const std::vector<Configuration2<double> >& points, const std::vector<bool> fixedAngles,  const real_type K);
  
  /* Templated functions in dp.cut
   * std::vector<Angle> bestAngles(DP::Cell* matrix, int discr, int size);
   * std::vector<Angle> solveDP (std::vector<Configuration2<double> > points, int discr, const std::vector<bool> fixedAngles, real_type* params, bool guessInitialAnglesVal=false);
   * __global__ void solveCell(DP::Cell* matrix, uint discr, uint size, std::vector<bool> fixedAngles, Configuration2<real_type>& c0, Configuration2<real_type>& c1, int& bestK, LEN_T& bestL, Angle& bestA, Angle& a00, Angle& a01, real_type* params, int i);
   */ 

  #include<dp.cut>

  template<class CURVE=Dubins<double> >
  __global__ void dubinsWrapper(Configuration2<real_type>c0, Configuration2<real_type> c1, double Kmax, double* L){
    CURVE c(c0, c1, Kmax);
    //printf("%.17f\n", c.l());
    L[0]+=c.l();
  }

  template<class CURVE>
  void solveDP(std::vector<Configuration2<real_type> >& points, int discr, const std::vector<bool> fixedAngles, std::vector<real_type> params, short type=2, bool guessInitialAnglesVal=false, uint nIter=1, uint threads=128, Angle _fullAngle=m_2pi){
    Angle fullAngle=_fullAngle;
    std::vector<Angle> angles; 
    //Passing the functions as pointers doesn't work for reasons I don't know
    //std::vector<Angle>(*func)(std::vector<Configuration2<real_type> > points, uint discr, const std::vector<bool> fixedAngles, std::vector<real_type> params, Angle fullAngle, bool halveDiscr, bool guessInitialAnglesVal)=NULL;
    for(uint i=0; i<nIter+1; ++i){
      //std::cout << "Refinement: " << i << std::endl;
      //std::cout << std::endl;
      switch(type){
        case 0: {
          angles=DP::solveDPFirstVersion<CURVE>(points, discr, fixedAngles, params, fullAngle, (i==0 ? false : true), guessInitialAnglesVal, threads);
          break;
        }
        case 1:{
          angles=DP::solveDPMatrixAllocator<CURVE>(points, discr, fixedAngles, params, fullAngle, (i==0 ? false : true), guessInitialAnglesVal, threads, i);
          break;
        }
        case 2: default:{
          angles=DP::solveDPAllIn1<CURVE>(points, discr, fixedAngles, params, fullAngle, (i==0 ? false : true), guessInitialAnglesVal, threads, i);
        }
      }

      for (uint j=0; j<angles.size(); j++){
        if (!fixedAngles[j]){
          points[j].th(angles[j]);
        }
      }
      //std::cout << "< ";
      //for (auto v : angles){
      //std::cout << std::setw(20) << std::setprecision(17) << mod2pi(v) << " ";
      //}
      //std::cout << ">" << std::endl;
      //
      //LEN_T* Length;
      //cudaMallocManaged(&Length, sizeof(LEN_T));
      //for (unsigned int h=points.size()-1; h>0; ){
      //  dubinsWrapper<<<1,1>>>(points[h-1], points[h], params[0], Length);
      //  cudaDeviceSynchronize();
      //  //std::cout << std::setw(20) << std::setprecision(17) << c.l() << std::endl;
      //  //Length+=c.l();
      //  h--;
      //}
      //std::cout << "Length: " << std::setw(20) << std::setprecision(17) << Length[0] << " " << std::endl; // setprecision(20) << (ABS<real_type>(Length, dLen)) << endl;
      //cudaFree(Length);

      if (i==0){
        fullAngle=fullAngle/(discr)*1.5;
        discr++; //This is because, yes.
      }
      else{
        fullAngle=fullAngle/(discr-1)*1.5;
      }
    }
  }

} //namespace DP


#endif //DP_HH

















