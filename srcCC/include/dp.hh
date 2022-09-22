#ifndef DP_HH
#define DP_HH

#include <settings.hh>
#include <utils.hh>
#include <typedefs.hh>
#include <configuration.hh>
#include <clothoidG1.hh>
#include <dubins.hh>

#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
//#include <cmath>

namespace DP {
  namespace {
    class Cell { //TODO please change name
    private:
      Angle _th;
      LEN_T _l; //Length of the curve
      Cell* _next;
      int _i, _j, _id;

    public:
      Cell() : _th(ANGLE::INVALID), _l(std::numeric_limits<LEN_T>::max()), _next(NULL) {}

      Cell(Angle th, LEN_T l, Cell* next, int i=0, int j=0, int id=0) : _th(th), _l(l),  _next(next), _i(i), _j(j), _id(id) {}

      Angle th() const { return this->_th; }

      LEN_T l() const { return this->_l; }
      
      Cell* next() const { return this->_next; }

      int i() const { return this->_i; }
      int j() const { return this->_j; }
      int id() const { return this->_id; }

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
        this->_id=d.id();

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
          out << "<" << (Angle) (this->th() * 1.0) << ", " << (LEN_T) (this->l()) << " (" << this->_i << ", " << this->_j << ", " << this->_id << ")" << ">";
        }
        return out;
      }

      friend std::ostream &operator<<(std::ostream &out, const Cell &data) {
        out << data.to_string().str();
        return out;
      }
      
      operator real_type(){ //TODO remove this function when finished debugging
        real_type val=this->l();
        if (false) val=this->th();
        return (val>1000000 ? 999999 : val);
      }
    };
  } //Anonymous namespace to hide information

  void guessInitialAngles(std::vector<Configuration2<double> >& points, const std::vector<bool> fixedAngles);
  std::vector<Angle> bestAngles(DP::Cell* matrix, int discr, int size);

template<class CURVE>
void solveCell(DP::Cell* matrix, uint discr, uint size, std::vector<bool> fixedAngles, Configuration2<real_type>& c0, Configuration2<real_type>& c1, int& bestK, LEN_T& bestL, Angle& bestA, Angle& a00, Angle& a01, real_type* params, int i){
  for (int k=0; k<discr; k++){
    LEN_T currL=std::numeric_limits<LEN_T>::max();
    if (!fixedAngles[i]){ c1.th(a01+m_2pi*k/(discr*1.0)); } //If angle is fixed I don't have to change it
    CURVE c;
    try{
      c=CURVE(c0, c1, params); 
    } catch (runtime_error* e){}
    DP::Cell* next=(i==size-1 ? NULL : &matrix[k*size+(i+1)]);
    if (c.l()>0){
      currL=c.l();
      if (next!=NULL){
        currL+=next->l();
      }  
      if (currL<bestL || bestL==std::numeric_limits<LEN_T>::max()){
        bestL=currL;
        bestA=c1.th();
        bestK=k;
      }
    }
    if (fixedAngles[i]){ k=discr; } //If the angle is fixed I don't have to change it
  }
}

template<class CURVE>
std::vector<Angle> solveDP (std::vector<Configuration2<double> > points, int discr, const std::vector<bool> fixedAngles, real_type* params, bool guessInitialAnglesVal=false){
  cout << "solveDP" << endl;
  uint size=points.size();
  if (points.size()!=fixedAngles.size()){
    cerr << "Number of points and number of fixed angles are not the same: " << points.size() << "!=" << fixedAngles.size() << endl;
    return std::vector<Angle>();
  }
  if (guessInitialAnglesVal){
    DP::guessInitialAngles(points, fixedAngles);
  }

  DP::Cell* matrix=new DP::Cell[discr*size];

  for (uint i=size-1; i>0; i--){
    Configuration2<double>c0=points[i-1];
    Configuration2<double>c1=points[i];
    Angle a00=c0.th(), a01=c1.th();
    for (int j=0; j<discr; j++){
      Angle bestA=0.0;
      LEN_T bestL=std::numeric_limits<LEN_T>::max(); 
      int bestK=0;
      if (!fixedAngles[i-1]){ c0.th(a00+m_2pi*j/(discr*1.0)); } //If angle is fixed I don't have to change it
      DP::solveCell<CURVE>(matrix, discr, size, fixedAngles, c0, c1, bestK, bestL, bestA, a00, a01, params, i);
      if (bestL!=std::numeric_limits<LEN_T>::max()){
        Cell* next=(i==size-1? NULL : &matrix[bestK*size+(i+1)]);
        matrix[j*size+i]=Cell(bestA, bestL, next, i, j, j*size+i);
      }
      if (i==1){
        matrix[size*j]=Cell(c0.th(), bestL, &matrix[size*j+i], 0, j, size*j);
      }
      if (fixedAngles[i-1]){ j=discr; } //If the angle is fixed I don't have to change it
    }
  }
#ifdef DEBUG
  cout << "Printing " << endl;
  printVM(matrix, discr, size)
  //Retrieve angles
  cout << "Computing best angles" << endl;
#endif
  std::vector<Angle> bestA=DP::bestAngles(matrix, discr, size);
#ifdef DEBUG
  printV(bestA)
#endif
  
  double Length=0.0;
  for (unsigned int i=bestA.size()-1; i>0; i--){
    points[i].th(bestA[i]);
    points[i-1].th(bestA[i-1]);
    CURVE c(points[i-1], points[i], params);
    Length+=c.l();
  }
#ifdef DEBUG
  cout << "Length: " << setprecision(20) << Length << endl;

  cout << "Printing for Matlab" << endl;
  cout << "X=[";
  for (unsigned int i=0; i<points.size(); i++){ cout << points[i].x() << (i!=points.size()-1 ? ", " : "];\n"); }
  cout << "Y=[";
  for (unsigned int i=0; i<points.size(); i++){ cout << points[i].y() << (i!=points.size()-1 ? ", " : "];\n"); }
  cout << "th=[";
  for (unsigned int i=0; i<bestA.size(); i++){ cout << bestA[i] << (i!=bestA.size()-1 ? ", " : "];\n"); }
  cout << "KMAX: " << params[0] << endl;
#endif

  delete [] matrix;
  return bestA;
}
} //namespace DP

#include<dp.tt>

#endif //DP_HH

















