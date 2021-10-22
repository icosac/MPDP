#ifndef CUDA_ON
#include <dp.hh>

static K_T Kmax=DUBINS_DEFAULT_KMAX;
#define MATRIX DP::matrix

#ifdef DEBUG
#define printMatrix(type)                                               \
  if (type==0) {std::cout << "Angle table: " << std::endl;}             \
  else if (type==1) {std::cout << "Length table: " << std::endl;}       \
  for (int i=0; i<MATRIX.size(); i++){                                  \
    for (int j=0; j<MATRIX[i].size(); j++){                             \
      if (type==0){                                                     \
        std::cout << std::setprecision(4) << MATRIX[i][j].th() << "\t"; \
      }                                                                 \
      else if (type==1){                                                \
        std::cout << std::setprecision(4) << MATRIX[i][j].l() << "\t";  \
      }                                                                 \
    }                                                                   \
    std::cout << std::endl;                                             \
  }                                                                     \
  std::cout << std::endl << std::endl;                                               
#else 
#define printMatrix(type)
#endif

// returns (up to) two circles through two points, given the radius
static inline
void circles(double x1, double y1, double x2, double y2, double r, std::vector<double> & XC, std::vector<double> & YC) 
{
  double TOL = 1e-8;
  
  double q = std::hypot(x2-x1, y2-y1);
  double x3 = 0.5*(x1+x2);
  double y3 = 0.5*(y1+y2);

  double delta = r*r-q*q/4.;
    
  XC.clear();
  YC.clear();

  if (delta < -TOL) return;
  
  if (delta < TOL) 
  {
    XC.push_back(x3);
    YC.push_back(y3);
  }
  else
  {
    double deltaS = std::sqrt(delta);
    XC.push_back(x3 + deltaS*(y1-y2)/q);
    YC.push_back(y3 + deltaS*(x2-x1)/q);
    XC.push_back(x3 - deltaS*(y1-y2)/q);
    YC.push_back(y3 - deltaS*(x2-x1)/q);
  }
}

void guessInitialAngles(const int i, std::vector<DP::Cell>& thPrev, std::vector<DP::Cell>& thCur, const std::vector<Configuration2>& points){
  thPrev.clear();
  thCur.clear();
  
  thPrev.reserve(5);
  thCur.reserve(5);

  // aligned on straight line
  double th = std::atan2(points[i].y()-points[i-1].y(), points[i].x()-points[i-1].x());
  thPrev.push_back(DP::Cell(th));
  thCur.push_back(DP::Cell(th));

  // aligned on circle
  std::vector<double> XC, YC;
  circles(points[i-1].x(), points[i-1].y(), points[i].x(), points[i].y(), 1./Kmax, XC, YC);
  for (int j=0; j<XC.size(); ++j) 
  {
    th = std::atan2(points[i-1].y()-YC[j], points[i-1].x()-XC[j]);
    thPrev.push_back(DP::Cell(th+M_PI/2.));
    thPrev.push_back(DP::Cell(th-M_PI/2.));
    th = std::atan2(points[i].y()-YC[j], points[i].x()-XC[j]);
    thCur.push_back(DP::Cell(th+M_PI/2.));
    thCur.push_back(DP::Cell(th-M_PI/2.));
  }
}

void setSamplingAngles(int discr, const std::vector<bool> fixedAngles, const std::vector<Configuration2> points){
  MATRIX.clear();
  MATRIX.resize(points.size());

  Angle dtheta=2*M_PI/discr;
  for (int i=0; i<points.size(); ++i){
    MATRIX.reserve(discr+10);
    for (int j=0; j<discr; ++j){
      MATRIX[i].push_back(DP::Cell(dtheta*j));
    }
    if (i>0){
      std::vector<DP::Cell> thPrev, thCur;
      guessInitialAngles(i, thPrev, thCur, points);
      MATRIX[i-1].insert(MATRIX[i-1].end(), thPrev.begin(), thPrev.end());
      MATRIX[i].insert(MATRIX[i].end(), thCur.begin(), thCur.end());
    }
  }
  for (int i=0; i<points.size(); ++i){
    if (fixedAngles[i]){
      MATRIX[i].clear();
      MATRIX[i]={ DP::Cell(points[i].th()) };
    }
  }
}

void setSamplingAngles(const std::vector<Configuration2> points, const std::vector<bool> fixedAngles, double hrange, int hn){
  MATRIX.clear();
  MATRIX.resize(points.size());

  Angle dtheta=hrange/hn;
  for (int i=0; i<points.size(); ++i){
    MATRIX.reserve(2*hn+11); // up to 10 "special" + one for thref + hn on both side of thref
    MATRIX[i].push_back(DP::Cell(points[i].th()));
    for (int j=1; j<=hn; ++j){
      MATRIX[i].push_back(DP::Cell(dtheta*j+points[i].th()));
      MATRIX[i].push_back(DP::Cell(-dtheta*j+points[i].th()));
    }
    if (i>=1){
      std::vector<DP::Cell> thPrev, thCur;
      guessInitialAngles(i, thPrev, thCur, points);
      MATRIX[i-1].insert(MATRIX[i-1].end(), thPrev.begin(), thPrev.end());
      MATRIX[i].insert(MATRIX[i].end(), thCur.begin(), thCur.end());
    }
  }
  for (int i=0; i<points.size(); ++i){
    if (fixedAngles[i]){
      MATRIX[i].clear();
      MATRIX[i]={ DP::Cell(points[i].th()) };
    }
  }
}

//Once the matrix is given, find the best angles
std::vector<real_type> bestAngles(int size, std::vector<Configuration2>* points=NULL){
  int bestIdx = -1;
  LEN_T bestL = std::numeric_limits<LEN_T>::max();
  //Find best path overall
  for (int i=0; i<MATRIX[0].size(); ++i){
    if (MATRIX[0][i].l()<bestL){
      bestL=MATRIX[0][i].l();
      bestIdx=i;
    }
  }  

  std::vector<Angle> vtheta;
  vtheta.push_back(MATRIX[0][bestIdx].th());
  for (int i=0; i<size-1; ++i){
    int nIdx=MATRIX[i][bestIdx].next();
    vtheta.push_back(MATRIX[i+1][nIdx].th());
    bestIdx=nIdx;
  }

  if (points!=NULL){
    for (int i=0; i<size; ++i){
      (*points)[i].th(vtheta[i]);
    }
  }

  vtheta.insert(vtheta.begin(), bestL);

  return vtheta;
}

std::vector<real_type> solveDPInner (std::vector<Configuration2>& points, real_type* params){
  printMatrix(0)
  for (int idx=(points.size()-1); idx>0; --idx){ //Cycle between all points starting from the last one
    Configuration2* c0=&points[idx-1];
    Configuration2* c1=&points[idx];
    
    #pragma omp parallel for
    for (int i=0; i<MATRIX[idx-1].size(); ++i){ //Consider the previous angle
      int bestJ=-1;
      double bestL = std::numeric_limits<LEN_T>::max();

      for (int j=0; j<MATRIX[idx].size(); ++j){ //Consider the next angle
        //Compute Dubins
        Dubins dub=Dubins(c0->x(), c0->y(), MATRIX[idx-1][i].th(), c1->x(), c1->y(), MATRIX[idx][j].th(), params);
        LEN_T curL=dub.l();
        if (idx==(points.size()-1) && i==1){
          COUT(MATRIX[idx-1][i])  
          COUT(MATRIX[idx][j])  
          COUT(dub)
          COUT(curL)
        }
    
        if (idx<(points.size()-1)){
          curL+=MATRIX[idx][j].l();
        }

        if (curL<bestL){
          bestL=curL;
          bestJ=j;
        }
      }
      MATRIX[idx-1][i].l(bestL);
      MATRIX[idx-1][i].next(bestJ);
    }
  }
  printMatrix(1)
  return bestAngles(points.size(), &points);
}


std::vector<real_type> DP::solveDP (std::vector<Configuration2> points, int discr, const std::vector<bool> fixedAngles, int nRefinements, real_type* params){  
  MATRIX.clear();
  if (params!=NULL){
    Kmax=params[0];
  }
  std::vector<real_type> bestA;
  
  //First round
  setSamplingAngles(discr, fixedAngles, points);
  bestA=solveDPInner(points, params);

  //Other refinements
  double hrange=2.0*M_PI;
  for (int ref=0; ref<nRefinements; ++ref){
    COUT(ref) 
    hrange=hrange/discr*1.5;
    setSamplingAngles(points, fixedAngles, hrange, discr/2);
    printV(bestA)
    
    bestA.clear();
    bestA=solveDPInner(points, params);
  }
  
  return bestA;
}

#endif 
