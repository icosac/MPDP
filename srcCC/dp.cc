#ifndef CUDA_ON
#include <dp.hh>

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

void guessInitialAngles(const int i, std::vector<Angle>& thPrev, std::vector<Angle>& thCur, std::vector<Configuration2>& points, const std::vector<bool> fixedAngles){
  thPrev.clear();
  thCur.clear();
  
  thPrev.reserve(5);
  thCur.reserve(5);

  // aligned on straight line
  double th = std::atan2(points[i].y()-points[i-1].y(), points[i].x()-points[i-1].x());
  thPrev.push_back(th);
  thCur.push_back(th);

  // aligned on circle
  std::vector<double> XC, YC;
  circles(points[i-1].x(), points[i-1].y(), points[i].x(), points[i].y(), 1./Kmax, XC, YC);
  for (int j=0; j<XC.size(); ++j) 
  {
    th = std::atan2(points[i-1].y()-YC[j], points[i-1].x()-XC[j]);
    thPrev.push_back(th+M_PI/2.);
    thPrev.push_back(th-M_PI/2.);
    th = std::atan2(points[i].y()-YC[j], points[i].x()-XC[j]);
    thCur.push_back(th+M_PI/2.);
    thCur.push_back(th-M_PI/2.);
  }
}

#define MATRIX DP::matrix
void setSamplingAngles(int n, const std::vector<bool> fixedAngles, const std::vector<Configuration2> points){
  MATRIX.clear();
  MATRIX.resize(points.size());

  Angle dtheta=2*M_PI/discr;
  for (int i=0; i<x.size(); ++i){
    if (!fixedAngles[i]){
      MATRIX.reserve(discr+10);
      for (int j=0; j<dicr; ++j){
        MATRIX[i].push_back(DP::Cell(dtheta*j));
      }
      if (i>=1){
        std::vector<DP::Cell> thPrev, thCur;
        guessInitialAngles(i, thPrev, thCur);
        MATRIX[i-1].inesert(MATRIX[i-1].end(), thPrev.begin(), thPrev.end());
        MATRIX[i].inesert(MATRIX[i].end(), thCur.begin(), thCur.end());
      }
    }
    else{
      MATRIX[i].reserve(1);
      MATRIX[i]={ DP::Cell(points[i].th()) };
    }
  }
}

void setSamplingAngles(const std::vector<Configuration2> points, const std::vector<bool> fixedAngles, double hrange, int hn){
  MATRIX.clear();
  MATRIX.resize(points.size());

  Angle dtheta=hrange/hn;
  for (int i=0; i<x.size(); ++i){
    if (!fixedAngles[i]){
      MATRIX.reserve(2*hn+11); // up to 10 "special" + one for thref + hn on both side of thref
      MATRIX[i].push_back(DP::Cell(points[i.th()]));
      for (int j=0; j<hn; ++j){
        MATRIX[i].push_back(DP::Cell(dtheta*j));
      }
      if (i>=1){
        std::vector<DP::Cell> thPrev, thCur;
        guessInitialAngles(i, thPrev, thCur);
        MATRIX[i-1].inesert(MATRIX[i-1].end(), thPrev.begin(), thPrev.end());
        MATRIX[i].inesert(MATRIX[i].end(), thCur.begin(), thCur.end());
      }
    }
    else{
      MATRIX[i].reserve(1);
      MATRIX[i]={ DP::Cell(points[i].th()) };
    }
  }
}

//Once the matrix is given, find the best angles
std::vector<Angle> bestAngles(int size, std::vector<Configuration2>* points=NULL){
  int bestIdx = -1;
  LEN_T bestL = std::numeric_limits<LEN_T>::max();
  //Find best path overall
  for (int i=0; i<MATRIX[0].size(); ++i){
    if (MATRIX[0][i].l()<bestL){
      bestL=cells[0][i].l();
      bestIdx=1;
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
      points[i].th(vtheta[i]);
    }
  }

  return vtheta;
}


void DP::solveDPInner (std::vector<Configuration2> points, real_type* params){
  for (int idx=(x.size()-1); idx>=1; --idx){ //Cycle between all points starting from the last one

    #pragma omp parallel for
    for (int i=0; i<MATRIX[idx-1].size(); ++i){ //Consider the previous angle
      int bestJ=-1;
      double bestL = std::numeric_limits<LEN_T>::max();

      for (int j=0; j<MATRIX[idx].size(); ++j){ //Consider the next angle
        
        //Compute Dubins
        Dubins dub=Dubins(points[idx-1], points[i], params);
        LEN_T curL=dub.l();
    
        if (idx<MAX_IDX){
          curL+=MATRIX[idx][j].l();
        }

        if (currL<bestL){
          bestL=curL;
          bestJ=j;
        }
      }
      MATRIX[idx-1][i].l(bestL);
      MATRIX[idx-1][i].next(bestJ);
    }
  }

  bestAngles(points.size(), points);
}


std::vector<Angle> DP::solveDP (std::vector<Configuration2> points, int discr, const std::vector<bool> fixedAngles, int nRefinements, real_type* params){
  cout << "solveDP" << endl;
  std::vector<Angle> bestA;
  setSamplingAngles(discr);
  solveDPInner(points, fixedAngles, params);
  for (int i=0; i<points.size(); ++i){
    bestA.push_back(points.th());
  }

  for (int i=0; i<nRefinements; ++i){
    hrange=hrange/n*1.5;
    setSamplingAngles(points, fixedAngles, hrange, discr/2);
    solveDPInner(points, fixedAngles, params);
    bestA.clear();
    for (int i=0; i<points.size(); ++i){
      bestA.push_back(points.th());
    }
  }
  
  return bestA;
}

#endif 
