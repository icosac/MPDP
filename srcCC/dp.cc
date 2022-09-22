#ifndef CUDA_ON
#include <dp.hh>

void DP::guessInitialAngles(std::vector<Configuration2<double> >& points, const std::vector<bool> fixedAngles){
  for (uint i=0; i<points.size()-1; i++){
    if (!fixedAngles[i]){
      double diffX=points[i].x()-points[i+1].x();
      double diffY=points[i].y()-points[i+1].y();
      points[i].th(atan2(diffY, diffX));
    }
  }
}

#endif 
