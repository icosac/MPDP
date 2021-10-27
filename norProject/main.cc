#include<iostream>
#include<stdio.h>
#include<cmath>

#include<dp.hh>

#include "IOUtils.hh"

int main(int argc, char const *argv[]) {
  std::vector<Configuration2> points=readConfigurationsFromFile("file1.txt", 1);
  
  std::vector<bool> fixedAngles; fixedAngles.resize(points.size(), false);
  fixedAngles.front()=true; fixedAngles.back()=true;
  
  std::cout << "Starting points" << std::endl;
  for (auto a : points){
    std::cout << a << std::endl;
  }

  real_type Kmax=0.5;
  std::vector<real_type> vtheta= DP::solveDP(points, 4, fixedAngles, 4, &Kmax);
  LEN_T len=vtheta[0];
  vtheta.erase(vtheta.begin());

  std::cout << "Finished points" << std::endl;

  for (int i=0; i<points.size(); i++){
    points[i].th(vtheta[i]);
    std::cout << points[i] << std::endl;
  }
  
  getMPMDInfo(points, Kmax, &vtheta, len);
  drawSolution(points, Kmax, "provaCC.asy");
  return 0;
}
