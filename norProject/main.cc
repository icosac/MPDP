#include<iostream>
#include<stdio.h>
#include<cmath>

#include<dp.hh>

using namespace std;

#include "IOUtils.hh"

int main(int argc, char const *argv[]) {
  std::vector<Configuration2> points=readConfigurationsFromFile("file1.txt", 1);
  
  std::vector<bool> fixedAngles; fixedAngles.resize(points.size(), false);
  fixedAngles.front()=true; fixedAngles.back()=true;
  
  for (auto a : points){
    std::cout << a << std::endl;
  }

  real_type Kmax=1;
  std::vector<real_type> vtheta= DP::solveDP(points, 4, fixedAngles, 4, &Kmax);
  LEN_T len=vtheta[0];
  vtheta.erase(vtheta.begin());

  for (auto a : points){
    std::cout << a << std::endl;
  }
  getMPMDInfo(points, Kmax, &vtheta, len);
  drawSolution(points, Kmax, "provaCC.asy");
  return 0;
}
