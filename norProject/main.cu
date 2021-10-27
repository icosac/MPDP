#include<iostream>
#include<stdio.h>
#include<cmath>

#include<dp.cuh>

using namespace std;

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
  std::vector<real_type> vtheta= DP::solveDP(points, 4, fixedAngles, std::vector<double>{Kmax}, 2, true, 4);
  
  std::cout << "Finished points" << std::endl;
  for (auto a : points){
    std::cout << a << std::endl;
  }
  getMPMDInfo(points, Kmax, &vtheta);
  drawSolution(points, Kmax, "provaCU.asy");
  return 0;
}
