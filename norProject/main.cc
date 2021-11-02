#include<iostream>
#include<dp.hh>
#include<IOUtils.hh>
#include<IOSettings.hh>

#include"generator.hh"

int main(int argc, char const *argv[]) {
  std::cout << __FILE__ << "\nWARNING did you remember to:\n - Change the function in this file?\n - Change the settings to the correct values?\n";
  std::vector<Configuration2> points(SET::nPoints);
  std::vector<bool> fixedAngles(SET::nPoints); 

  if (SET::generator){
    fixedAngles.front()=true; 
    fixedAngles.back()=true;    
    
    points=getPoints(SET::nPoints, SET::xSize, SET::ySize, &fixedAngles);
  }
  else {
    std::pair<std::vector<Configuration2>, std::vector<bool> > confAnglesPair=readConfigurationsFromFile(SET::filename, SET::type);
    points=confAnglesPair.first;
    fixedAngles=confAnglesPair.second;
  }
  
  std::cout << "\nStarting points" << std::endl;
  for (auto a : points){
    std::cout << a << std::endl;
  }

  std::pair<LEN_T, std::vector<Angle> >ret= DP::solveDP(points, fixedAngles, std::vector<real_type>(1, SET::Kmax), SET::discr, SET::nRef, SET::saveAngles);
  LEN_T len=ret.first;
  std::vector<Angle> vtheta=ret.second;

  std::cout << "Finished points" << std::endl;

  for (int i=0; i<points.size(); i++){
    std::cout << points[i] << std::endl;
  }
  
  SET::outfile="provaCC.asy";

  getMPMDInfo(points, SET::Kmax, &vtheta, len);
  drawSolution(points, SET::Kmax, SET::outfile);

  return 0;
}
