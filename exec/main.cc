#include <dubins.hh>
#include <dp.hh>
#include <timeperf.hh>
#include<tests.hh>

#include<iostream>
#include<fstream>
#include<cmath>
#include<random>
#include<vector>
#include<utility>
#include<iomanip>
#include<algorithm>

void PrintScientific1D(real_type d){
  if (d == 0)
  {
    printf ("%*d", 6, 0);
    return;
  }

  int exponent  = (int)floor(log10( fabs(d)));  // This will round down the exponent
  real_type base   = d * pow(10, -1.0*exponent);

  printf("%1.1lfe%+01d", base, exponent);
}

void PrintScientific2D(real_type d){
  if (d == 0)
  {
    printf ("%*d", 7, 0);
    return;
  }

  int exponent  = (int)floor(log10( fabs(d)));  // This will round down the exponent
  real_type base   = d * pow(10, -1.0*exponent);

  printf("%1.1lfe%+02d", base, exponent);
}


std::vector<Configuration2> example1 = {
  Configuration2(0,0,-2.0*M_PI/8.0),
  Configuration2(2,2,ANGLE::INVALID),
  Configuration2(6,-1,ANGLE::INVALID),
  Configuration2(8,1,2.0*M_PI/8.0)
};

std::vector<std::string> testsNames = { 
  "Kaya Example 1",
  "Kaya Example 2",
  "Kaya Example 3",
  "Kaya Example 4",
  "Omega",
  "Circuit"
}; 

std::vector<std::vector<Configuration2> > Tests = {
  kaya1, kaya2, kaya3, kaya4, omega, spa
};

std::vector<K_T> Ks = {3.0, 3.0, 5.0, 3.0, 3.0, 3.0};
std::vector<uint> discrs = {4, 16, 90, 360};
std::vector<uint> refins = {1, 2, 4, 8, 16};
std::vector<LEN_T> exampleLenghts={3.41557885807514871601142658619, 6.27803455030931356617429628386, 11.9162126542854860389297755319, 7.46756219733842652175326293218, 41.0725016438839318766440555919, 6988.66098639942993031581863761}; //the last length is SPA

std::string nameTest(std::string name, std::string add="", std::string conc=" "){
  if (add==""){
    return name;
  }
  else{
    return name+conc+add;
  }
}


int main (){
  cout << "C++" << endl;

  for (int testID=0; testID<Tests.size(); testID++){
    // if (testID!=3){continue;}
    real_type dLen=exampleLenghts[testID];

    std::vector<bool> fixedAngles;
    for (uint i=0; i<Tests[testID].size(); i++){
      if (i==0 || i==Tests[testID].size()-1) {
        fixedAngles.push_back(true);
      }
      else {
        fixedAngles.push_back(false);
      }
    }
    std::vector<real_type> curveParamV={Ks[testID]};
    real_type* curveParam=curveParamV.data();
    
    for (auto DISCR :  discrs){
      // if (DISCR!=4){continue;}
      for (auto r : refins){
        // if (r!=4){continue;}
        //std::cout << DISCR << " " << r << " ";
        TimePerf tp, tp1;
        std::vector<Configuration2>points=Tests[testID];

        tp.start();
        std::vector<real_type> vtheta=DP::solveDP(points, DISCR, fixedAngles, r, curveParam); 
        auto time1=tp.getTime();

        LEN_T ComLength=vtheta[0];
        vtheta.erase(vtheta.begin());

        LEN_T Length;
        for (unsigned int idjijij=points.size()-1; idjijij>0; idjijij--){
          points[idjijij-1].th(vtheta[idjijij-1]);
          points[idjijij].th(vtheta[idjijij]);
          Dubins c(points[idjijij-1], points[idjijij], Ks[testID]);
          // std::cout << c << std::endl;
          Length+=c.l();
        }
        // std::cout << Length << " " << exampleLenghts[testID] << std::endl;
        printf("%3d & %2d & ", DISCR, r);
        PrintScientific2D((ComLength-exampleLenghts[testID])*1000.0);
        printf(" & ");
        PrintScientific1D(time1);
        printf("\\\\\n");
      }
    }
    printf("\n\n\n\n");
  }
  
  return 0;
}
