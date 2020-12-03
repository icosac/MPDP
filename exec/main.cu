#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include<stdlib.h>
#include<unistd.h>
using namespace std;

// #define DEBUG

#include<utils.cuh>
#include<dubins.cuh>
#include<dp.cuh>
#include<timeperf.hh>
#include<constants.cuh>

#include<tests.hh>

vector<Configuration2<real_type> > example1 = {
  Configuration2<real_type>(0,0,-2.0*M_PI/8.0),
  Configuration2<real_type>(2,2,ANGLE::INVALID),
  Configuration2<real_type>(6,-1,ANGLE::INVALID),
  Configuration2<real_type>(8,1,2.0*M_PI/8.0)
};

vector<std::string> testsNames = { 
  "Kaya Example 1",
  "Kaya Example 2",
  "Kaya Example 3",
  "Kaya Example 4",
  "Omega",
  "Circuit"
}; 

vector<vector<Configuration2<double> > > Tests = {
  kaya1, kaya2, kaya3, kaya4, omega, albert
};

vector<K_T> Ks = {3.0, 3.0, 5.0, 3.0, 3.0, 0.1};
vector<uint> discrs = {4, 120, 360, 720, 1440, 2880};

#define DISCR 1440

__global__ void kernel(){
  printf("ciao\n");
}

int main (){
  cout << "CUDA" << endl;
  cudaFree(0);

  int devicesCount;
  cudaGetDeviceCount(&devicesCount);
  cudaDeviceProp deviceProperties;
  cudaGetDeviceProperties(&deviceProperties, 0);
  printf("[%d] %s\n", 0, deviceProperties.name);

#if true 
  int testI=0;
  fstream json_out; json_out.open("tests.json", std::fstream::app);
  // std::cout << "\t\t        \tMatrix\t\tCol\tCol-Matrix" << std::endl;
  for (uint discr : discrs){
    cout << "Discr: " << discr << endl;
    for (uint j=0; j<Tests.size(); j++){
      std::vector<bool> fixedAngles;
      vector<Configuration2<double> > v=Tests[j];
      for (int i=0; i<v.size(); i++){
        if (i==0 || i==v.size()-1) {
          fixedAngles.push_back(true);
        }
        else {
          fixedAngles.push_back(false);
        }
      }
      std::vector<real_type> curveParamV={Ks[j]};
      real_type* curveParam=curveParamV.data();

      system((std::string("tegrastats --interval 50 --start --logfile ")+std::to_string(testI)+".log").c_str());
      sleep(2);
      
      TimePerf tp, tp1;
      
      tp.start();
      DP::solveDPMatrix<Dubins<double> >(v, discr, fixedAngles, curveParamV, false);
      auto time1=tp.getTime();
      Run r1(deviceProperties.name, discr, time1, testsNames[j]);
      r1.write(json_out);
      
      // tp1.start();
      // DP::solveDP<Dubins<double> >(v, discr, fixedAngles, curveParamV, false);
      // auto time2=tp1.getTime();
      // Run r2("Xavier", discr, time2, testsNames[j]);
      // r2.write(json_out);
      
      sleep(2);
      system("tegrastats --stop");
      testI++;
      // cout << "\tExample " << j+1 << std::setw(20) << std::setprecision(5) << time1 << "ms\t" << std::setw(20) << std::setprecision(5) <<  time2 << "ms\t" << std::setw(10) << (time2-time1) << "ms" << endl;
    }
  }
  json_out << "]}\n";
  json_out.close();
  
#else
  #define KAYA omega
  std::vector<bool> fixedAngles;
  for (int i=0; i<KAYA.size(); i++){
    if (i==0 || i==KAYA.size()-1) {
      fixedAngles.push_back(true);
    }
    else {
      fixedAngles.push_back(false);
    }
  }
  std::vector<real_type> curveParamV={2.0};
  real_type* curveParam=curveParamV.data();
  
  TimePerf tp, tp1;
  tp.start();
  DP::solveDPMatrix<Dubins<double> >(KAYA, DISCR, fixedAngles, curveParamV, false);
  auto time1=tp.getTime();
  tp1.start();
  DP::solveDP<Dubins<double> >(KAYA, DISCR, fixedAngles, curveParamV, false);
  auto time2=tp1.getTime();
  cout << "Elapsed: " << std::setw(10) << time1 << "ms\t" << std::setw(10) << time2 << "ms\t" << std::setw(10) << (time2-time1) << "ms" << endl;
#endif
  return 0;
}

