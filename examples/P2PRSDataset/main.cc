// System includes
#include <random>
#include <fstream>
#include <iomanip>
#include <iostream>

// Library includes
#include <rs.hh>
#include <timeperf.hh>

int generateDatasetRS(){
  std::string datasetName = "DS_RS.csv";

  int n = 101; // discretisation for the angles
  int m = 51;  // discretisation for the curvature
  double th_min = -m_pi;
  double th_max = m_pi;
  double kappa_min = 0.1;
  double kappa_max = 9.0;

  std::vector<double> KAPPA (m);
  std::vector<double> THI (n);
  std::vector<double> THF (n);

  double v = kappa_min - (kappa_max-kappa_min)/(m-1);
  std::generate(KAPPA.begin(), KAPPA.end(), [&v, m, kappa_min, kappa_max]{ return v+=(kappa_max-kappa_min)/(m-1); });
  v = th_min - (th_max-th_min)/(n-1);
  std::generate(THI.begin(), THI.end(), [&v, n, th_min, th_max]{ return v+=(th_max-th_min)/(n-1); });
  v = th_min - (th_max-th_min)/(n-1);
  std::generate(THF.begin(), THF.end(), [&v, n, th_min, th_max]{ return v+=(th_max-th_min)/(n-1); });

  std::ofstream out (datasetName);

  int counter = 0;
  int secondCounter = 0;

  for (auto kappa : KAPPA){
    for (auto thi : THI){
      for (auto thf : THF){
        if ((thi >= 0) && (std::abs(thf) <= thi)) {
          secondCounter ++;
          if (secondCounter % 1000 == 0){
            std::cout << "Second counter: " << secondCounter << std::endl;
          }

          // generate only triangle
          Configuration2 ci(-1, 0, thi);
          Configuration2 cf(1, 0, thf);

          RS myRS = RS(ci, cf, { kappa });
          myRS.solve();


          out << std::setprecision(16) << std::fixed;
          out << thi << " " << thf << " " << kappa << " "
              << myRS.l() << " " << myRS.getNman() << " "
              << myRS.getManTypeStr() << " "
              << myRS.getNseg() << std::endl;
        }
      }
    }
  }
  out.close();

  std::cout << "Found " << counter << " alternatives" << std::endl;

  return 0;
}

int main(int argc, char** argv){
  return generateDatasetRS();
}