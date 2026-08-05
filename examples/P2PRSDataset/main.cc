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

  int n = 371; // discretisations for the angles
  int m = 51;  // discretisations for the curvature
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

  for (auto kappa : KAPPA){
    for (auto thi : THI){
      for (auto thf : THF){
        // generate only triangle
        if ((thi >= 0) && (std::abs(thf) <= thi)) {
          counter ++;
          if (counter % 1000 == 0){
            std::cout << "Second counter: " << counter << std::endl;
          }

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

  return 0;
}


int draw_rs(){
  // std::ofstream file ("rs_draw_failing_1.asy");
  // Configuration2 ci(-1, 0, 0.909844);
  // Configuration2 cf( 1, 0, -0.604671);
  // double kmax = 0.825669;
  // double pred_man = 27.0;

  // std::ofstream file ("rs_draw_failing_2.asy");
  // Configuration2 ci(-1, 0,  1.57809);
  // Configuration2 cf( 1, 0, -1.56737);
  // double kmax = 1.08065;
  // double pred_man = 31.0;

  // std::ofstream file ("rs_draw_failing_3.asy");
  // Configuration2 ci(-1, 0, 0.831105);
  // Configuration2 cf( 1, 0, -0.746854);
  // double kmax = 0.719206;
  // double pred_man = 6.0;

  // std::ofstream file ("rs_draw_failing_4.asy");
  // Configuration2 ci(-1, 0, 1.6906986599898821);
  // Configuration2 cf( 1, 0, 1.4288992721907325);
  // double kmax = 0.625;
  // double pred_man = 6.0;

  std::ofstream file ("rs_draw_failing_5.asy");
  Configuration2 ci(-1, 0, 0.831105);
  Configuration2 cf( 1, 0, -0.746854);
  double kmax = 0.719206;
  double pred_man = 6.0;

  std::cout << "===============================" << std::endl;
  std::cout << "Computing man with predicted maneuver: " << pred_man << std::endl;

  RS myRS = RS(ci, cf, { kmax, pred_man });
  myRS.solve();
  myRS.draw(file, 800, 800, false, false, true, true);

  std::cout << "RS length: " << myRS.l() << std::endl;
  std::cout << "Manoeuvre: " << myRS.getNman() << ", Type: " << myRS.getManTypeStr() << std::endl;
  std::cout << "Segments: " << myRS.getNseg() << std::endl;
  for (int i = 0; i < myRS.getNseg(); ++i){
    std::cout << "Segment " << i << ": x0 = " << myRS.getX()[i] << ", y0 = " << myRS.getY()[i] << ", L = " << myRS.getL()[i] << ", K = " << myRS.getK()[i] << ", D = " << myRS.getD()[i] << std::endl;
  }

  std::cout << "===============================" << std::endl;
  std::cout << "Computing optimal man" << std::endl;
  
  myRS = RS(ci, cf, { kmax });
  myRS.solve();
  myRS.draw(file, 800, 800, false, false, false, false, {"black+1bp", "purple+1bp"}, {"red", "orange"}, {"purple"});

  std::cout << "RS length: " << myRS.l() << std::endl;
  std::cout << "Manoeuvre: " << myRS.getNman() << ", Type: " << myRS.getManTypeStr() << std::endl;
  std::cout << "Segments: " << myRS.getNseg() << std::endl;
  for (int i = 0; i < myRS.getNseg(); ++i){
    std::cout << "Segment " << i << ": x0 = " << myRS.getX()[i] << ", y0 = " << myRS.getY()[i] << ", L = " << myRS.getL()[i] << ", K = " << myRS.getK()[i] << ", D = " << myRS.getD()[i] << std::endl;
  }
  std::cout << "==============================" << std::endl;

  file.close();

  return 0;
}

int main(int argc, char** argv){
  // return generateDatasetRS();
  return draw_rs();
}






