/**
 * @file main.cc
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Main file for the Dubins and Reed-Shepp paths computation.
 */

#include <dubins.hh>
#include <dp.hh>
#include <timeperf.hh>
#include <tests.hh>

#include<iostream>
#include<fstream>
#include<cmath>
#include<random>
#include<vector>
#include<utility>
#include<iomanip>
#include<algorithm>
#include<random>

using namespace mpdp;
using namespace mpdp::cpu;


int main(int argc, char** argv){
  TimePerf tp;
  tp.start();

  for (int i = 0; i < 1000; i++){
    for (int j = 1; j < 49; j++){
      RS rs (Configuration2(-1.0,  1.0, M_PI / 4.0),
            Configuration2( 1.0, -0.5, M_PI / 3.0),
            {0.5, static_cast<double>(j)});
      rs.solve();
    }
  }

  double elapsed = tp.getTime<std::micro>() / 1000.0;
  std::cout << "Elapsed time: " << elapsed << " us\n";

  // 0.70 kmax 1.90 xm 0.00 ym 0.10 th_i 0.79 th_f -0.79 th_m 6.06 true: 12 pred: 17 (0.7393) 11 (0.0938) 10 (0.0742) 16 (0.0684) 15 (0.0118)
  // std::vector<Configuration2> points = {
  //   Configuration2(-0.7, 0.0, 0.79), 
  //   Configuration2(0.0,  0.1, ANGLE::FREE), 
  //   Configuration2(0.7,  0.0, -0.79)
  // };

  // std::vector<bool> fixed_angles (points.size(), false);
  // fixed_angles[0] = true;
  // fixed_angles[fixed_angles.size()-1] = true;

  // float kmax = 1.90;

  // std::vector<double> params = {kmax};

  // DP dp;
  // auto res = dp.solveDP<Dubins>(points, fixed_angles, params, 90, 4, true);
  // std::cout << "Length: " << res.first << std::endl;
  // std::cout << "Angles: ";
  // for(auto a : res.second) std::cout << a << " ";
  // std::cout << std::endl;

  // dp.exportVisualizationData("dp_snapshot.json");

  // auto file = std::ofstream("dubins_path.asy");
  // for (size_t i=0; i<points.size() - 1; i++){
  //   Dubins dub (points[i], points[i+1], params);
  //   dub.draw(file, "dubins_path_" + std::to_string(i), 8, 8, false, false, i==0);
  //   std::cout << "Dubins " << i << " type: " << dub.dtype() << " " << (int)(dub.dtype()) << " " << dub.D_TYPE_STR[(int)(dub.dtype())] << std::endl;
  //   std::cout << "Dubins " << i << " length: " << dub.l() << std::endl;
  //   std::cout << "Dubins " << i << " s1: " << dub.L(1) << ", s2: " << dub.L(2) << ", s3: " << dub.L(3) << std::endl;
  //   std::cout << "Dubins " << i << " k1: " << dub.k(1) << ", k2: " << dub.k(2) << ", k3: " << dub.k(3) << std::endl;
  //   auto [first_intermediate, second_intermediate] = dub.get_intermediate_configurations();
  //   std::cout << "Dubins " << i << " first intermediate: " << first_intermediate << std::endl;
  //   std::cout << "Dubins " << i << " second intermediate: " << second_intermediate << std::endl;
  // }
  // file.close();


  // auto file2 = std::ofstream("dubins_marco.asy");
  // Dubins dub_marco1 (points[0], Configuration2(points[1].x(), points[1].y(), 1.0503), {kmax});
  // Dubins dub_marco2 (Configuration2(points[1].x(), points[1].y(), 1.0503), points[2], {kmax});
  // dub_marco1.draw(file2, "dubins_marco_1", 8, 8, false, false, true);
  // std::cout << "Dubins Marco 1 length: " << dub_marco1.l() << std::endl;
  // std::cout << "Dubins Marco 1 type: " << dub_marco1.dtype() << " " << (int)(dub_marco1.dtype()) << " " << dub_marco1.D_TYPE_STR[(int)(dub_marco1.dtype())] << std::endl;
  // std::cout << "Dubins Marco 1 s1: " << dub_marco1.L(1) << ", s2: " << dub_marco1.L(2) << ", s3: " << dub_marco1.L(3) << std::endl;
  // std::cout << "Dubins Marco 1 k1: " << dub_marco1.k(1) << ", k2: " << dub_marco1.k(2) << ", k3: " << dub_marco1.k(3) << std::endl;
  // Configuration2 dubins1_first_intermediate = dub_marco1.next_configuration(*(dub_marco1.ci()), dub_marco1.L(1), dub_marco1.k(1));
  // Configuration2 dubins1_second_intermediate = dub_marco1.next_configuration(dubins1_first_intermediate, dub_marco1.L(2), dub_marco1.k(2));
  // std::cout << "Dubins Marco 1 first intermediate: " << dubins1_first_intermediate << std::endl;
  // std::cout << "Dubins Marco 1 second intermediate: " << dubins1_second_intermediate << std::endl;

  // dub_marco2.draw(file2, "dubins_marco_2", 8, 8, false, false, false);
  // std::cout << "Dubins Marco 2 length: " << dub_marco2.l() << std::endl;
  // std::cout << "Dubins Marco 2 type: " << dub_marco2.dtype() << " " << (int)(dub_marco2.dtype()) << " " << dub_marco2.D_TYPE_STR[(int)(dub_marco2.dtype())] << std::endl;
  // std::cout << "Dubins Marco 2 s1: " << dub_marco2.L(1) << ", s2: " << dub_marco2.L(2) << ", s3: " << dub_marco2.L(3) << std::endl;
  // std::cout << "Dubins Marco 2 k1: " << dub_marco2.k(1) << ", k2: " << dub_marco2.k(2) << ", k3: " << dub_marco2.k(3) << std::endl;
  // auto [dubins2_first_intermediate, dubins2_second_intermediate] = dub_marco2.get_intermediate_configurations();
  // std::cout << "Dubins Marco 2 first intermediate: " << dubins2_first_intermediate << std::endl;
  // std::cout << "Dubins Marco 2 second intermediate: " << dubins2_second_intermediate << std::endl;
  // file2.close();

  // std::cout << "Dubins Marco combined length: " << dub_marco1.l() + dub_marco2.l() << std::endl;
  // std::cout << "Dubins Marco 1 type: " << dub_marco1.dtype() << " " << (int)(dub_marco1.dtype()) << " " << dub_marco1.D_TYPE_STR[(int)(dub_marco1.dtype())] << std::endl;
  // std::cout << "Dubins Marco 2 type: " << dub_marco2.dtype() << " " << (int)(dub_marco2.dtype()) << " " << dub_marco2.D_TYPE_STR[(int)(dub_marco2.dtype())] << std::endl;

  return 0;
}

