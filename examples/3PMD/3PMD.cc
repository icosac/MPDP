/**
 * @file 3PMD.cpp
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License Agero 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Main file for the 3 Point Markov-Dubins Problem.
 */

#include <3PMD.hh>

int main3PMDBruteForce(){
  // LRL-RLR case
  //Configuration2 Pi (0, 0, m_pi);
  //Configuration2 Pf (1, 1, m_pi);
  //Configuration2 Pm (0.5, 0.5, 0);

  // RSR-RLR case
  //Configuration2 Pi(-1, 0, m_pi / 3);
  //Configuration2 Pf(1, 0, -m_pi / 2+m_pi/8);
  //Configuration2 Pm(0.5, 0.5, 0);

  // Casse 2
  Configuration2 Pi(-1, 0, m_pi / 3);
  Configuration2 Pf(1, 0,  m_pi / 2+m_pi/8);
  Configuration2 Pm(0.5, 0.5, 0);

  K_T kmax =1;
  std::cout << "Test........." << std::endl;
  double thInc = 0.001;
  LEN_T bestLen = std::numeric_limits<LEN_T>::infinity();
  Angle bestAngle = 0.0;
  int n = 360;sud

  TimePerf time;
  time.start();
  int npts = 1000;
  std::string bestMan = "";
  double L0, L1, L2, L3, L4, L5;

  for (int j = 0; j < npts; ++j) {
    //for (double thm = -0*2 * m_pi; thm < 2* m_pi; thm+=thInc){
    //Pi.x(j);
    /*for (int i = 0; i < n; ++i) {
      double thm = 2.0*m_pi / n * i;
      Pm.th(thm);
      LEN_T currLen = Dubins(Pi, Pm, { kmax }).l() + Dubins(Pm, Pf, { kmax }).l();
          printf("%.2f, ", currLen);
      if (bestLen > currLen) {
        bestLen = currLen;
        bestAngle = thm;
      }
    }
      */


    for (int i = 0; i < n; ++i) {
      double thm = 2.0 * m_pi / n * i;
      Pm.th(thm);
      Dubins d1(Pi, Pm, { kmax });
      Dubins d2(Pm, Pf, { kmax });
      LEN_T currLen = d1.l() + d2.l();
      //printf("%.2f, ", currLen);
      if (bestLen > currLen) {
        bestLen = currLen;
        bestAngle = thm;
        bestMan = d1.man_to_string() + " " + d2.man_to_string();
        L0 = d1.s1();
        L1 = d1.s2();
        L2 = d1.s3();
        L3 = d2.s1();
        L4 = d2.s2();
        L5 = d2.s3();
      }
    }
  }


  printf("Shortest path with angle %.8f (%.2fpi) and total length %.8f with man %s\n", bestAngle, (bestAngle/ m_pi), bestLen, bestMan.c_str());
  std::cout << "L0 = " << L0 << std::endl;
  std::cout << "L1 = " << L1 << std::endl;
  std::cout << "L2 = " << L2 << std::endl;
  std::cout << "L3 = " << L3 << std::endl;
  std::cout << "L4 = " << L4 << std::endl;
  std::cout << "L5 = " << L5 << std::endl;
  std::cout << "Took: " << time.getTime()/npts << "ms" << std::endl;


  // plot
  std::ofstream file("3pointDubins.asy");
  file << "import graph; \n include \"clothoidLib.asylib\";\n size(8cm, 8cm); " << std::endl;

  // Solve the Reed-Shepp considering all the possible cases since curveParam only contains the curvature
  Pm.th(bestAngle);
  std::vector<real_type> curveParam = { kmax };
  Dubins myRS1(Pi, Pm, curveParam);
  Dubins myRS2(Pm, Pf, curveParam);

  std::string str = "";
  file << "path p;" << std::endl;

  Configuration2 c = myRS1.ci()[0];
  for (int ii = 1; ii < 4; ++ii) {
    str = "p = clothoidPoints((" + std::to_string(c.x()) + "," + std::to_string(c.y()) + "), " + std::to_string(c.th())
          + "," + std::to_string(myRS1.k(ii)) + ", 0, " + std::to_string(myRS1.L(ii)) + ");";
    std::cout << " k= " << myRS1.k(ii) << "     kmax = " << myRS1.kmax() << std::endl;
    file << str << std::endl;
    //plot.verbatim(str);
    str = "royalblue";

    //plot.verbatim("draw(p," + str + ");");

    file << "draw(p," << str << ");" << std::endl;

    //plot.dot(myRS.getX()[ii], myRS.getY()[ii], "red");
    file << "dot((" << c.x() << "," << c.y() << "), red);" << std::endl;

    c = circleLine(myRS1.L(ii), myRS1.k(ii), c);
  }
  file << "dot((" << myRS1.ci()->x() << "," << myRS1.ci()->y() << "), black);" << std::endl;
  file << "dot((" << myRS1.cf()->x() << "," << myRS1.cf()->y() << "), purple+3bp);" << std::endl;

  c = myRS2.ci()[0];
  for (int ii = 1; ii < 4; ++ii) {
    str = "p = clothoidPoints((" + std::to_string(c.x()) + "," + std::to_string(c.y()) + "), " + std::to_string(c.th())
          + "," + std::to_string(myRS2.k(ii)) + ", 0, " + std::to_string(myRS2.L(ii)) + ");";
    std::cout << " k= " << myRS2.k(ii) << "     kmax = " << myRS2.kmax() << std::endl;
    file << str << std::endl;
    //plot.verbatim(str);
    str = "red";

    //plot.verbatim("draw(p," + str + ");");

    file << "draw(p," << str << ");" << std::endl;

    //plot.dot(myRS.getX()[ii], myRS.getY()[ii], "red");
    file << "dot((" << c.x() << "," << c.y() << "), red);" << std::endl;

    c = circleLine(myRS2.L(ii), myRS2.k(ii), c);
  }
  file << "dot((" << myRS2.ci()->x() << "," << myRS2.ci()->y() << "), black);" << std::endl;
  file << "dot((" << myRS2.cf()->x() << "," << myRS2.cf()->y() << "), purple+3bp);" << std::endl;


  system("pause");

  return 0;
}

int main3PMDBruteForceWithPlot() {
  Configuration2 Pi(-1, 0, m_pi / 2);
  Configuration2 Pf(1, 1, 0);
  Configuration2 Pm(0.2, 0, m_pi / 2 + m_pi / 8);
  K_T kmax = 1;

  double thInc = 0.001; //std::numeric_limits<double>::epsilon();
  LEN_T bestLen = std::numeric_limits<LEN_T>::infinity();
  Angle bestAngle = 0.0;
  int n = 360;

  TimePerf time;
  time.start();
  int npts = 1;
  std::string bestMan = "";
  double L0, L1, L2, L3, L4, L5;

  for (int j = 0; j < npts; ++j) {
    //for (double thm = -0*2 * m_pi; thm < 2* m_pi; thm+=thInc){
    Pi.x(j);
    /*for (int i = 0; i < n; ++i) {
        double thm = 2.0*m_pi / n * i;
        Pm.th(thm);
        LEN_T currLen = Dubins(Pi, Pm, { kmax }).l() + Dubins(Pm, Pf, { kmax }).l();
        printf("%.2f, ", currLen);
        if (bestLen > currLen) {
            bestLen = currLen;
            bestAngle = thm;
        }
    }*/

    for (int i = 0; i < n; ++i) {
      double thm = 2.0 * m_pi / n * i;
      Pm.th(thm);
      Dubins d1(Pi, Pm, { kmax });
      Dubins d2(Pm, Pf, { kmax });
      LEN_T currLen = d1.l() + d2.l();
      //printf("%.2f, ", currLen);
      if (bestLen > currLen) {
        bestLen = currLen;
        bestAngle = thm;
        bestMan = d1.man_to_string() + " " + d2.man_to_string();
        L0 = d1.s1();
        L1 = d1.s2();
        L2 = d1.s3();
        L3 = d2.s1();
        L4 = d2.s2();
        L5 = d2.s3();
      }
    }
  }


  printf("Shortest path with angle %.8f (%.2fpi) and total length %.8f with man %s\n", bestAngle, (bestAngle / m_pi), bestLen, bestMan.c_str());
  std::cout << "L0 = " << L0 << std::endl;
  std::cout << "L1 = " << L1 << std::endl;
  std::cout << "L2 = " << L2 << std::endl;
  std::cout << "L3 = " << L3 << std::endl;
  std::cout << "L4 = " << L4 << std::endl;
  std::cout << "L5 = " << L5 << std::endl;
  std::cout << "Took: " << time.getTime() / npts << "ms" << std::endl;


  // plot
  std::ofstream file("3pointDubins.asy");
  file << "import graph; \n include \"clothoidLib.asylib\";\n size(8cm, 8cm); " << std::endl;

  // Solve the Reed-Shepp considering all the possible cases since curveParam only contains the curvature
  Pm.th(bestAngle);
  std::vector<real_type> curveParam = { kmax };
  Dubins myRS1(Pi, Pm, curveParam);
  Dubins myRS2(Pm, Pf, curveParam);

  std::string str = "";
  file << "path p;" << std::endl;

  Configuration2 c = myRS1.ci()[0];
  for (int ii = 1; ii < 4; ++ii) {
    str = "p = clothoidPoints((" + std::to_string(c.x()) + "," + std::to_string(c.y()) + "), " + std::to_string(c.th())
          + "," + std::to_string(myRS1.k(ii)) + ", 0, " + std::to_string(myRS1.L(ii)) + ");";
    std::cout << " k= " << myRS1.k(ii) << "     kmax = " << myRS1.kmax() << std::endl;
    file << str << std::endl;
    //plot.verbatim(str);
    str = "royalblue";

    //plot.verbatim("draw(p," + str + ");");

    file << "draw(p," << str << ");" << std::endl;

    //plot.dot(myRS.getX()[ii], myRS.getY()[ii], "red");
    file << "dot((" << c.x() << "," << c.y() << "), red);" << std::endl;

    c = circleLine(myRS1.L(ii), myRS1.k(ii), c);
  }
  file << "dot((" << myRS1.ci()->x() << "," << myRS1.ci()->y() << "), black);" << std::endl;
  file << "dot((" << myRS1.cf()->x() << "," << myRS1.cf()->y() << "), purple+3bp);" << std::endl;

  c = myRS2.ci()[0];
  for (int ii = 1; ii < 4; ++ii) {
    str = "p = clothoidPoints((" + std::to_string(c.x()) + "," + std::to_string(c.y()) + "), " + std::to_string(c.th())
          + "," + std::to_string(myRS2.k(ii)) + ", 0, " + std::to_string(myRS2.L(ii)) + ");";
    std::cout << " k= " << myRS2.k(ii) << "     kmax = " << myRS2.kmax() << std::endl;
    file << str << std::endl;
    //plot.verbatim(str);
    str = "red";

    //plot.verbatim("draw(p," + str + ");");

    file << "draw(p," << str << ");" << std::endl;

    //plot.dot(myRS.getX()[ii], myRS.getY()[ii], "red");
    file << "dot((" << c.x() << "," << c.y() << "), red);" << std::endl;

    c = circleLine(myRS2.L(ii), myRS2.k(ii), c);
  }
  file << "dot((" << myRS2.ci()->x() << "," << myRS2.ci()->y() << "), black);" << std::endl;
  file << "dot((" << myRS2.cf()->x() << "," << myRS2.cf()->y() << "), purple+3bp);" << std::endl;


  system("pause");
  return 0;
}

void main3PDP(){
  Configuration2 pi(0                 , 0                 , m_pi/3.0         );
  Configuration2 pm(10                , 5                 , 0                );
  Configuration2 pf(15                , 20                , m_pi/6.0         );

  K_T kmax = 1.0;

//  Dubins dub1 = Dubins(pi, pm, {kmax});
//  std::cout << std::endl << std::endl;
//  Dubins dub2 = Dubins(pm, pf, {kmax});
//  std::cout << std::endl << std::endl;
//
//  std::ofstream file("Dubins.asy");
//  initAsyFile(file);
//  dub1.draw(file, "P_i");
//  dub2.draw(file, "P_m");
//  file.close();

  std::vector<bool> fixedAngles = {true, false, true};
  std::vector<Configuration2> points = {pi, pm, pf};
  std::vector<double> curveParam = {kmax};
  TimePerf time1;
  time1.start();
  std::pair<LEN_T, std::vector<Angle> >ret=DP().solveDP(points, fixedAngles, curveParam, 360, 1);
  std::cout << "ms: " << time1.getTime() << std::endl;
//  std::cout << std::setprecision(12) << "Dub1: " << dub1.man_to_string() << " " << dub1.l() << std::endl;
//  std::cout << std::setprecision(12) << "Dub2: " << dub2.man_to_string() << " " << dub2.l() << std::endl;
//  std::cout << std::setprecision(12) << "Total length " << (dub1.l()+dub2.l()) << std::endl;
  std::cout << std::endl << std::endl << "MPDP len: " << ret.first << std::endl;
  std::cout << "#angles: " << ret.second.size() << std::endl;
  for (auto angle : ret.second){
    std::cout << std::setprecision(12) << angle << " ";
  }
  std::cout << std::endl;
  pm.th(ret.second[1]);
  std::cout << pm.th() << std::endl;
  Dubins curve1 = Dubins(pi, pm, {kmax});
  std::cout << std::setprecision(12) << "Curve1: " << curve1.man_to_string() << " " << curve1.l() << std::endl;
  std::cout << std::endl;
  Dubins curve2 = Dubins(pm, pf, {kmax});
  std::cout << std::setprecision(12) << "Curve2: " << curve2.man_to_string() << " " << curve2.l() << std::endl;
  std::cout << std::endl;
  std::cout << "Sum: " << (curve1.l()+curve2.l()) << std::endl;

  std::ofstream file1("Dubins1.asy");
  initAsyFile(file1);
  curve1.draw(file1, "P_i");
  curve2.draw(file1, "P_m");
  file1.close();

//  std::cout << "BRUTE FORCE" << std::endl;
//
//  LEN_T bestLen = std::numeric_limits<LEN_T>::infinity();
//  std::string bestMan = "";
//  int DISCR = 360;
//  Angle ang = 0.0, bestAngle = 0.0;
//  for (int i=0; i<DISCR; i++){
//    Dubins dub1 = Dubins(pi, Configuration2(pm.x(), pm.y(), ang), {kmax});
//    Dubins dub2 = Dubins(Configuration2(pm.x(), pm.y(), ang), pf, {kmax});
//    LEN_T currLen = dub1.l() + dub2.l();
//    if (bestLen > currLen) {
//      bestLen = currLen;
//      bestMan = dub1.man_to_string() + " " + dub2.man_to_string();
//      bestAngle = ang;
//    }
//    ang += 2.0*m_pi/DISCR;
//  }
//  std::cout << "Shortest path with angle " << bestAngle << " and total length " << bestLen << " given man: " << bestMan << std::endl;
}

void generateDataset3PDPCircle(int argc, char** argv){
  int kmax_min = 1;
  int kmax_max = 1;
  int k_discr = 1;
  int angle_discr = 5;

  if (argc == 4) {
    kmax_min = std::stoi(argv[1]);
    kmax_max = std::stoi(argv[2]);
    k_discr  = std::stoi(argv[3]);
  }
  else if (argc == 5) {
    kmax_min = std::stoi(argv[1]);
    kmax_max = std::stoi(argv[2]);
    k_discr  = std::stoi(argv[3]);
    angle_discr = std::stoi(argv[4]);
  }

  double theta_i;
  double theta_f;
  double alpha_m;
  double alpha_f;
  double kmax;

  // Open file named 3PDS.csv
  std::string filename_base = "3PDS_Circle" + std::to_string(angle_discr) + "_" + std::to_string(kmax_min) + "_"  + std::to_string(kmax_max) + "_"  + std::to_string(k_discr);
  std::string filename = filename_base + ".csv";
  std::string filename_log = filename_base + ".log";

  uint64_t tot_counter = angle_discr*angle_discr*angle_discr*angle_discr*k_discr;
  uint64_t counter = 0;
  uint64_t prev_counter = 0;
  uint64_t actual_counter = 0;

  std::cout << "Generating " << PrintScientificLargeInt(tot_counter) << " tests" << std::endl;

  std::cout << "Writing entries to " << filename << std::endl;
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cout << "Error opening db " << filename << std::endl;
    return;
  }

  std::cout << "Writing log to " << filename_log << std::endl;
  std::ofstream log_file(filename_log);
  std::streambuf* coutbuf = nullptr;
  if (!log_file.is_open()) {
    std::cout << "Error opening log file " << filename_log << std::endl;
    return;
  }
  else {
    coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(log_file.rdbuf());
  }

  std::cout << "Generating " << PrintScientificLargeInt(tot_counter) << " tests" << std::endl;

  kmax = kmax_max;
  for (int k = 0; k<k_discr; k++){
    TimePerf time1; time1.start();
    theta_i=m_pi;
    for (int g = 0; g < angle_discr; g++){
      theta_f=m_pi;
      for (int h = 0; h < angle_discr; h++) {
        alpha_m=m_pi;
        for (int i = 0; i < angle_discr; i++) {
          alpha_f=m_pi;
          for (int j = 0; j < angle_discr; j++) {
            Configuration2 pi = Configuration2(1, 0, theta_i);
            Configuration2 pm = Configuration2(cos(alpha_m), sin(alpha_m), 0);
            Configuration2 pf = Configuration2(cos(alpha_f), sin(alpha_f), theta_f);

            if (pm.x() != pi.x() && pm.y() != pi.y() && pm.x() != pf.x() && pm.y() != pf.y()){
              actual_counter ++;
               // Solve multipoint problem
               std::vector<Configuration2> points = {pi, pm, pf};
               std::vector<bool> fixedAngles = {true, false, true};
               std::vector<double> curveParam = {kmax};
               int discr = 90;
               int refinements = 4;
               TimePerf time;
               time.start();
               std::pair<LEN_T, std::vector<Angle> > ret = DP().solveDP(points, fixedAngles, curveParam, discr,
                                                                        refinements);
               if (ret.first == 0.0) {
                 std::cout << pi << std::endl << pm << std::endl << pf << std::endl;
                 throw std::runtime_error("Zero length");
               }
               auto dtime = time.getTime();
//            std::cout << "Took " << dtime << " ms to find multi-point" << std::endl;

               // Set angle for intermediate problem and compute the two Dubins
               pm.th(ret.second[1]);
               time.start();
               Dubins dub1 = Dubins(pi, pm, {kmax});
               Dubins dub2 = Dubins(pm, pf, {kmax});
               dtime = time.getTime();
//            std::cout << "Took " << dtime << " ms to find Dubins" << std::endl;
               LEN_T len = dub1.l() + dub2.l();

               // Get the manoeuvre combination, and if it's not in the 18 valid ones, search for an alternative
               std::string man_comb = dub1.man_to_string() + dub2.man_to_string();
               int id_man_comb = 19;
               time.start();
               auto search = P3DP_DICT.find(man_comb);
               if (search == P3DP_DICT.end()) {
                 for (auto man: P3DP_DICT) {
                   Dubins::D_TYPE dub1_man = std::get<1>(man.second);
                   Dubins::D_TYPE dub2_man = std::get<2>(man.second);
                   try {
                     Dubins dub1 = Dubins(pi, pm, {kmax}, dub1_man);
                     Dubins dub2 = Dubins(pm, pf, {kmax}, dub2_man);
                     if (std::abs(dub1.l() + dub2.l() - len) < 1e-8) {
                       id_man_comb = std::get<0>(man.second);
                       break;
                     }
                   }
                   catch (std::runtime_error &e) {
                     continue;
                   }
                 }
               } else {
                 id_man_comb = std::get<0>(search->second);
               }
               dtime = time.getTime();
//            std::cout << "Took " << dtime << " ms to find alternative" << std::endl;

               // Write data to file
               file << std::setprecision(5) << kmax << " " << theta_i << " " << theta_f << " " << alpha_m << " "
                    << alpha_f << " " << pm.th() << " " << id_man_comb << " " << len << std::endl;
             }

            // Update the angle
            alpha_f -= 2.0*m_pi/angle_discr;

            // Print time
            auto dtime1 = time1.getTime();
          	auto part = tot_counter > 100 ? tot_counter/100 : 1;
            if (counter % part == 0) {
              std::cout << 1.0 * counter / tot_counter * 100.0 << "% " << counter << " in " << dtime1 << "ms, avg " << (dtime1/(1.0*(counter-prev_counter))) << "ms" << std::endl;
              prev_counter = counter;
              time1.start();
            }
            counter ++;
          }
          alpha_m -= 2.0*m_pi/angle_discr;
        }
        theta_f -= 2.0*m_pi/angle_discr;
      }
      theta_i -= 2.0*m_pi/angle_discr;
    }
    kmax -= 1.0*(kmax_max-kmax_min)/k_discr;
  }

  std::cout << "Generated " << PrintScientificLargeInt(actual_counter) << " entries to " << filename << std::endl;

  if (coutbuf != nullptr){
    std::cout.rdbuf(coutbuf);
  }

  file.close();
}

void generateDataset3PDPRect(int argc, char** argv){
  int kmax_min = 1;
  int kmax_max = 1;
  int k_discr = 1;
  int angle_discr = 5;
  int xf_discr = 20;
  int xm_discr = 10;
  int ym_discr = 20;

  if (argc == 4) {
    kmax_min = std::stoi(argv[1]);
    kmax_max = std::stoi(argv[2]);
    k_discr  = std::stoi(argv[3]);
  }
  else if (argc == 8) {
    kmax_min = std::stoi(argv[1]);
    kmax_max = std::stoi(argv[2]);
    k_discr  = std::stoi(argv[3]);
    angle_discr = std::stoi(argv[4]);
    xf_discr = std::stoi(argv[5]);
    xm_discr = std::stoi(argv[6]);
    ym_discr = std::stoi(argv[7]);
  }

  // Open file named 3PDS.csv
  std::string filename_base = "3PDSRect_" + std::to_string(angle_discr) + "_" + std::to_string(kmax_min) + "_"  + std::to_string(kmax_max) + "_"  + std::to_string(k_discr);
  std::string filename = filename_base + ".csv";
  std::string filename_log = filename_base + ".log";

  uint64_t tot_counter = angle_discr*angle_discr*xf_discr*xm_discr*ym_discr*k_discr;
  uint64_t counter = 0;
  uint64_t prev_counter = 0;
  uint64_t actual_counter = 0;

  std::cout << "Generating " << PrintScientificLargeInt(tot_counter) << " tests" << std::endl;

  std::cout << "Writing entries to " << filename << std::endl;
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cout << "Error opening db " << filename << std::endl;
    return;
  }

  std::cout << "Writing log to " << filename_log << std::endl;
  std::ofstream log_file(filename_log);
  std::streambuf* coutbuf = nullptr;
  if (!log_file.is_open()) {
    std::cout << "Error opening log file " << filename_log << std::endl;
    return;
  }
  else {
    coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(log_file.rdbuf());
  }

  std::cout << "Generating " << PrintScientificLargeInt(tot_counter) << " tests" << std::endl;

  double kmax;

  std::ofstream points_file("3PDSRect_points.txt");

  double thi = -m_pi;
  double thf = -m_pi;
  double xf = 1;
  double xm = 0;
  double ym = 0;

  double dth = 2.0*m_pi/angle_discr;
  double dxf = 2.0/xf_discr;
  double dxm = 1.0/xm_discr;
  double dym = 2.0/ym_discr;
  double dkmax = 1.0*(kmax_max-kmax_min)/k_discr;

  std::cout << "dth: " << dth << std::endl;
  std::cout << "dxf: " << dxf << std::endl;
  std::cout << "dxm: " << dxm << std::endl;
  std::cout << "dym: " << dym << std::endl;
  std::cout << "dkmax: " << dkmax << std::endl;

  std::set<std::vector<double>> points_set;

  kmax = kmax_min;
  while(kmax <= (double)kmax_max)
  {
    TimePerf time1; time1.start();
    thi = -m_pi;
    while(thi < m_pi){
      thf = -m_pi;
      while(thf < m_pi)
      {
        xf = 1;
        while(xf > -1.0)
        {
          xm = 0;
          while(xm < 1)
          {
            ym = 0;
            while (ym <= 2)
            {
              if (ym == 0 && std::abs(xf-xm) < 1e-8) { ym += dym; continue; }
              Configuration2 pi = Configuration2(-1, 0, thi);
              Configuration2 pm = Configuration2(xm, ym, 0);
              Configuration2 pf = Configuration2(xf, 0, thf);

              points_file <<
                "(" << pi.x() << ", " << pi.y() << ", g)\n" <<
                "(" << pm.x() << ", " << pm.y() << ", r)\n" <<
                "(" << pf.x() << ", " << pf.y() << ", b)" << std::endl;
              if (pm.x() != pi.x() && pm.y() != pi.y() && pm.x() != pf.x() && pm.y() != pf.y()){
                actual_counter ++;
              // Solve multipoint problem
              std::vector<Configuration2> points = {pi, pm, pf};
              std::vector<bool> fixedAngles = {true, false, true};
              std::vector<double> curveParam = {kmax};
              int discr = 16;
              int refinements = 4;
              TimePerf time;
              time.start();
              std::pair<LEN_T, std::vector<Angle> > ret = DP().solveDP(points, fixedAngles, curveParam, discr, refinements);
                if (ret.first == 0.0) {
                  std::cout << pi << std::endl << pm << std::endl << pf << std::endl;
                  throw std::runtime_error("Zero length");
                }
                auto dtime = time.getTime();
                //            std::cout << "Took " << dtime << " ms to find multi-point" << std::endl;

                // Set angle for intermediate problem and compute the two Dubins
                pm.th(ret.second[1]);
                time.start();
                Dubins dub1 = Dubins(pi, pm, {kmax});
                Dubins dub2 = Dubins(pm, pf, {kmax});
                dtime = time.getTime();
                //            std::cout << "Took " << dtime << " ms to find Dubins" << std::endl;
                LEN_T len = dub1.l() + dub2.l();

                // Get the manoeuvre combination, and if it's not in the 18 valid ones, search for an alternative
                std::string man_comb = dub1.man_to_string() + dub2.man_to_string();
                int id_man_comb = 19;
                time.start();
                auto search = P3DP_DICT.find(man_comb);
                if (search == P3DP_DICT.end()) {
                  for (auto man: P3DP_DICT) {
                    Dubins::D_TYPE dub1_man = std::get<1>(man.second);
                    Dubins::D_TYPE dub2_man = std::get<2>(man.second);
                    try {
                      Dubins dub1 = Dubins(pi, pm, {kmax}, dub1_man);
                      Dubins dub2 = Dubins(pm, pf, {kmax}, dub2_man);
                      if (std::abs(dub1.l() + dub2.l() - len) < 1e-8) {
                        id_man_comb = std::get<0>(man.second);
                        break;
                      }
                    }
                    catch (std::runtime_error &e) {
                      continue;
                    }
                  }
                } else {
                  id_man_comb = std::get<0>(search->second);
                }
                dtime = time.getTime();
                //            std::cout << "Took " << dtime << " ms to find alternative" << std::endl;

                // Write data to file
                file << std::setprecision(5) << kmax << " " << thi << " " << xf << " " << thf << " "
                     << xm << " " << ym << " " << id_man_comb << " " << pm.th() << " " << len << std::endl;
              }

              // Print time
              auto dtime1 = time1.getTime();
              if (counter % (tot_counter/100000) == 0) {
                std::cout << 1.0 * counter / tot_counter * 100.0 << "% " << counter << " in " << dtime1 << "ms, avg " << (dtime1/(1.0*(counter-prev_counter))) << "ms" << std::endl;
                prev_counter = counter;
                time1.start();
              }
              counter ++;
              ym += dym;
            }
            xm += dxm;
          }
          xf -= dxf;
        }
        thf += dth;
      }
      thi += m_pi/angle_discr;
    }
    kmax += (dkmax > 0 ? dkmax : kmax_max);
  }

  std::cout << "Generated " << PrintScientificLargeInt(actual_counter) << " entries to " << filename << std::endl;

  if (coutbuf != nullptr){
    std::cout.rdbuf(coutbuf);
  }

  points_file.close();
  file.close();
}