/**
 * @file 3PMD.cpp
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Main file for the 3 Point Markov-Dubins Problem.
 */

#include <3PMD.hh>
#include <cstring>

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
  int n = 360;

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
      LEN_T currLen = Dubins(Pi, Pm, kmax).l() + Dubins(Pm, Pf, kmax).l();
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
      Dubins d1(Pi, Pm, kmax);
      Dubins d2(Pm, Pf, kmax);
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

  #ifdef _WIN32
  system("pause");
  #endif


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
        LEN_T currLen = Dubins(Pi, Pm, kmax).l() + Dubins(Pm, Pf, kmax).l();
        printf("%.2f, ", currLen);
        if (bestLen > currLen) {
            bestLen = currLen;
            bestAngle = thm;
        }
    }*/

    for (int i = 0; i < n; ++i) {
      double thm = 2.0 * m_pi / n * i;
      Pm.th(thm);
      Dubins d1(Pi, Pm, kmax);
      Dubins d2(Pm, Pf, kmax);
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


  #ifdef _WIN32
  system("pause");
  #endif

  return 0;
}

void compute3Pman(std::string ThreePman){
  Configuration2 Pi(1, 0, 4/9*m_pi);
  Configuration2 Pm(cos(-2/9*m_pi), sin(-2/9*m_pi), 0.0);
  Configuration2 Pf(cos(-5/9*m_pi), sin(-5/9*m_pi), -m_pi/3);

  if (ThreePman!=""){
    Dubins::D_TYPE man1 = std::get<1>(P3DP_DICT.at(ThreePman));
    Dubins::D_TYPE man2 = std::get<2>(P3DP_DICT.at(ThreePman));
  }

  std::vector<double> curveParam = {2.0};

  TimePerf time1;
  time1.start();

  // Dubins dub1 = Dubins(Pi, Pm, curveParam);
  // Dubins dub2 = Dubins(Pm, Pf, curveParam);

  // std::cout << "Dub1: " << dub1 << std::endl;
  // std::cout << dub1.s1() << " " << dub1.s2() << " " << dub1.s3() << std::endl;
  // std::cout << "Dub2: " << dub2 << std::endl;
  // std::cout << dub2.s1() << " " << dub2.s2() << " " << dub2.s3() << std::endl;

  // std::cout << "Len: " << dub1.l() + dub2.l() << std::endl;
  // std::cout << "ms: " << time1.getTime() << std::endl;

  // Compute with DP
  std::vector<Configuration2> points = {Pi, Pm, Pf};
  std::vector<bool> fixedAngles = {true, false, true};

  time1.start();
  std::pair<LEN_T, std::vector<Angle> >ret=DP().solveDP(points, fixedAngles, curveParam, 360, 4);
  std::cout << "ms: " << time1.getTime() << std::endl;
  std::cout << "MPDP len: " << ret.first << std::endl;
  std::cout << "Angles" << std::endl;
  for (auto angle : ret.second){
    std::cout << angle << " ";
  }
}

double main3PDP(Configuration2& pi, Configuration2& pm, Configuration2& pf, K_T kmax){
  // Configuration2 pi(-1.0  , 0.0  ,  m_pi/2.0 );
  // Configuration2 pm( 0.25 , 0.75 ,  0.0      );
  // Configuration2 pf( 1.0  , 0.0  , -m_pi/2.0 );

  // K_T kmax = 2.0;

//  Dubins dub1 = Dubins(pi, pm, kmax);
//  std::cout << std::endl << std::endl;
//  Dubins dub2 = Dubins(pm, pf, kmax);
//  std::cout << std::endl << std::endl;
//
//  std::ofstream file("Dubins3PSquares.asy");
//  initAsyFile(file);
//  dub1.draw(file, "P_i");
//  dub2.draw(file, "P_m");
//  file.close();

  std::vector<bool> fixedAngles = {true, false, true};
  std::vector<Configuration2> points = {pi, pm, pf};
  std::vector<double> curveParam = { kmax };
  TimePerf time1;
  time1.start();
  std::pair<LEN_T, std::vector<Angle> >ret=DP().solveDP(points, fixedAngles, curveParam, 1440, 1);
  std::cout << "ms: " << time1.getTime() << std::endl;
  // std::cout << std::setprecision(12) << "Dub1: " << dub1.man_to_string() << " " << dub1.l() << std::endl;
  // std::cout << std::setprecision(12) << "Dub2: " << dub2.man_to_string() << " " << dub2.l() << std::endl;
  // std::cout << std::setprecision(12) << "Total length " << (dub1.l()+dub2.l()) << std::endl;
  std::cout << std::endl << std::endl << "MPDP len: " << ret.first << std::endl;
  std::cout << "#angles: " << ret.second.size() << std::endl;
  for (auto angle : ret.second){
    std::cout << std::setprecision(12) << angle << " ";
  }
  std::cout << std::endl;

  pm.th(ret.second[1]);
  std::cout << pm.th() << std::endl;
  Dubins curve1 = Dubins(pi, pm, kmax);
  std::cout << std::setprecision(12) << "Curve1: " << curve1.man_to_string() << " " << curve1.l() << std::endl;
  std::cout << std::setprecision(12) << "Curve1: " << curve1 << std::endl;
  std::cout << "s1: " << curve1.s1() << " " << "s2: " << curve1.s2() << " " << "s3: " << curve1.s3() << std::endl;
  std::cout << std::endl;

  Dubins curve2 = Dubins(pm, pf, kmax);
  std::cout << std::setprecision(12) << "Curve2: " << curve2.man_to_string() << " " << curve2.l() << std::endl;
  std::cout << std::setprecision(12) << "Curve2: " << curve2 << std::endl;
  std::cout << "s1: " << curve2.s1() << " " << "s2: " << curve2.s2() << " " << "s3: " << curve2.s3() << std::endl;
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
//    Dubins dub1 = Dubins(pi, Configuration2(pm.x(), pm.y(), ang), kmax);
//    Dubins dub2 = Dubins(Configuration2(pm.x(), pm.y(), ang), pf, kmax);
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


void main3PDPConfigurations(){
  Configuration2 pi( 2.0    , 0.0    , -3.1416 );
  Configuration2 pm( 3.8637 , 2.4289 ,  0.0      );
  Configuration2 pf(-2.0    , 0.0    , -4.1887  );

  // xi = 2., yi = 0., xm = 3.8637, ym = 2.4289, xf = -2., yf = 0., thi = -3.1416, thf = -4.1887, r = 1.

  K_T kmax = 1.0;

  main3PDP(pi, pm, pf, kmax);
}


void main3PDPCircle(){
  // double thi = 2.3562;
  // double thf = 1.5708;
  // double alpham = -0.7854;
  // double alphaf = 0.0;

  double kmax = 1.5;
  double thi = 1.88495559215;
  double thm = 0.357692328704;
  double thf = 0.628318530718;
  double alpham = -3.14159265359; 
  double alphaf = 1.88495559215;
  double len = 6.19613;

  Configuration2 pi(1, 0, thi);
  Configuration2 pm(cos(alpham), sin(alpham), thm);
  Configuration2 pf(cos(alphaf), sin(alphaf), thf);

	std::string D_TYPE_STR__[7] = {"INVALID", "LRL", "RLR", "LSL", "LSR", "RSL", "RSR"};
  for (int i=1; i<7; i++){
    Dubins::D_TYPE man1 = static_cast<Dubins::D_TYPE>(i);
    for (int j=1; j<7; j++){
      Dubins::D_TYPE man2 = static_cast<Dubins::D_TYPE>(j);
      std::cout << "Testing man1 " << D_TYPE_STR__[man1] << " man2 " << D_TYPE_STR__[man2] << std::endl;
      try{
        Dubins dub1_ = Dubins(pi, pm, {kmax}, man1);
        Dubins dub2_ = Dubins(pm, pf, {kmax}, man2);
        if (std::abs(dub1_.l() + dub2_.l() - len) < 1e-4){
          std::cout << 
            dub1_.man_to_string() << " " << dub2_.man_to_string() << " " << 
            std::setprecision(12) <<
            (std::abs(dub1_.l() + dub2_.l()) - len) << " " <<
            (dub1_.l() + dub2_.l()) << " " << 
            std::setprecision(4) <<
            dub1_.s1() << " " << dub1_.s2() << " " << dub1_.s2() << " " <<
            dub2_.s1() << " " << dub2_.s2() << " " << dub2_.s3() <<
            std::endl;
        }

      } catch (std::runtime_error &e){
        // std::cout << "Error in " << D_TYPE_STR__[man1] << " " << D_TYPE_STR__[man2] << std::endl;
      }
    }
  }
}

void retesting_man_19(std::string filename){
  std::ifstream file(filename);

  double thi, thf, alpham, alphaf, kmax, len, thm;
  int id_man;

  // Skip first line
  std::string line;
  std::getline(file, line);

  int total = 0, wrong = 0;

  while(file >> kmax >> thi >> thf >> alpham >> alphaf >> thm >> id_man >> len){
    if (id_man == 19){
      Configuration2 pi(1, 0, thi);
      Configuration2 pm(cos(alpham), sin(alpham), 0);
      Configuration2 pf(cos(alphaf), sin(alphaf), thf);

      std::vector<double> curveParam = {kmax};
      std::vector<bool> fixedAngles = {true, false, true};
      std::vector<Configuration2> points = {pi, pm, pf};

      std::pair<LEN_T, std::vector<Angle> >ret=DP().solveDP(points, fixedAngles, curveParam, 2880, 5);
      if (std::abs(ret.first - len) > 1e-4){
        wrong++;
        std::cout << "Error in " << kmax << " " << thi << " " << thf << " " << alpham << " " << alphaf << " " << thm << " " << id_man << " " << len << " " << ret.first << std::endl;
      }
      else {
        std::cout << "Could not fix " << kmax << " " << thi << " " << thf << " " << alpham << " " << alphaf << " " << thm << " " << id_man << " " << len << " " << ret.first << std::endl;
      }
      total++;
    }
  }

  std::cout << wrong << "/" << total << std::endl;
}

std::vector<double>
find_best_circle(
	Configuration2& pi,
	Configuration2& pm,
	Configuration2& pf,
	std::vector<bool> fixedAngles,
	std::vector<double> curveParam,
	int discr,
	int refinements
){
	std::vector<Configuration2> points = {pi, pm, pf};
	double kmax = curveParam[0];
	std::pair<LEN_T, std::vector<Angle> > ret = DP().solveDP(points, fixedAngles, curveParam, discr, refinements);
	if (ret.first == 0.0) {
		std::cout << pi << std::endl << pm << std::endl << pf << std::endl;
		throw std::runtime_error("Zero length");
	}
	// Set angle for intermediate problem and compute the two Dubins
	pm.th(ret.second[1]);
	Dubins dub1 = Dubins(pi, pm, kmax);
	Dubins dub2 = Dubins(pm, pf, kmax);
	// std::cout << "Took " << dtime << " ms to find Dubins" << std::endl;
	LEN_T len = dub1.l() + dub2.l();
	// Get the manoeuver combination, and if it's not in the 18 valid ones, search for an alternative
	std::string man_comb = dub1.man_to_string() + dub2.man_to_string();
	int id_man_comb = 19;
	auto search = P3DP_DICT.find(man_comb);
	if (search == P3DP_DICT.end()) {
    std::pair<std::string, double> shortest;
    shortest.first = "";
    shortest.second = std::numeric_limits<double>::infinity();
	  std::string D_TYPE_STR_[7] = {"INVALID", "LRL", "RLR", "LSL", "LSR", "RSL", "RSR"};
		for (auto man: P3DP_DICT) {
			Dubins::D_TYPE dub1_man = std::get<1>(man.second);
			Dubins::D_TYPE dub2_man = std::get<2>(man.second);
			try {
				Dubins dub1 = Dubins(pi, pm, { kmax }, dub1_man);
				Dubins dub2 = Dubins(pm, pf, { kmax }, dub2_man);
				if (std::abs(dub1.l() + dub2.l() - len) < 1e-8) {
					id_man_comb = std::get<0>(man.second);
					break;
				}
        else if (shortest.second > std::abs(dub1.l() + dub2.l() - len)) {
          shortest.first = dub1.man_to_string() + dub2.man_to_string();
          shortest.second = std::abs(dub1.l() + dub2.l() - len);
        }
			}
			catch (std::runtime_error &e) {
				continue;
			}
		}
    if (shortest.first != "") {
      std::cout << "Maneuver combination not found in the dictionary. Closest is " << shortest.first << " with error " << shortest.second << " out of " << len << std::endl;
    }
	} else {
		id_man_comb = std::get<0>(search->second);
	}
	return {static_cast<double>(id_man_comb), len};
}

/**
 * @brief Generates a dataset of 3PDP problems with the circle constraint.
 *
 * @param argc The number of arguments, either 1, 4 or 5. Since they are passed directly from the command line, argc is always at least 1.
 *             If argc is 4, the arguments are kmax_min, kmax_max, k_discr. If argc is 5, the arguments are kmax_min, kmax_max, k_discr, angle_discr.
 * @param argv
 */
void generateDataset3PDPCircle(int argc, char** argv){
  double kmax_min = 1;
  double kmax_max = 1;
  double k_step = 1;
  int angle_discr = 5;

  if (argc == 4) {
    kmax_min = std::stof(argv[1]);
    kmax_max = std::stof(argv[2]);
    k_step  = std::stof(argv[3]);
  }
  else if (argc == 5) {
    kmax_min = std::stof(argv[1]);
    kmax_max = std::stof(argv[2]);
    k_step  = std::stof(argv[3]);
    angle_discr = std::stoi(argv[4]);
  }

  // Open file named 3PDS.csv
  std::string filename_base = "3PDS_Circle" + std::to_string(angle_discr) + "_" + std::to_string(kmax_min) + "_"  + std::to_string(kmax_max) + "_"  + std::to_string(k_step);
  std::string filename = filename_base + ".csv";
  std::string filename_log = filename_base + ".log";

  uint64_t counter = 0;
  uint64_t prev_counter = 0;
  uint64_t actual_counter = 0;

	uint64_t k_discr = (kmax_max-kmax_min)/k_step;
	std::vector<double> k_discrs (k_discr, 0.0);
	double dth = 2.0 * m_pi / angle_discr;
	std::vector<double> th_discrs (angle_discr, dth/2.0);

	std::generate(k_discrs.begin(), k_discrs.end(), [k_step, kmax_tmp = kmax_min]() mutable {
		return (kmax_tmp += k_step);
	});

	std::generate(th_discrs.begin(), th_discrs.end(), [dth, th = m_pi]() mutable{
		return (th -= dth);
	});

  uint64_t tot_counter = th_discrs.size()*th_discrs.size()*th_discrs.size()*th_discrs.size()*k_discrs.size();
  std::cout << "Generating " << PrintScientificLargeInt(tot_counter) << " tests." << std::endl;

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

  file << "kmax" << " " << "theta_i" << " " << "theta_f" << " " << "alpha_m" << " "
       << "alpha_f" << " " << "th_m" << " " << "id_man_comb" << " " << "len" << std::endl;

  std::cout << "k_discrs: " << k_discrs.size() << std::endl;
  for(auto kmax_tmp : k_discrs){
    std::cout << kmax_tmp << " ";
  }
  std::cout << std::endl;
  std::cout << "th_discrs: " << th_discrs.size() << std::endl;
  for(auto th : th_discrs){
    std::cout << th << " ";
  }
  std::cout << "Total: " << k_discrs.size()*th_discrs.size()*th_discrs.size()*th_discrs.size()*th_discrs.size() << std::endl;

  for (auto kmax : k_discrs){
    TimePerf time1; time1.start();
    for (double theta_i : th_discrs){
      for (double theta_f : th_discrs) {
        for (double alpha_m : th_discrs) {
          for (double alpha_f : th_discrs) {
            Configuration2 pi = Configuration2(1, 0, theta_i);
            Configuration2 pm = Configuration2(cos(alpha_m), sin(alpha_m), 0);
            Configuration2 pf = Configuration2(cos(alpha_f), sin(alpha_f), theta_f);

            if (pm.x() != pi.x() && pm.y() != pi.y() && pm.x() != pf.x() && pm.y() != pf.y()){
              counter ++;
              // Solve multipoint problem
              std::vector<Configuration2> points = {pi, pm, pf};
              std::vector<bool> fixedAngles = {true, false, true};
              std::vector<double> curveParam = { kmax };
              int discr = 90;
              int refinements = 4;
              TimePerf time;
              time.start();

            	std::vector<double> res = find_best_circle(pi, pm, pf, fixedAngles, curveParam, discr, refinements);

              auto dtime = time.getTime();

            	int id_man_comb = static_cast<int>(res[0]);
            	double len = res[1];

              // Write data to file
              file << std::setprecision(5) << kmax << " " << theta_i << " " << theta_f << " " << alpha_m << " "
                    << alpha_f << " " << pm.th() << " " << id_man_comb << " " << len << std::endl;
            }

            // Print time
            auto dtime1 = time1.getTime();
          	auto part = tot_counter > 100 ? tot_counter/100 : 1;
            if (counter % part == 0) {
              std::cout << 100.0 * counter / tot_counter << "% " << counter << " in " << dtime1 << "ms, avg " << (dtime1/(1.0*(counter-prev_counter))) << "ms" << std::endl;
              prev_counter = counter;
              time1.start();
            }
            counter ++;
          }
        }
      }
    }
  }

  std::cout << "Generated " << PrintScientificLargeInt(actual_counter) << " entries to " << filename << std::endl;

  if (coutbuf != nullptr){
    std::cout.rdbuf(coutbuf);
  }

  file.close();
}


/**
 * @brief Generates a dataset of 3PDP problems with the circle constraint.
 *
 * @param argc The number of arguments, either 1, 4 or 5. Since they are passed directly from the command line, argc is always at least 1.
 *             If argc is 4, the arguments are kmax_min, kmax_max, k_discr. If argc is 5, the arguments are kmax_min, kmax_max, k_discr, angle_discr.
 * @param argv
 */
void generateDataset3PDPCircleRandom(int argc, char** argv){
  if (argc != 5) {
    std::cout << "Usage: " << argv[0] << " kmax_min kmax_max k_step angle_discr" << std::endl;
    return;
  }
  double kmax_min = std::atof(argv[1]);
  double kmax_max = std::atof(argv[2]);
  int k_discr_in  = std::atoi(argv[3]);
  int angle_discr = std::atoi(argv[4]);

  std::uniform_real_distribution<double> th_distribution(-m_pi, m_pi);
  std::uniform_real_distribution<double> k_distribution(kmax_min, kmax_max);

  std::mt19937 rng;
  rng.seed(41);

  // new_argv[1] = strdup ("-0.1");
  // new_argv[2] = strdup ("8.3");
  // new_argv[3] = strdup ("0.3");
  // new_argv[4] = strdup ("18.0");

  std::string filename_base = "3PDS_Circle_random" + std::to_string(angle_discr) + "_" + std::to_string(kmax_min) + "_"  + std::to_string(kmax_max) + "_"  + std::to_string(k_discr_in);
  std::string filename = filename_base + ".csv";
  std::string filename_log = filename_base + ".log";

  uint64_t counter = 0;
  uint64_t actual_counter = 0;
  uint64_t comp_counter = 0;

  uint64_t k_discr = static_cast<uint64_t>(k_discr_in);
  std::vector<double> k_discrs;
  std::vector<double> thi_discrs;
  std::vector<double> thf_discrs;
  std::vector<double> alpham_discrs;
  std::vector<double> alphaf_discrs;

  for (uint64_t i = 0; i < k_discr; i++){
    k_discrs.push_back(k_distribution(rng));
  }
  for (uint64_t i = 0; i < angle_discr; i++){
    thi_discrs.push_back(th_distribution(rng));
    thf_discrs.push_back(th_distribution(rng));
    alpham_discrs.push_back(th_distribution(rng));
    alphaf_discrs.push_back(th_distribution(rng));
  }

  uint64_t tot_counter = thi_discrs.size()*thf_discrs.size()*alpham_discrs.size()*alphaf_discrs.size()*k_discrs.size();
  std::cout << "Generating " << PrintScientificLargeInt(tot_counter) << " tests." << std::endl;

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

  file << "kmax" << " " << "theta_i" << " " << "theta_f" << " " << "alpha_m" << " "
      << "alpha_f" << " " << "th_m" << " " << "id_man_comb" << " " << "len" << std::endl;

  std::cout << "k_discrs: " << k_discrs.size() << std::endl;
  for(auto kmax_tmp : k_discrs){
    std::cout << kmax_tmp << " ";
  }
  std::cout << std::endl;
  std::cout << "th_discrs: " << thi_discrs.size() << std::endl;
  std::cout << "thi: " << std::endl;
  for(auto th : thi_discrs){
    std::cout << th << " ";
  }
  std::cout << "thf: " << std::endl;
  for(auto th : thf_discrs){
    std::cout << th << " ";
  }
  std::cout << "alpham: " << std::endl;
  for(auto th : alpham_discrs){
    std::cout << th << " ";
  }
  std::cout << "alphaf: " << std::endl;
  for(auto th : alphaf_discrs){
    std::cout << th << " ";
  }
  std::cout << "Total: " << k_discrs.size()*thi_discrs.size()*thf_discrs.size()*alpham_discrs.size()*alphaf_discrs.size() << std::endl;

  unsigned long time_sd = 0;
  unsigned long time_ps = 0;

  int part = tot_counter/100;

  auto start = std::chrono::high_resolution_clock::now();
  for (auto kmax : k_discrs){
    for (double theta_i : thi_discrs){
      for (double theta_f : thf_discrs) {
        for (double alpha_m : alpham_discrs) {
          for (double alpha_f : alphaf_discrs) {
            Configuration2 pi = Configuration2(1, 0, theta_i);
            Configuration2 pm = Configuration2(cos(alpha_m), sin(alpha_m), 0);
            Configuration2 pf = Configuration2(cos(alpha_f), sin(alpha_f), theta_f);

            if (pm.x() != pi.x() && pm.y() != pi.y() && pm.x() != pf.x() && pm.y() != pf.y()){
              counter ++;
              // Solve multipoint problem
              std::vector<Configuration2> points = {pi, pm, pf};
              std::vector<bool> fixedAngles = {true, false, true};
              std::vector<double> curveParam = { kmax };
              int discr = 90;
              int refinements = 4;

            	std::vector<double> res = find_best_circle(pi, pm, pf, fixedAngles, curveParam, discr, refinements);

            	int id_man_comb = static_cast<int>(res[0]);
            	double len = res[1];

              if (id_man_comb < 19 && id_man_comb > 0){
                actual_counter ++;
                // Write data to file
                file << std::setprecision(5) << kmax << " " << theta_i << " " << theta_f << " " << alpha_m << " "
                      << alpha_f << " " << pm.th() << " " << id_man_comb << " " << len << std::endl;
              }
            }

            auto dtime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();

            // Print time
          	auto part = tot_counter > 100 ? tot_counter/100 : 1;
            if (counter % part == 0) {
              std::cout << 100.0 * counter / tot_counter << "% " << counter << " in " << dtime << "ms, avg " << (dtime/(1.0*counter)) << "ms" << std::endl;
            }
            counter ++;
          }
        }
      }
    }
  }

  std::cout << "Generated " << PrintScientificLargeInt(actual_counter) << " entries to " << filename << std::endl;

  if (coutbuf != nullptr){
    std::cout.rdbuf(coutbuf);
  }

  file.close();
}


/**
 * @brief Generates a dataset of 3PDP problems with the circle constraint.
 *
 * @param argc The number of arguments, either 1, 4 or 5. Since they are passed directly from the command line, argc is always at least 1.
 *             If argc is 4, the arguments are kmax_min, kmax_max, k_discr. If argc is 5, the arguments are kmax_min, kmax_max, k_discr, angle_discr.
 * @param argv
 */
void generateDataset3PDPCircleWithAllLabels(int argc, char** argv){
  // double kmax_min = 1;
  // double kmax_max = 1;
  // double k_step = 1;
  // int angle_discr = 5;

  // if (argc == 4) {
  //   kmax_min = std::stof(argv[1]);
  //   kmax_max = std::stof(argv[2]);
  //   k_step  = std::stof(argv[3]);
  // }
  // else if (argc == 5) {
  //   kmax_min = std::stof(argv[1]);
  //   kmax_max = std::stof(argv[2]);
  //   k_step  = std::stof(argv[3]);
  //   angle_discr = std::stoi(argv[4]);
  // }

  // // Open file named 3PDS.csv
  // std::string filename_base = "3PDS_CircleLabels" + std::to_string(angle_discr) + "_" + std::to_string(kmax_min) + "_"  + std::to_string(kmax_max) + "_"  + std::to_string(k_step);
  // std::string filename = filename_base + ".csv";
  // std::string filename_log = filename_base + ".log";

  // uint64_t counter = 0;
  // uint64_t prev_counter = 0;
  // uint64_t actual_counter = 0;

	// uint64_t k_discr = (kmax_max-kmax_min)/k_step;
	// std::vector<double> k_discrs (k_discr, 0.0);
	// double dth = 2.0 * m_pi / angle_discr;
	// std::vector<double> th_discrs (angle_discr, dth/2.0);

	// std::generate(k_discrs.begin(), k_discrs.end(), [k_step, kmax_tmp = kmax_min]() mutable {
	// 	return (kmax_tmp += k_step);
	// });

	// std::generate(th_discrs.begin(), th_discrs.end(), [dth, th = m_pi]() mutable{
	// 	return (th -= dth);
	// });

  // uint64_t tot_counter = th_discrs.size()*th_discrs.size()*th_discrs.size()*th_discrs.size()*k_discrs.size();
  // std::cout << "Generating " << PrintScientificLargeInt(tot_counter) << " tests." << std::endl;

  // std::cout << "Writing entries to " << filename << std::endl;
  // std::ofstream file(filename);
  // if (!file.is_open()) {
  //   std::cout << "Error opening db " << filename << std::endl;
  //   return;
  // }

  // std::cout << "Writing log to " << filename_log << std::endl;
  // std::ofstream log_file(filename_log);
  // std::streambuf* coutbuf = nullptr;
  // if (!log_file.is_open()) {
  //   std::cout << "Error opening log file " << filename_log << std::endl;
  //   return;
  // }
  // else {
  //   coutbuf = std::cout.rdbuf();
  //   std::cout.rdbuf(log_file.rdbuf());
  // }

  // std::cout << "Generating " << PrintScientificLargeInt(tot_counter) << " tests" << std::endl;

  // file << "kmax" << " " << "theta_i" << " " << "theta_f" << " " << "alpha_m" << " "
  //      << "alpha_f" << " " << "th_m" << " " << "id_man_comb" << " " << "len" << std::endl;

  // std::cout << "k_discrs: " << k_discrs.size() << std::endl;
  // for(auto kmax_tmp : k_discrs){
  //   std::cout << kmax_tmp << " ";
  // }
  // std::cout << std::endl;
  // std::cout << "th_discrs: " << th_discrs.size() << std::endl;
  // for(auto th : th_discrs){
  //   std::cout << th << " ";
  // }
  // std::cout << "Total: " << k_discrs.size()*th_discrs.size()*th_discrs.size()*th_discrs.size()*th_discrs.size() << std::endl;

  // for (auto kmax : k_discrs){
  //   TimePerf time1; time1.start();
  //   for (double theta_i : th_discrs){
  //     for (double theta_f : th_discrs) {
  //       for (double alpha_m : th_discrs) {
  //         for (double alpha_f : th_discrs) {
  //           Configuration2 pi = Configuration2(1, 0, theta_i);
  //           Configuration2 pm = Configuration2(cos(alpha_m), sin(alpha_m), 0);
  //           Configuration2 pf = Configuration2(cos(alpha_f), sin(alpha_f), theta_f);

  //           if (pm.x() != pi.x() && pm.y() != pi.y() && pm.x() != pf.x() && pm.y() != pf.y()){
  //             counter ++;
  //             // Solve multipoint problem
  //             std::vector<Configuration2> points = {pi, pm, pf};
  //             std::vector<bool> fixedAngles = {true, false, true};
  //             std::vector<double> curveParam = { kmax };
  //             int discr = 90;
  //             int refinements = 4;
  //             TimePerf time;
  //             time.start();

  //             std::vector<std::pair<double, int>> results;
  //             for(int man1 = 1; man1 < 7; man1++){
  //               for(int man2 = 1; man2 < 7; man2++){
  //                 Dubins::D_TYPE man1_type = static_cast<Dubins::D_TYPE>(man1);
  //                 Dubins::D_TYPE man2_type = static_cast<Dubins::D_TYPE>(man2);
  //                 Dubins dub1 = Dubins(pi, pm, curveParam, man1_type);
  //                 Dubins dub2 = Dubins(pm, pf, curveParam, man2_type);
  //                 LEN_T len = dub1.l() + dub2.l();
  //                 file << std::setprecision(5) << kmax << " " << theta_i << " " << theta_f << " " << alpha_m << " "
  //                      << alpha_f << " " << pm.th() << " " << man1 << " " << man2 << " " << len << std::endl;
  //               }
  //             }


  //             auto dtime = time.getTime();

  //           	int id_man_comb = static_cast<int>(res[0]);
  //           	double len = res[1];

  //             // Write data to file
  //             file << std::setprecision(5) << kmax << " " << theta_i << " " << theta_f << " " << alpha_m << " "
  //                   << alpha_f << " " << pm.th() << " " << id_man_comb << " " << len << std::endl;
  //           }

  //           // Print time
  //           auto dtime1 = time1.getTime();
  //         	auto part = tot_counter > 100 ? tot_counter/100 : 1;
  //           if (counter % part == 0) {
  //             std::cout << 100.0 * counter / tot_counter << "% " << counter << " in " << dtime1 << "ms, avg " << (dtime1/(1.0*(counter-prev_counter))) << "ms" << std::endl;
  //             prev_counter = counter;
  //             time1.start();
  //           }
  //           counter ++;
  //         }
  //       }
  //     }
  //   }
  // }

  // std::cout << "Generated " << PrintScientificLargeInt(actual_counter) << " entries to " << filename << std::endl;

  // if (coutbuf != nullptr){
  //   std::cout.rdbuf(coutbuf);
  // }

  // file.close();
}

/**
 * 1) con curvatura che prende valori a metà di quelli che abbiamo usato, cioè se abbiamo
 * usato 0.2 0.4 0.6 ecc fammi un set con 0.3 0.5 0.7 ecc. e angoli possibilmente anche
 * sfasati, tipo prendi gli angoli del train e sfasali di mezzo delta_angolo 2) uno meno
 * cattivo, con dati più vicini al train set
 * @param argc
 * @param argv
 */
void
generateDataset3PDPCircleTest (int argc, char** argv){
	if (argc != 2){
		throw std::runtime_error("Invalid number of arguments, expected 1");
	}
	char** new_argv = new char*[5];
	for (int i = 0; i < 5; i++){
		new_argv[i] = new char[10];
	}
	new_argv[0] = "";
	if (std::stoi(argv[1]) == 1 || std::stoi(argv[1]) == 3)
	{
		new_argv[1] = strdup ("0.1");
		new_argv[2] = strdup ("8.1");
		new_argv[3] = strdup ("0.2");
		new_argv[4] = strdup ("36.0");
		generateDataset3PDPCircle (5, new_argv);
	}
	else if (std::stoi(argv[1]) == 2 || std::stoi(argv[1]) == 3){
		new_argv[1] = strdup ("-0.1");
		new_argv[2] = strdup ("8.3");
		new_argv[3] = strdup ("0.3");
		new_argv[4] = strdup ("18.0");
		generateDataset3PDPCircle (5, new_argv);
	}
	else {
		throw std::runtime_error("Invalid argument");
	}
	delete[] new_argv;
}



/** @brief Generates a dataset of 3PDP problems with the rectangle constraint.
*
* @param argc The number of arguments, either 1, 4 or 5. Since they are passed directly from the command line, argc is always at least 1.
*             If argc is 4, the arguments are kmax_min, kmax_max, k_discr. If argc is 5, the arguments are kmax_min, kmax_max, k_discr, angle_discr.
* @param argv
*/
void generateDataset3PDPRect(int argc, char** argv){
  if (argc != 9) {
    std::cout << "Usage: " << argv[0] << " xi_discr xm_discr ym_discr xf_discr kmax_min kmax_max k_discr angle_discr" << std::endl;
    return ;
  }
  int xi_discr    = std::atoi(argv[1]);
  int xm_discr    = std::atoi(argv[2]);
  int ym_discr    = std::atoi(argv[3]);
  int xf_discr    = std::atoi(argv[4]);
  double kmax_min = std::atof(argv[5]);
  double kmax_max = std::atof(argv[6]);
  int k_discr_in  = std::atoi(argv[7]);
  int angle_discr = std::atoi(argv[8]);

  std::uniform_real_distribution<double> th_distribution(-m_pi, m_pi);
  std::uniform_real_distribution<double> k_distribution(kmax_min, kmax_max);
  std::uniform_real_distribution<double> xc_distribution(0.5, 1.0+0.5);
  std::uniform_real_distribution<double> m_distribution(-0.5, 1.0+0.5);

  std::mt19937 rng;
  rng.seed(41);
 
  std::string filename_base = "3PDS_Rect_random" + 
            std::to_string(angle_discr) + "_" + std::to_string(xi_discr) + "_" + 
            std::to_string(xm_discr) + "_" + std::to_string(ym_discr) + "_" +
            std::to_string(xf_discr) + "_" + std::to_string(kmax_min) + "_"  + 
            std::to_string(kmax_max) + "_"  + std::to_string(k_discr_in);
  std::string filename = filename_base + ".csv";
  std::string filename_log = filename_base + ".log";

  uint64_t tot_counter = xi_discr*xm_discr*ym_discr*xf_discr*angle_discr*angle_discr*k_discr_in;
  std::cout << "Generating at most " << PrintScientificLargeInt(tot_counter) << " tests." << std::endl;

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

  file << "kmax" << " " << "xi" << " " << "xm" << " " << "ym" << " " << "xf" << " " << 
          "theta_i" << " " << "theta_f" << " " << "th_m" << " " << "id_man_comb" << " " << "len" << std::endl;
 
  uint64_t counter = 0;
  uint64_t actual_counter = 0;
  uint64_t comp_counter = 0;

  unsigned long time_sd = 0;
  unsigned long time_ps = 0;

  int part = tot_counter/100;

  std::vector<double> k_discrs (k_discr_in+1, 0.0);
  std::vector<double> thi_discrs (angle_discr, 0.0);
  std::vector<double> thf_discrs (angle_discr, 0.0);
  std::vector<double> xi_discrs (xi_discr, 0.0);
  std::vector<double> xm_discrs (xm_discr+1, 0.0);
  std::vector<double> ym_discrs (ym_discr+1, 0.0);
  std::vector<double> xf_discrs (xf_discr, 0.0);
  double k_step = (kmax_max - kmax_min) / (double)(k_discr_in);
  std::generate(k_discrs.begin(), k_discrs.end(), [kmax_min, kmax_tmp = kmax_min, kmax_max = kmax_max, k_step=k_step]() mutable {
    double k = kmax_tmp;
    kmax_tmp += k_step;
    return k;
  });
  std::generate(thi_discrs.begin(), thi_discrs.end(), [angle_discr, th = -m_pi]() mutable{
    double th_i = th;
    th += (2.0 * m_pi) / (double)angle_discr;
    return th_i;
  });
  std::generate(thf_discrs.begin(), thf_discrs.end(), [angle_discr, th = -m_pi]() mutable{
    double th_f = th;
    th += (2.0 * m_pi) / (double)angle_discr;
    return th_f;
  });
  double xi_step = 1.0 / (double)xi_discr;
  std::generate(xi_discrs.begin(), xi_discrs.end(), [xi_step, xi_tmp = xi_step]() mutable {
    double xi = xi_tmp;
    xi_tmp += xi_step;
    return xi;
  });
  double xm_step = 1.0 / (double)xm_discr;
  double first_step = 0;
  std::generate(xm_discrs.begin(), xm_discrs.end(), [xm_step, xm_tmp = first_step]() mutable {
    double xm = xm_tmp;
    xm_tmp += xm_step;
    return xm;
  });
  double ym_step = 1.0 / (double)ym_discr;
  std::generate(ym_discrs.begin(), ym_discrs.end(), [ym_step, ym_tmp = first_step]() mutable {
    double ym = ym_tmp;
    ym_tmp += ym_step;
    return ym;
  });
  double xf_step = 1.0 / (double)xf_discr;
  std::generate(xf_discrs.begin(), xf_discrs.end(), [xf_step, xf_tmp = xf_step]() mutable {
    double xf = xf_tmp;
    xf_tmp += xf_step;
    return xf;
  });
  std::cout << "k_discrs: " << k_discrs.size() << std::endl;
  for(auto kmax_tmp : k_discrs){
    std::cout << kmax_tmp << " ";
  }
  std::cout << std::endl;
  std::cout << "th_discrs: " << thi_discrs.size() << std::endl;
  for(auto th : thi_discrs){
    std::cout << th << " ";
  }
  std::cout << std::endl;
  std::cout << "th_discrs: " << thf_discrs.size() << std::endl;
  for(auto th : thf_discrs){
    std::cout << th << " ";
  }
  std::cout << std::endl;
  std::cout << "xi_discrs: " << xi_discrs.size() << std::endl;
  for(auto xi : xi_discrs){
    std::cout << xi << " ";
  }
  std::cout << std::endl;
  std::cout << "xm_discrs: " << xm_discrs.size() << std::endl;
  for(auto xm : xm_discrs){
    std::cout << xm << " ";
  }
  std::cout << std::endl;
  std::cout << "ym_discrs: " << ym_discrs.size() << std::endl;
  for(auto ym : ym_discrs){
    std::cout << ym << " ";
  }
  std::cout << std::endl;
  std::cout << "xf_discrs: " << xf_discrs.size() << std::endl;
  for(auto xf : xf_discrs){
    std::cout << xf << " ";
  }
  std::cout << std::endl;
  std::cout << "Total: " << k_discrs.size()*xi_discrs.size()*xm_discrs.size()*ym_discrs.size()*xf_discrs.size()*thi_discrs.size()*thf_discrs.size() << std::endl;

  std::cout << "Part: " << part << std::endl;

  auto start = std::chrono::high_resolution_clock::now();
  for(auto k_max : k_discrs){
    for(auto xi : xi_discrs){
      for(auto xm : xm_discrs){
        for(auto ym : ym_discrs){
          for(auto xf : xf_discrs){
            for(auto theta_i : thi_discrs){
              for(auto theta_f : thf_discrs){
                Configuration2 pi = Configuration2(xi, 0, theta_i);
                Configuration2 pm = Configuration2(xm, ym, 0);
                Configuration2 pf = Configuration2(xf, 0, theta_f);

                if (pm.x() != pi.x() && pm.y() != pi.y() && pm.x() != pf.x() && pm.y() != pf.y()){
                  counter ++;
                  // Solve multipoint problem
                  std::vector<Configuration2> points = {pi, pm, pf};
                  std::vector<bool> fixedAngles = {true, false, true};
                  std::vector<double> curveParam = { k_max };
                  int discr = 90;
                  int refinements = 4;

                  std::vector<double> res = find_best_circle(pi, pm, pf, fixedAngles, curveParam, discr, refinements);

                  int id_man_comb = static_cast<int>(res[0]);
                  double len = res[1];

                  if (id_man_comb < 19 && id_man_comb > 0){
                    actual_counter ++;
                    // Write data to file
                    file << std::setprecision(5) << k_max << " " << xi << " " << xm << " " << ym << " " << xf << " " << 
                            theta_i << " " << theta_f << " " << pm.th() << " " << id_man_comb << " " << len << std::endl;
                  }
                }

                auto dtime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();

                // Print time
                auto part = tot_counter > 100 ? tot_counter/100 : 1;
                if (counter % part == 0) {
                  std::cout << 100.0 * counter / tot_counter << "% " << counter << " in " << dtime << "ms, avg " << (dtime/(1.0*counter)) << "ms" << std::endl;
                }
                counter ++;
              }
            }
          }
        }
      }
    }
  }

  std::cout << "Generated " << PrintScientificLargeInt(actual_counter) << " entries to " << filename << std::endl;

  if (coutbuf != nullptr){
    std::cout.rdbuf(coutbuf);
  }

  file.close();
 }




std::vector<std::string>
find_all_best(
	Configuration2& pi,
	Configuration2& pm,
	Configuration2& pf,
	std::vector<bool> fixedAngles,
	std::vector<double> curveParam,
	int discr,
	int refinements
){
	std::vector<Configuration2> points = {pi, pm, pf};
	double kmax = curveParam[0];
	std::pair<LEN_T, std::vector<Angle> > ret = DP().solveDP(points, fixedAngles, curveParam, discr, refinements);
	if (ret.first == 0.0) {
		std::cout << pi << std::endl << pm << std::endl << pf << std::endl;
		throw std::runtime_error("Zero length");
	}
  double best_len = ret.first;
	// Set angle for intermediate problem and compute all possible combinations of Dubins
	pm.th(ret.second[1]);

  std::vector<Dubins::D_TYPE> man_types = {
    Dubins::D_TYPE::LRL,
    Dubins::D_TYPE::RLR,
    Dubins::D_TYPE::LSL,
    Dubins::D_TYPE::LSR,
    Dubins::D_TYPE::RSL,
    Dubins::D_TYPE::RSR};
  
  std::vector<int> best_ids = {};
  std::vector<std::string> best_man = {};

  for (auto man_type_d1 : man_types){
    for (auto man_type_d2 : man_types){
      try {
        Dubins dub1 = Dubins(pi, pm, curveParam, man_type_d1);
        Dubins dub2 = Dubins(pm, pf, curveParam, man_type_d2);
        LEN_T len = dub1.l() + dub2.l();
        if (std::abs(len-best_len) < 1e-8) {
          std::string man_comb = dub1.man_to_string() + dub2.man_to_string();
          best_man.push_back(man_comb);
          // int id_man_comb = 19;
          // auto search = P3DP_DICT.find(man_comb);
          // if (search != P3DP_DICT.end()) {
          //   best_ids.push_back(std::get<0>(search->second));
          // }
        }
      } catch (...){}
    }
  }

  // return best_ids;
  return best_man;
}


 /** @brief Generates a dataset of 3PDP problems with the rectangle constraint setting for each entry an array of classes
*
* @param argc The number of arguments, either 1, 4 or 5. Since they are passed directly from the command line, argc is always at least 1.
*             If argc is 4, the arguments are kmax_min, kmax_max, k_discr. If argc is 5, the arguments are kmax_min, kmax_max, k_discr, angle_discr.
* @param argv
*/
void generateDataset3PDPRectMulti(int argc, char** argv){
  if (argc != 9) {
    std::cout << "Usage: " << argv[0] << " xi_discr xm_discr ym_discr xf_discr kmax_min kmax_max k_discr angle_discr" << std::endl;
    return ;
  }
  int xi_discr    = std::atoi(argv[1]);
  int xm_discr    = std::atoi(argv[2]);
  int ym_discr    = std::atoi(argv[3]);
  int xf_discr    = std::atoi(argv[4]);
  double kmax_min = std::atof(argv[5]);
  double kmax_max = std::atof(argv[6]);
  int k_discr_in  = std::atoi(argv[7]);
  int angle_discr = std::atoi(argv[8]);

  std::uniform_real_distribution<double> th_distribution(-m_pi, m_pi);
  std::uniform_real_distribution<double> k_distribution(kmax_min, kmax_max);
  std::uniform_real_distribution<double> xc_distribution(0.5, 1.0+0.5);
  std::uniform_real_distribution<double> m_distribution(-0.5, 1.0+0.5);

  std::mt19937 rng;
  rng.seed(41);
 
  std::string filename_base = "3PDS_Rect_multi" + 
            std::to_string(angle_discr) + "_" + std::to_string(xi_discr) + "_" + 
            std::to_string(xm_discr) + "_" + std::to_string(ym_discr) + "_" +
            std::to_string(xf_discr) + "_" + std::to_string(kmax_min) + "_"  + 
            std::to_string(kmax_max) + "_"  + std::to_string(k_discr_in);
  std::string filename = filename_base + ".csv";
  std::string filename_log = filename_base + ".log";

  uint64_t tot_counter = xi_discr*xm_discr*ym_discr*xf_discr*angle_discr*angle_discr*k_discr_in;
  std::cout << "Generating at most " << PrintScientificLargeInt(tot_counter) << " tests." << std::endl;

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

  file << "kmax" << " " << "xi" << " " << "xm" << " " << "ym" << " " << "xf" << " " << 
          "theta_i" << " " << "theta_f" << " " << "th_m" << " " << "id_man_comb" << std::endl;
 
  uint64_t counter = 0;
  uint64_t actual_counter = 0;
  uint64_t comp_counter = 0;

  unsigned long time_sd = 0;
  unsigned long time_ps = 0;

  int part = tot_counter/100;

  std::vector<double> k_discrs (k_discr_in+1, 0.0);
  std::vector<double> thi_discrs (angle_discr, 0.0);
  std::vector<double> thf_discrs (angle_discr, 0.0);
  std::vector<double> xi_discrs (xi_discr, 0.0);
  std::vector<double> xm_discrs (xm_discr+1, 0.0);
  std::vector<double> ym_discrs (ym_discr+1, 0.0);
  std::vector<double> xf_discrs (xf_discr, 0.0);
  double k_step = (kmax_max - kmax_min) / (double)(k_discr_in);
  std::generate(k_discrs.begin(), k_discrs.end(), [kmax_min, kmax_tmp = kmax_min, kmax_max = kmax_max, k_step=k_step]() mutable {
    double k = kmax_tmp;
    kmax_tmp += k_step;
    return k;
  });
  std::generate(thi_discrs.begin(), thi_discrs.end(), [angle_discr, th = -m_pi]() mutable{
    double th_i = th;
    th += (2.0 * m_pi) / (double)angle_discr;
    return th_i;
  });
  std::generate(thf_discrs.begin(), thf_discrs.end(), [angle_discr, th = -m_pi]() mutable{
    double th_f = th;
    th += (2.0 * m_pi) / (double)angle_discr;
    return th_f;
  });
  double xi_step = 1.0 / (double)xi_discr;
  std::generate(xi_discrs.begin(), xi_discrs.end(), [xi_step, xi_tmp = xi_step]() mutable {
    double xi = xi_tmp;
    xi_tmp += xi_step;
    return xi;
  });
  double xm_step = 1.0 / (double)xm_discr;
  double first_step = 0;
  std::generate(xm_discrs.begin(), xm_discrs.end(), [xm_step, xm_tmp = first_step]() mutable {
    double xm = xm_tmp;
    xm_tmp += xm_step;
    return xm;
  });
  double ym_step = 1.0 / (double)ym_discr;
  std::generate(ym_discrs.begin(), ym_discrs.end(), [ym_step, ym_tmp = first_step]() mutable {
    double ym = ym_tmp;
    ym_tmp += ym_step;
    return ym;
  });
  double xf_step = 1.0 / (double)xf_discr;
  std::generate(xf_discrs.begin(), xf_discrs.end(), [xf_step, xf_tmp = xf_step]() mutable {
    double xf = xf_tmp;
    xf_tmp += xf_step;
    return xf;
  });
  std::cout << "k_discrs: " << k_discrs.size() << std::endl;
  for(auto kmax_tmp : k_discrs){
    std::cout << kmax_tmp << " ";
  }
  std::cout << std::endl;
  std::cout << "th_discrs: " << thi_discrs.size() << std::endl;
  for(auto th : thi_discrs){
    std::cout << th << " ";
  }
  std::cout << std::endl;
  std::cout << "th_discrs: " << thf_discrs.size() << std::endl;
  for(auto th : thf_discrs){
    std::cout << th << " ";
  }
  std::cout << std::endl;
  std::cout << "xi_discrs: " << xi_discrs.size() << std::endl;
  for(auto xi : xi_discrs){
    std::cout << xi << " ";
  }
  std::cout << std::endl;
  std::cout << "xm_discrs: " << xm_discrs.size() << std::endl;
  for(auto xm : xm_discrs){
    std::cout << xm << " ";
  }
  std::cout << std::endl;
  std::cout << "ym_discrs: " << ym_discrs.size() << std::endl;
  for(auto ym : ym_discrs){
    std::cout << ym << " ";
  }
  std::cout << std::endl;
  std::cout << "xf_discrs: " << xf_discrs.size() << std::endl;
  for(auto xf : xf_discrs){
    std::cout << xf << " ";
  }
  std::cout << std::endl;
  std::cout << "Total: " << k_discrs.size()*xi_discrs.size()*xm_discrs.size()*ym_discrs.size()*xf_discrs.size()*thi_discrs.size()*thf_discrs.size() << std::endl;

  std::cout << "Part: " << part << std::endl;

  auto start = std::chrono::high_resolution_clock::now();
  for(auto k_max : k_discrs){
    for(auto xi : xi_discrs){
      for(auto xm : xm_discrs){
        for(auto ym : ym_discrs){
          for(auto xf : xf_discrs){
            for(auto theta_i : thi_discrs){
              for(auto theta_f : thf_discrs){
                Configuration2 pi = Configuration2(xi, 0, theta_i);
                Configuration2 pm = Configuration2(xm, ym, 0);
                Configuration2 pf = Configuration2(xf, 0, theta_f);

                if (pm.x() != pi.x() && pm.y() != pi.y() && pm.x() != pf.x() && pm.y() != pf.y()){
                  counter ++;
                  // Solve multipoint problem
                  std::vector<Configuration2> points = {pi, pm, pf};
                  std::vector<bool> fixedAngles = {true, false, true};
                  std::vector<double> curveParam = { k_max };
                  int discr = 90;
                  int refinements = 4;

                  std::vector<std::string> res = find_all_best(pi, pm, pf, fixedAngles, curveParam, discr, refinements);

                  if (res.size() > 0){
                    actual_counter ++;
                    // Write data to file
                    file << std::setprecision(5) << k_max << " " << xi << " " << xm << " " << ym << " " << xf << " " << 
                            theta_i << " " << theta_f << " " << pm.th() << " " << res[0] << std::endl;
                  }
                }

                auto dtime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();

                // Print time
                auto part = tot_counter > 100 ? tot_counter/100 : 1;
                if (counter % part == 0) {
                  std::cout << 100.0 * counter / tot_counter << "% " << counter << " in " << dtime << "ms, avg " << (dtime/(1.0*counter)) << "ms" << std::endl;
                }
                counter ++;
              }
            }
          }
        }
      }
    }
  }

  std::cout << "Generated " << PrintScientificLargeInt(actual_counter) << " entries to " << filename << std::endl;

  if (coutbuf != nullptr){
    std::cout.rdbuf(coutbuf);
  }

  file.close();
 }