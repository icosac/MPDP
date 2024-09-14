/**
 * @file main.cpp
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License Agero 3.0.
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

// #include<RSPredict.h>
// #include<SVMQuadraticPredict.h>
// #include<NNWideCompactPredict.h>

std::vector<Configuration2> example1 = {
    Configuration2(0,0,-2.0*m_pi/8.0),
    Configuration2(2,2,ANGLE::FREE),
    Configuration2(6,-1,ANGLE::FREE),
    Configuration2(8,1,2.0*m_pi/8.0)
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


int allexamples (){
  std::cout << "DISCR & ref & dl & t\\" << std::endl;
  for (uint testID=0; testID<Tests.size(); testID++){
    std::cout << "Test " << testID << std::endl;
    if (testID!=2){continue;}
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
    std::vector<real_type> curveParam={Ks[testID]};

    for (auto DISCR :  discrs){
//      if (DISCR!=4){continue;}
      for (auto r : refins){
//        if (r!=4){continue;}
        //std::cout << DISCR << " " << r << " ";
        TimePerf tp, tp1;
        std::vector<Configuration2>points=Tests[testID];

        tp.start();
        std::pair<LEN_T, std::vector<Angle> >ret=DP().solveDP(points, fixedAngles, curveParam, DISCR, r);
        auto time1=tp.getTime();
        LEN_T ComLength=ret.first;
        std::vector<Angle> vtheta=ret.second;

        LEN_T Length = 0.0;

        for (unsigned int idjijij=points.size()-1; idjijij>0; idjijij--){
          points[idjijij-1].th(vtheta[idjijij-1]);
          points[idjijij].th(vtheta[idjijij]);
          RS c(points[idjijij-1], points[idjijij], {Ks[testID]});
          // std::cout << c << std::endl;
          Length+=c.l();
        }
        // std::cout << Length << " " << exampleLenghts[testID] << std::endl;
        printf("%3d & %2d & ", DISCR, r);
        PrintScientific2D((ComLength-exampleLenghts[testID])*1000.0);
        printf(" & ");
        PrintScientific1D(time1);
        printf("     %.16f %.16f", ComLength, exampleLenghts[testID]);
        printf("\\\\\n");
      }
    }
    printf("\n\n\n\n");
  }

  return 0;
}


int main3PMDBruteForce(){
  // caso LRL-RLR
  //Configuration2 Pi (0, 0, m_pi);
  //Configuration2 Pf (1, 1, m_pi);
  //Configuration2 Pm (0.5, 0.5, 0);

  // caso 1: RSR-RLR
  //Configuration2 Pi(-1, 0, m_pi / 3);
  //Configuration2 Pf(1, 0, -m_pi / 2+m_pi/8);
  //Configuration2 Pm(0.5, 0.5, 0);

  // caso 2:
  Configuration2 Pi(-1, 0, m_pi / 3);
  Configuration2 Pf(1, 0,  m_pi / 2+m_pi/8);
  Configuration2 Pm(0.5, 0.5, 0);


  K_T kmax =1;
  std::cout << "Prova........." << std::endl;
  double thInc = 0.001; //std::numeric_limits<double>::epsilon();
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

/**
 * @brief Function to generate a random dataset of Dubins for testing.
 * @param fix_initial_pos Whether the initial position should be fixed to (0,0,0) or not
 * @return
 */
int genDSDubinsP2P(bool fix_initial_pos = false){
  std::seed_seq seed{1,2,3,4,5,6,7,8,9};
  std::mt19937 eng(seed);

  std::string filename = "DSdubinsP2P6.csv";
  if (!fix_initial_pos) {
    filename = "DSdubinsP2P6Rand.csv";
  }
  std::ofstream file(filename.c_str());
  file << "x0 y0 th0 x1 y1 th1 kmax l type time" << std::endl;

  K_T kmax = 1;
  int nTests = 1000000;

  for (; nTests > 0; nTests--){
    std::uniform_real_distribution<> cooDist(0, 1000.0);
    std::uniform_real_distribution<> piDist(0, m_pi*2.0);
    Configuration2 start (0, 0, 0);
    if (!fix_initial_pos) {
      start = Configuration2(cooDist(eng), cooDist(eng), piDist(eng));
    }
    Configuration2 end (cooDist(eng), cooDist(eng), piDist(eng));

    TimePerf time;
    time.start();
    Dubins d(start, end, {kmax});
    auto t = time.getTime();

    file  << std::setprecision(12)
          << start.x() << " " << start.y() << " " << start.th() << " "
          << end.x() << " " << end.y() << " " << end.th() << " "
          << d.kmax() << " " << d.l() << " " << d.man_to_string() << " "
          << t << std::endl;
  }
  file.close();
  return 0;
}

void test_times(){
  bool fix_initial_pos = false;
  std::seed_seq seed{1,2,3,4,5,6,7,8,9};
  std::mt19937 eng(seed);

  std::string filename = "DSdubinsP2P6.csv";
  if (!fix_initial_pos) {
    filename = "DSdubinsP2P6Rand.csv";
  }
  std::ofstream file(filename.c_str());
  file << "x0 y0 th0 x1 y1 th1 kmax l type time" << std::endl;

  K_T kmax = 1;
  int nTests = 1000000;

  for (; nTests > 0; nTests--){
    std::uniform_real_distribution<> cooDist(0, 1000.0);
    std::uniform_real_distribution<> piDist(0, m_pi*2.0);
    Configuration2 start (0, 0, 0);
    if (!fix_initial_pos) {
      start = Configuration2(cooDist(eng), cooDist(eng), piDist(eng));
    }
    Configuration2 end (cooDist(eng), cooDist(eng), piDist(eng));

    TimePerf time;
    time.start();
    Dubins d(start, end, {kmax});
    auto t = time.getTime();

    file  << std::setprecision(12)
          << start.x() << " " << start.y() << " " << start.th() << " "
          << end.x() << " " << end.y() << " " << end.th() << " "
          << d.kmax() << " " << d.l() << " " << d.man_to_string() << " "
          << t << std::endl;
  }
  file.close();
}

void findThose2Manoeuvres(std::string filename){
  std::ofstream file(filename.c_str());
  file << "x0 y0 th0 x1 y1 th1 kmax ";
  for (int i=0; i<52; i++){
    file << i << " ";
  }
  std::vector<double> THIs = {};
  std::vector<double> THFs = {};
  std::vector<double> Ks = {};


  for (double thi : THIs){
    for (double thf : THFs) {
      for (double k: Ks) {
        Configuration2 initPos = Configuration2(0, 0, thi);
        Configuration2 finalPos = Configuration2(1, 0, thf);
        RS myRS(initPos, finalPos, {k});
      }
    }
  }

  file.close();
}

int checkAnotherMan(RS myRS){
  int oldMan = myRS.getNman();
  int newMan = -1;
  LEN_T oldLen = myRS.l();
  LEN_T newLen = 0.0;

  for(int i=1; i<49; ++i){
    if (i!=2 && i!=4){
      RS newRS = RS(*myRS.ci(), *myRS.cf(), {myRS.getKmax(), (double)i});
      newRS.solve();
      if (std::abs(newRS.l()-oldLen) < 1e-8){
        newMan = i;
        newLen = myRS.l();
      }
    }
  }

//  if (newMan != 1 && newMan != 3) {
    std::cout << std::setprecision(12);
    std::cout << "Found alternative: "
              << newMan << " " << newLen << " to: "
              << oldMan << " " << oldLen << " diff: "
              << std::abs(newLen - oldLen) << std::endl;
//  }

  return newMan;
}

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

// struct PerfData{
//   double avgTime = 0.0;
//   double maxTime = 0.0;
//   double minTime = 100000000.0;
// };

// int predictRS(short modelType){
//   if (modelType == 1){
//     std::cout << "Calling predictRS with NNWide model" << std::endl;
//   }
//   else if(modelType == 2) {
//     std::cout << "Calling predictRS with SVMQuadratic model" << std::endl;
//   }
//   else{
//     std::cout << "Invalid model type" << std::endl;
//     return 1;
//   }

//   double idx[22];
//   double dv[3];
//   int tmp[2];

//   int nTests = 100;
//   int nRepet = 100;

//   double predTime = 0.0;
//   double fullTime = 0.0;
//   double direTime = 0.0;

//   std::uniform_real_distribution<double> randTh(0, m_2pi);
//   std::uniform_real_distribution<double> randKmax(0, 10.0);
//   std::default_random_engine re(13);

//   // Getting a random double value

//   for (int iTest = 0; iTest < nTests; iTest++) {
//     dv[0] = randTh(re);
//     dv[1] = randTh(re);
//     dv[2] = randKmax(re);

//     for (int iRep = 0; iRep < nRepet; iRep++) {
//       TimePerf time;
//       time.start();
//       RS myRS = RS(Configuration2(0.0, 0.0, dv[0]), Configuration2(1.0, 0.0, dv[1]), {dv[2]});
//       myRS.solve();
//       auto delta = time.getTime();
//       fullTime += delta;

//       time.start();
//       RSPredict(modelType, dv, idx, tmp);
//       delta = time.getTime();
//       predTime += delta;

//       time.start();
//       RS myRS2 = RS(Configuration2(0, 0, dv[0]), Configuration2(1, 0, dv[1]), {dv[2], (double) myRS.getNman()});
//       myRS2.solve();
//       delta = time.getTime();
//       direTime += delta;
//     }
//   }

//   std::cout << "AVG brute-force: " << fullTime/(nTests*nRepet) << std::endl;
//   std::cout << "AVG prediction: " << predTime/(nTests*nRepet) << std::endl;
//   std::cout << "AVG one-shot: " << direTime/(nTests*nRepet) << std::endl;

//   return 0;
// }

// std::vector<int> labels = {
//     1, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24
// };

// int SVMQuadraticRS(){
//   std::cout << "Calling SVMQuadraticRS" << std::endl;
//   double idx[22];
//   double dv[3];

//   int nTests = 100;
//   int nRepet = 100;

//   std::uniform_real_distribution<double> randTh(0, m_2pi);
//   std::uniform_real_distribution<double> randKmax(0, 10.0);
//   std::default_random_engine re(13);

//   PerfData predData; // Prediction data
//   PerfData bfData;   // Brute-force data
//   PerfData oneData;  // One-shot data

//   for (int iTest = 0; iTest < nTests; iTest++) {
//     dv[0] = randTh(re);
//     dv[1] = randTh(re);
//     dv[2] = randKmax(re);

//     for (int iRep = 0; iRep < nRepet; iRep++) {
//       TimePerf time;
//       RS myRS = RS(Configuration2(0.0, 0.0, dv[0]), Configuration2(1.0, 0.0, dv[1]), {dv[2]});
//       time.start();
//       myRS.solve();
//       auto delta = time.getTime();
//       bfData.avgTime += delta;
//       bfData.maxTime = std::max(bfData.maxTime, delta);
//       bfData.minTime = std::min(bfData.minTime, delta);

//       time.start();
//       SVMQuadraticPredict(dv, idx);
//       delta = time.getTime();
//       predData.avgTime += delta;
//       predData.maxTime = std::max(predData.maxTime, delta);
//       predData.minTime = std::min(predData.minTime, delta);

//       RS myRS2 = RS(Configuration2(0, 0, dv[0]), Configuration2(1, 0, dv[1]), {dv[2], idx[0]});
//       time.start();
//       myRS2.solve();
//       delta = time.getTime();

// //      if (myRS.l() != myRS2.l()){// throw exception saying that the two lengths are different and printing the two lengths
// //        std::cout << "Lengths are different: " << myRS.l() << " " << myRS2.l() << std::endl;
// //        throw std::runtime_error("Lengths are different");
// //      }

//       oneData.avgTime += delta;
//       oneData.maxTime = std::max(oneData.maxTime, delta);
//       oneData.minTime = std::min(oneData.minTime, delta);
//     }
//   }

//   std::cout << "AVG brute-force: " << bfData.avgTime/(nTests*nRepet)  << " MIN: " << bfData.minTime << " MAX: " <<  bfData.maxTime << std::endl;
//   std::cout << "AVG prediction: " << predData.avgTime/(nTests*nRepet) << " MIN: " << predData.minTime << " MAX: " <<  predData.maxTime << std::endl;
//   std::cout << "AVG one-shot: " << oneData.avgTime/(nTests*nRepet)    << " MIN: " << oneData.minTime << " MAX: " <<  oneData.maxTime << std::endl;

//   return 0;
// }

// int NNWideRS(){
//   std::cout << "Calling NNWideRS" << std::endl;
//   double idx[22];
//   double dv[3];
//   int tmp[2];

//   int nTests = 100;
//   int nRepet = 100;

//   std::uniform_real_distribution<double> randTh(0, m_2pi);
//   std::uniform_real_distribution<double> randKmax(0, 10.0);
//   std::default_random_engine re(13);

//   PerfData predData; // Prediction data
//   PerfData bfData;   // Brute-force data
//   PerfData oneData;  // One-shot data

//   for (int iTest = 0; iTest < nTests; iTest++) {
//     dv[0] = randTh(re);
//     dv[1] = randTh(re);
//     dv[2] = randKmax(re);

//     for (int iRep = 0; iRep < nRepet; iRep++) {
//       TimePerf time;
//       time.start();
//       double bestL = 1e100;
//       int bestMan = 0;
//       for (int i = 1; i < 49; i++) {
//         if(i == 2 || i == 4) { continue; }
//         RS myRS = RS(Configuration2(0.0, 0.0, dv[0]), Configuration2(1.0, 0.0, dv[1]), {dv[2], (double)(i)});
//         myRS.solve();
//         if (myRS.l() < bestL) {
//           bestL = myRS.l();
//           bestMan = i;
//         }
//       }
//       auto delta = time.getTime();
//       bfData.avgTime += delta;
//       bfData.maxTime = std::max(bfData.maxTime, delta);
//       bfData.minTime = std::min(bfData.minTime, delta);

//       time.start();
//       NNWideCompactPredict(dv, idx, tmp);
//       delta = time.getTime();
//       predData.avgTime += delta;
//       predData.maxTime = std::max(predData.maxTime, delta);
//       predData.minTime = std::min(predData.minTime, delta);

//       RS myRS2 = RS(Configuration2(0, 0, dv[0]), Configuration2(1, 0, dv[1]), {dv[2], idx[0]});
//       time.start();
//       myRS2.solve();
//       delta = time.getTime();
//       oneData.avgTime += delta;
//       oneData.maxTime = std::max(oneData.maxTime, delta);
//       oneData.minTime = std::min(oneData.minTime, delta);

// //      if (myRS.l() != myRS2.l()){// throw exception saying that the two lengths are different and printing the two lengths
// //        std::cout << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;
// //        std::cout << myRS.getNman() << " " << idx[0] << " " << myRS.getManTypeStr() << " " << myRS.getSegmentsData()[0].l << " " << myRS.getSegmentsData()[1].l << " " << myRS.getSegmentsData()[2].l << " " << myRS.getSegmentsData()[3].l << " " << std::endl;
// //        std::cout << "Lengths are different: " << myRS.l() << " " << myRS2.l() << std::endl;
// //        throw std::runtime_error("Lengths are different");
// //      }
//     }
//   }

//   std::cout << "AVG brute-force: " << bfData.avgTime/(nTests*nRepet)  << " MIN: " << bfData.minTime << " MAX: " <<  bfData.maxTime << std::endl;
//   std::cout << "AVG prediction: " << predData.avgTime/(nTests*nRepet) << " MIN: " << predData.minTime << " MAX: " <<  predData.maxTime << std::endl;
//   std::cout << "AVG one-shot: " << oneData.avgTime/(nTests*nRepet)    << " MIN: " << oneData.minTime << " MAX: " <<  oneData.maxTime << std::endl;

//   return 0;
// }


void drawRS(){
//  RS curve = RS(Configuration2(0, 0, m_pi_2), Configuration2(1, 1, 0), {1});
//  curve.solve();
//  std::ofstream file("RS.asy");
  Dubins curve = Dubins(Configuration2(1, 1, -m_pi), Configuration2(0, 0, -m_pi_2), {1});
  std::ofstream file("Dubins.asy");

  std::cout << "Length: " << curve.l() << std::endl;

  initAsyFile(file);
  file << "path p;" << std::endl;
  curve.draw(file);
  file.close();

//   system("asy -f pdf RS.asy && xdg-open RS.pdf");
}

int testDubins(){
  std::ifstream file2("/home/enrico/Projects/mpdp/tmp.txt");
  if (!file2.is_open()) {
    std::cout << "Error opening file" << std::endl;
    return 1;
  }
  std::string line;
  double sx, sy, sth, ex, ey, eth, k, l, t;
  double totDiff = 0.0, maxPosDiff = 0.0, maxNegDiff = 0.0;
  int i = 0;
  while (file2 >> sx >> sy >> sth >> ex >> ey >> eth >> k >> l >> t) {
    Configuration2 start(sx, sy, sth);
    Configuration2 end(ex, ey, eth);

    std::vector<double> tmp = {k};

    TimePerf tp; tp.start();
    Dubins d(start, end, {k});
//    Dubins d = DP::solveP2P<Dubins>(start, end, tmp);
    auto time1=tp.getTime();

    totDiff += time1 - t;
    maxPosDiff = std::max(maxPosDiff, time1 - t);
    maxNegDiff = std::min(maxNegDiff, time1 - t);

    // If the two lengths differs for more than 1e-10 print an error with the same precision
    if (std::abs(d.l() - l) > 1e-9) {
      std::cout << "Error at line " << i << " with values: " << std::setprecision(12) << d.l() << " " << l << std::endl;
    }
    i++;
  }
  std::cout << "Total difference: " << totDiff << std::endl;
  std::cout << "Average difference: " << totDiff / i << std::endl; // should be "0
  std::cout << "Max positive difference: " << maxPosDiff << std::endl;
  std::cout << "Max negative difference: " << maxNegDiff << std::endl;

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
#include <map>
#include <tuple>

static std::map<std::string, std::tuple<int, Dubins::D_TYPE, Dubins::D_TYPE>> P3DP_DICT = {
  {"RLRRLR", {1, Dubins::D_TYPE::RLR, Dubins::D_TYPE::RLR}},
  {"LRLLRL", {2, Dubins::D_TYPE::LRL, Dubins::D_TYPE::LRL}},
  {"RLRRSR", {3, Dubins::D_TYPE::RLR, Dubins::D_TYPE::RSR}},
  {"RLRRSL", {4, Dubins::D_TYPE::RLR, Dubins::D_TYPE::RSL}},
  {"LRLLSL", {5, Dubins::D_TYPE::LRL, Dubins::D_TYPE::LSL}},
  {"LRLLSR", {6, Dubins::D_TYPE::LRL, Dubins::D_TYPE::LSR}},
  {"RSRRLR", {7, Dubins::D_TYPE::RSR, Dubins::D_TYPE::RLR}},
  {"LSRRLR", {8, Dubins::D_TYPE::LSR, Dubins::D_TYPE::RLR}},
  {"RSLLRL", {9, Dubins::D_TYPE::RSL, Dubins::D_TYPE::LRL}},
  {"LSLLRL", {10, Dubins::D_TYPE::LSL, Dubins::D_TYPE::LRL}},
  {"RSRRSR", {11, Dubins::D_TYPE::RSR, Dubins::D_TYPE::RSR}},
  {"LSRRSR", {12, Dubins::D_TYPE::LSR, Dubins::D_TYPE::RSR}},
  {"RSRRSL", {13, Dubins::D_TYPE::RSR, Dubins::D_TYPE::RSL}},
  {"LSRRSL", {14, Dubins::D_TYPE::LSR, Dubins::D_TYPE::RSL}},
  {"LSLLSL", {15, Dubins::D_TYPE::LSL, Dubins::D_TYPE::LSL}},
  {"RSLLSL", {16, Dubins::D_TYPE::RSL, Dubins::D_TYPE::LSL}},
  {"LSLLSR", {17, Dubins::D_TYPE::LSL, Dubins::D_TYPE::LSR}},
  {"RSLLSR", {18, Dubins::D_TYPE::RSL, Dubins::D_TYPE::LSR}}
};


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
  std::string filename_base = "3PDS_" + std::to_string(angle_discr) + "_" + std::to_string(kmax_min) + "_"  + std::to_string(kmax_max) + "_"  + std::to_string(k_discr);
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
            if (counter % (tot_counter/100000) == 0) {
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

#include <set>

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

int main(int argc, char** argv){
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Dubins stuff
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// testDubins();
	// genDSDubinsP2P(true);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// MPDP Examples
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// allexamples();
	// tentaclesFigDubins();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 3 points Dubins stuff
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// main3PMDBruteForce();
  // main3PDP();
  generateDataset3PDPCircle(argc, argv);
  // generateDataset3PDPRect(argc, argv);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// RS stuff
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//  generateDatasetRS();
	//  NNWideRS();
	//  SVMQuadraticRS();
	// drawRS();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  return 0;
}


