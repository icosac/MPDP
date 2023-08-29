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
//    if (testID!=0){continue;}
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
        std::pair<LEN_T, std::vector<Angle> >ret=DP::solveDP(CURVE_TYPE::DUBINS, points, fixedAngles, curveParam, DISCR, r);
        auto time1=tp.getTime();
        LEN_T ComLength=ret.first;
        std::vector<Angle> vtheta=ret.second;

        LEN_T Length = 0.0;
        for (unsigned int index=points.size()-1; index>0; index--){
          points[index-1].th(vtheta[index-1]);
          points[index].th(vtheta[index]);
          Dubins c(points[index-1], points[index], {Ks[testID]});
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

Configuration2 circleLine(double s, double kur, Configuration2 c) {
  double sigmaDir = 1,  sign = kur > 0? 1: -1;
  double xEnd, yEnd, thetaEnd;
  xEnd = c.x() - f(s, kur , mod2pi(c.th() + sigmaDir * m_pi));
  yEnd = c.y() - g(s, kur , mod2pi(c.th() + sigmaDir * m_pi));
  thetaEnd = mod2pi(c.th() + kur  * s);

  Configuration2 res(xEnd, yEnd, thetaEnd, kur );

  return res;
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
        bestMan = d1.type_to_string() + " " + d2.type_to_string();
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
        bestMan = d1.type_to_string() + " " + d2.type_to_string();
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
    std::uniform_real_distribution<> piDist(0, M_PI*2.0);
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
          << d.kmax() << " " << d.l() << " " << d.type_to_string() << " "
          << t << std::endl;
  }
  file.close();
  return 0;
}

int example() {
  // Set problem data
  std::vector<Configuration2> points = {
      Configuration2(0, 0, 0),
      Configuration2(1, 0, 0),
      Configuration2(1, 1, 0),
      Configuration2(0, 1, 0),
      Configuration2(0, 0, 0)
  };
  K_T kmax = 1;
  int refinements = 2;
  int discretizations = 360;
  std::vector<bool> fixedAngles = { true, false, false, false, true };
  std::vector<double> parameters = {kmax};

  // Solve the multi-point problem
  std::pair<LEN_T, std::vector<Angle> >ret=DP::solveDP(CURVE_TYPE::RS, points, fixedAngles, parameters, discretizations, refinements);

  std::cout << "Total length: " << ret.first << std::endl;

  // Set actual computed angles
  for (int i=0; i<ret.second.size(); ++i){
    points[i].th(ret.second[i]);
  }

  // Print final points
  std::cout << "Final points with angles:" << std::endl;
  for (auto point : points) {
    std::cout << point << std::endl;
  }

  return 0;
}

int main() {
  return
//  tentaclesFigDubins();
//  main3PMDBruteForce();
//  allexamples();
//  genDSDubinsP2P(true);
  example();
}
