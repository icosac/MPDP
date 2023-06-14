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
        std::pair<LEN_T, std::vector<Angle> >ret=DP::solveDP(CURVE_TYPE::RS, points, fixedAngles, curveParam, DISCR, r);
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

Configuration2 circleLine(double s, double kur, Configuration2 c) {
  double sigmaDir = 1,  sign = kur > 0? 1: -1;
  double xEnd, yEnd, thetaEnd;
  xEnd = c.x() - f(s, kur , mod2pi(c.th() + sigmaDir * m_pi));
  yEnd = c.y() - g(s, kur , mod2pi(c.th() + sigmaDir * m_pi));
  thetaEnd = mod2pi(c.th() + kur  * s);

  Configuration2 res(xEnd, yEnd, thetaEnd, kur );

  return res;
}


int tentaclesFigDubins() {
  // Define the data
  double Kmax = 3;
  double thi = - m_pi / 4;
  double thf = 0 * m_pi / 3;
  //Configuration2 CI(-1, 1, thi), CF(1, -0.5, thf);
  std::vector<real_type> curveParam = { Kmax };
  std::vector<double> debug = {};
  std::vector<double> xf = { 1, 1, 1, 0 };
  std::vector<double> yf = { 1, -1, 0, 1 };
  std::vector<double> thf_i = { 0, 0 , -m_pi / 2, 0 };
  std::vector<double> thf_f = { 1 * m_pi, m_pi, m_pi / 2, m_pi };
  int Npts = xf.size();

  int Nangles = 16;


  TimePerf time;
  time.start();
  //G2lib::AsyPlot plot("TentaclesDubins.asy", false, 10, 5, false);
  std::ofstream file("TentaclesDubins.asy");
  file << "import graph; \n include \"clothoidLib.asylib\";\n size(8cm, 8cm); " << std::endl;
  for (int i = 0; i < Npts; ++i) {
    // loop over final points
    std::cout << "Nangles = " << Nangles << std::endl;
    for (int j = 0; j < Nangles; ++j) {
      //loop over final angles
      thf = thf_i[i] + ((thf_f[i] - thf_i[i]) / ((double)Nangles)) * j;
      std::cout << "thf = " << thf << std::endl;
      Configuration2 CI(-1, 0, thi), CF(xf[i], yf[i], thf);
      // Solve the Reed-Shepp considering all the possible cases since curveParam only contains the curvature
      Dubins myRS(CI, CF, curveParam);

      //plot
      std::string str = "";
      file << "path p;" << std::endl;
      //plot.verbatim("path p;");
      int ii = 0;
      Configuration2 c = myRS.ci()[0];
      for (ii = 1; ii < 4; ++ii) {
        str = "p = clothoidPoints((" + std::to_string(c.x()) + "," + std::to_string(c.y()) + "), " + std::to_string(c.th())
              + "," + std::to_string(myRS.k(ii)) + ", 0, " + std::to_string(myRS.L(ii)) + ");";
        std::cout <<  " k= " << myRS.k(ii) << "     kmax = "<< myRS.kmax() << std::endl;
        file << str << std::endl;
        //plot.verbatim(str);
        str = "royalblue";

        //plot.verbatim("draw(p," + str + ");");

        file << "draw(p," << str << ");" << std::endl;

        //plot.dot(myRS.getX()[ii], myRS.getY()[ii], "red");
        file << "dot((" << c.x() << "," << c.y() << "), red);" << std::endl;

        c = circleLine(myRS.L(ii), myRS.k(ii), c);
      }
      file << "dot((" << myRS.ci()->x() << "," << myRS.ci()->y() << "), black);" << std::endl;
      file << "dot((" << myRS.cf()->x() << "," << myRS.cf()->y() << "), purple+3bp);" << std::endl;

      //plot.dot(myRS.getX()[0], myRS.getY()[0], "black");
      //plot.dot(myRS.getX()[ii], myRS.getY()[ii], "red");

    }
  }

  std::cout << "Took: " << time.getTime() << "ms" << std::endl;

  double lim = 1.5;
  //plot.displayAxes("$x$", "$y$", -lim, lim, -1, 1);

  return 0;
}

struct Best {
  LEN_T len, L0, L1, L2, L3, L4, L5;
  Angle angle;
  std::string man;
};

void test3PMD(const Configuration2& Pi, const Configuration2& Pf, Configuration2& Pm,
              size_t disc, size_t ref, double kmax, double halfInt,
              struct Best & best) {
  std::vector<double> testTh;
  for (int i = 0; i < disc; ++i) {
    testTh.push_back(2.0 * m_pi / disc * i);
  }

  for (int r = 0; r < ref; ++r) {
    for (double thm: testTh) {
      Pm.th(thm);
      Dubins d1(Pi, Pm, {kmax});
      Dubins d2(Pm, Pf, {kmax});
      LEN_T currLen = d1.l() + d2.l();
      if (best.len > currLen) {
        best.len = currLen;
        best.angle = thm;
        best.man = d1.type_to_string() + " " + d2.type_to_string();
        best.L0 = d1.s1();
        best.L1 = d1.s2();
        best.L2 = d1.s3();
        best.L3 = d2.s1();
        best.L4 = d2.s2();
        best.L5 = d2.s3();
      }
    }

    double thtmp = best.angle - halfInt;
    testTh.clear();
    testTh.push_back(thtmp);
    for (int j = 1; j < disc; ++j) {
      thtmp += halfInt * 2.0 / disc * j;
      testTh.push_back(thtmp);
    }
  }
}

int main3PMDBruteForce(){
  // Vecto containing the number of discretizations to move through
  std::vector<int> discs(14); std::iota(discs.begin(), discs.end(), 3);
  discs.push_back(90); discs.push_back(360);
  // Vector containing the number of refinements to try
  std::vector<int> refs = { 1, 2, 4, 16};

  // Number of repetitions for each test
  size_t repetitions = 33;

  Configuration2 Pi (0, 0, m_pi/2);
  Configuration2 Pf (1, 1, 0);
  Configuration2 Pm (0.2, 0, m_pi / 2+ m_pi/8);
  K_T kmax = 1;
  // Half the interval around the best found angle.
  // So let's say the \theta* = \pi, then the tested interval will be [3/5\pi, 7/5\pi]
  double halfInt = m_pi*2/5;

  // Struct that contains variable for analysis
  struct Best best;
  best.len = std::numeric_limits<LEN_T>::infinity();
  best.angle = 0.0;
  // While `best` will be updated for each test, these values will store the best combination
  struct Best bestest;
  bestest.len = std::numeric_limits<LEN_T>::infinity();
  bestest.angle = 0.0;
  bestest.man = "";
  double bestestTime = 0.0;
  size_t bestestDisc = 0;
  size_t bestestRef = 0;

  test3PMD(Pi, Pf, Pm, 360, 16, kmax, halfInt, best);
  LEN_T refLen = best.len;

  // Matrix that will contain the error, time tuples to print
  std::vector<std::vector<std::pair<double, double>>> values (discs.size());
  for (int i = 0; i < discs.size(); ++i) {
    values[i].resize(refs.size());
  }

  for (auto disc : discs) {
    for (auto ref : refs){
      std::cout << "testing disc: " << disc << " ref: " << ref << std::endl;
      best.len = std::numeric_limits<LEN_T>::infinity();
      best.angle = 0.0;

      TimePerf time1;
      time1.start();

      test3PMD(Pi, Pf, Pm, disc, ref, kmax, halfInt, best);

      double elapsed = time1.getTime()/(double)(repetitions);

      if (best.len < bestest.len || (best.len == bestest.len && elapsed < bestestTime)) {
        bestest.len = best.len;
        bestest.angle = best.angle;
        bestest.man = best.man;
        bestestTime = elapsed;
        bestestDisc = disc;
        bestestRef = ref;
      }

      int distanceId = std::distance(discs.begin(), std::find(discs.begin(), discs.end(), disc));
      int refId = std::distance(refs.begin(), std::find(refs.begin(), refs.end(), ref));
      values[distanceId][refId] = {best.len, elapsed};
    }
  }

  // Print LaTeX-like table
  std::cout << "k&m&L&Time\\\\\n\\hline"<< std::endl;
  for (int i = 0; i < discs.size(); ++i) {
    for (int j = 0; j < refs.size(); ++j) {
      std::cout << std::setw(3) << discs[i] << "&";
      std::cout << std::setw(2) << refs[j] << "&";
      PrintScientific1D(values[i][j].first - refLen);
//      std::cout << std::fixed << std::setprecision(8) << std::setw(11) << values[i][j].first - refLen << "&";
      std::cout << "&" << std::fixed << std::setprecision(3) << std::setw(6) << values[i][j].second << "\\\\" << std::endl;
    }
  }

  // Print info for best global values
  printf("Shortest path with angle %.8f (%.2fpi) and total length %.8f with man %s\n", bestest.angle, (bestest.angle/ m_pi), bestest.len, bestest.man.c_str());
  std::cout << "L0 = " << best.L0 << std::endl;
  std::cout << "L1 = " << best.L1 << std::endl;
  std::cout << "L2 = " << best.L2 << std::endl;
  std::cout << "L3 = " << best.L3 << std::endl;
  std::cout << "L4 = " << best.L4 << std::endl;
  std::cout << "L5 = " << best.L5 << std::endl;
  std::cout << "Took: " << bestestTime << "ms with disc: " << bestestDisc << " and ref: " << bestestRef << std::endl;

  // Plot
  std::ofstream file("3pointDubins.asy");
  file << "import graph; \n include \"clothoidLib.asylib\";\n size(8cm, 8cm); " << std::endl;
  file << "path p;" << std::endl;

  Pm.th(best.angle);
  std::vector<real_type> curveParam = { kmax };
  Dubins myDub1(Pi, Pm, curveParam);
  Dubins myDub2(Pm, Pf, curveParam);

  myDub1.draw(file);
  myDub2.draw(file);

  file.close();

  system("asy -f pdf 3pointDubins.asy");

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
  system("pause");
#endif

  return 0;
}

int main() {
  //tentaclesFigDubins();
  main3PMDBruteForce();
}