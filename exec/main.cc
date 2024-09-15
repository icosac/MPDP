/**
 * @file main.cc
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License Agero 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Main file for the Dubins and Reed-Shepp paths computation.
 */

#include <dubins.hh>
#include <dp.hh>
#include <timeperf.hh>
#include <tests.hh>

#include <iostream>
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


