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
  // Dubins curve = Dubins(Configuration2(1, 1, -m_pi), Configuration2(0, 0, -m_pi_2), {1});
  // std::ofstream file("Dubins.asy");

  // std::cout << "Length: " << curve.l() << std::endl;

  // initAsyFile(file);
  // file << "path p;" << std::endl;
  // curve.draw(file);
  // file.close();

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
    Dubins d(start, end, k);
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
  // generateDataset3PDPCircle(argc, argv);
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


