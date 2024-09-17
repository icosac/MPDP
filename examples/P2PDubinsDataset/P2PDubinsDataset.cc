#include "P2PDubinsDataset.hh"

/**
 * @brief Function to generate a random dataset of Dubins for testing.
 * @param fix_initial_pos Whether the initial position should be fixed to (0,0,0) or not
 * @return
 */
int genDSDubinsP2P(bool fix_initial_pos){
  std::seed_seq seed{1,2,3,4,5,6,7,8,9};
  std::mt19937 eng(seed);

  std::string filename = "DSdubinsP2P6.csv";
  if (!fix_initial_pos) {
    filename = "DSdubinsP2P6Rand.csv";
  }
  std::ofstream file(filename.c_str());
  file << "x0 y0 th0 x1 y1 th1 kmax l type time" << std::endl;

  K_T kmax = 1;

  for (int nTests = 1000000; nTests > 0; nTests--){
    std::uniform_real_distribution<> cooDist(0, 1000.0);
    std::uniform_real_distribution<> piDist(0, m_pi*2.0);
    Configuration2 start (0, 0, 0);
    if (!fix_initial_pos) {
      start = Configuration2(cooDist(eng), cooDist(eng), piDist(eng));
    }
    Configuration2 end (cooDist(eng), cooDist(eng), piDist(eng));

    TimePerf time;
    time.start();
    Dubins d(start, end, kmax);
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
