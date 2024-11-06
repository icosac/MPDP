// The tests are defined in examples/MPDP/MPDP.hh
#include <MPMD.hh>

//! Tests names
std::vector<std::string> testsNames = {
    "Kaya Example 1",
    "Kaya Example 2",
    "Kaya Example 3",
    "Kaya Example 4",
    "Omega",
    "Circuit"
};

//! Tests descriptions
std::vector<std::vector<Configuration2> > Tests = {
    kaya1, kaya2, kaya3, kaya4, omega, spa
};

//! The tests's curvature
std::vector<K_T> Ks = {3.0, 3.0, 5.0, 3.0, 3.0, 3.0};

//! The tests's lengths
std::vector<LEN_T> exampleLenghts={3.41557885807514871601142658619, 6.27803455030931356617429628386, 11.9162126542854860389297755319, 7.46756219733842652175326293218, 41.0725016438839318766440555919, 6988.66098639942993031581863761}; //the last length is SPA

//! The number of discretizations to tests
std::vector<uint> discrs = {4, 16, 90, 360};
//! The number of refinements to tests
std::vector<uint> refins = {1, 2, 4, 8, 16};

/**
 * @brief This functions runs all the examples and prints the results in a LaTeX-like table
 * 
 * @return 0 if everything is OK
 */
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

int generateDataset(){
    std::string filename = "MP_dubins_dataset.csv";
    std::ofstream file;
    file.open(filename);

    int n_points_max = 20; 
    int n_points_min = 3;
    int n_tests_per_k = 1e3;
    std::vector<K_T> Ks = {1.0, 3.0, 5.0};

    double x_min = 0.0; 
    double x_max = 10.0;
    double y_min = 0.0;
    double y_max = 10.0;

    int discr = 360;
    int refin = 4;

    // Define random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> n_points_dist(n_points_min, n_points_max);
    std::uniform_real_distribution<double> x_dist(x_min, x_max);
    std::uniform_real_distribution<double> y_dist(y_min, y_max);
    std::uniform_real_distribution<double> theta_dist(0.0, m_pi);

    std::cout << "Computing a total of " << Ks.size()*n_tests_per_k << " tests" << std::endl;
    int counter = 0;

    for (auto k : Ks){
        for (int i=0; i<n_tests_per_k; i++){
            int n_points = n_points_dist(gen);
            // std::cout << "Gen test file with " << n_points << " points and k=" << k << std::endl;
            std::vector<Configuration2> points;
            std::vector<bool> fixedAngles;
            for (int j=0; j<n_points; j++){
                real_type x = x_dist(gen);
                real_type y = y_dist(gen);
                if (j==0 || j==n_points-1){
                    points.push_back(Configuration2(x, y, theta_dist(gen)));
                    fixedAngles.push_back(true);
                }
                else{
                    points.push_back(Configuration2(x, y, ANGLE::FREE));
                    fixedAngles.push_back(false);
                }
            }

            std::vector<real_type> curveParam = {k};

            std::pair<LEN_T, std::vector<Angle>>ret=DP().solveDP(points, fixedAngles, curveParam, discr, refin);

            file << n_points << " ";
            file << k << " ";
            file << ret.first << " ";
            file << points[0].th() << " ";
            file << points.back().th() << " ";
            for (auto point : points){
                file << point.x() << " ";
                file << point.y() << " ";
            }
            file << std::endl;
            counter++;
            std::cout << "Completed: " << 100.0*counter/(Ks.size()*n_tests_per_k) << std::endl;
        }
    }

    file.close();
}