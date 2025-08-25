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
    kaya1, kaya2, kaya3, kaya4 //, omega, spa
};

//! The tests's curvature
std::vector<K_T> Ks = {3.0, 3.0, 5.0, 3.0, 3.0, 3.0};
// std::vector<K_T> Ks = {3.0/2.0, 3.0/2.0, 5.0/2.0, 3.0/2.0, 3.0/2.0, 3.0/2.0};

//! The tests's lengths
std::vector<LEN_T> exampleLenghts={3.41557885807514871601142658619, 6.27803455030931356617429628386, 11.9162126542854860389297755319, 7.46756219733842652175326293218, 41.0725016438839318766440555919, 6988.66098639942993031581863761}; //the last length is SPA

//! The number of discretizations to tests
// std::vector<uint> discrs = {4, 16, 90, 360};
std::vector<uint> discrs = {360};
//! The number of refinements to tests
// std::vector<uint> refins = {1, 2, 4, 8, 16};
std::vector<uint> refins = {4};

/**
 * @brief This functions runs all the examples and prints the results in a LaTeX-like table
 *
 * @return 0 if everything is OK
 */
int allexamples (){
  std::cout << "DISCR & ref & dl & t\\" << std::endl;
  for (uint testID=0; testID<Tests.size(); testID++){
    // if (testID!=3) continue;
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
      for (auto r : refins){
        TimePerf tp, tp1;
        std::vector<Configuration2>points=Tests[testID];

        tp.start();
        std::pair<LEN_T, std::vector<Angle> >ret=DP().solveDP(points, fixedAngles, curveParam, DISCR, r);
        auto time1=tp.getTime();
        LEN_T ComLength=ret.first;
        std::vector<Angle> vtheta=ret.second;

        LEN_T Length = 0.0;

        std::ofstream file;
        file.open("dubins_" + std::to_string(testID) + ".asy");
        for (unsigned int point_id = 0; point_id < points.size() - 1; point_id++) {
          Dubins d(points[point_id].x(), points[point_id].y(), vtheta[point_id],
                   points[point_id + 1].x(), points[point_id + 1].y(), vtheta[point_id + 1],
                   curveParam);
          Length += d.l();
          d.draw(file, "Dubins" + std::to_string(point_id), 800, 600, true, false, point_id==0);
        }
        file.close();
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

    // Define random number generator with seed 13
    std::mt19937 gen(13);
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

            TimePerf time;
            time.start();
            std::pair<LEN_T, std::vector<Angle>>ret=DP().solveDP(points, fixedAngles, curveParam, discr, refin);
            auto dt = time.getTime();

            file << n_points << " ";
            file << k << " ";
            file << std::setprecision(10) << dt << " ";
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

    return 1;
}

int test_from_file(const std::string& filename, bool set_th0, bool set_thf, std::string fig_filename){
		std::cout << "Should I be reading th0 and thf? " << set_th0 << " " << set_thf << std::endl;

		std::ifstream file(filename);
		if (!file.is_open()) {
			std::cout << "Error opening file" << std::endl;
			return 1;
		}

		int n_lines = 0;
		file >> n_lines;

		std::vector<Configuration2> points;

		double x, y;
		uint64_t counter = 0;
		while(file >> x >> y){
			std::cout << counter << " " << x << " " << y <<  std::endl;
			if (counter == 0 && set_th0){
				std::cout << "Reading th0" << std::endl;
				double th0;
				file >> th0;
				points.push_back(Configuration2(x, y, th0));
			}
			else if (counter == n_lines - 1 && set_thf){
				std::cout << "Reading th1" << std::endl;
				double thf;
				file >> thf;
				points.push_back(Configuration2(x, y, thf));
			}
			else { points.push_back (Configuration2 (x, y, ANGLE::FREE)); }
			counter++;
		}

		if (counter != n_lines){
			std::cout << "Error reading file" << std::endl;
			return 1;
		}

		if (!set_th0)
		{
			points[0].th (atan2 (points[1].y() - points[0].y(), points[1].x() - points[0].x()));
		}
		if (!set_thf)
		{
			points.back().th(atan2(points[points.size()-1].y()-points[points.size()-2].y(), points[points.size()-1].x()-points[points.size()-2].x()));
		}

		std::vector<bool> fixedAngles;
		for (uint i=0; i<points.size(); i++){
			std::cout << points[i] << std::endl;
			if (i==0 || i==points.size()-1) {
				fixedAngles.push_back(true);
			}
			else {
				fixedAngles.push_back(false);
			}
		}

		K_T kmax = 1.0;
		std::vector<real_type> curveParam={kmax};

		int discr = 360;
		int refin = 4;

		TimePerf tp;
		tp.start();
		std::pair<LEN_T, std::vector<Angle> >ret=DP().solveDP(points, fixedAngles, curveParam, discr, refin);
		auto time1=tp.getTime();
		LEN_T len=ret.first;

		for (auto angle : ret.second){
			std::cout << angle << std::endl;
		}

		std::cout << "Computed path in " << time1 << "ms " << len << std::endl;

		for (size_t i = 0; i < points.size(); i++){
			points[i].th(ret.second[i]);
		}

		std::cout << std::endl;
		std::ofstream draw_file(fig_filename);
		for (size_t i = 0; i < points.size()-1; i++){
			Dubins d(points[i], points[i+1], {kmax});
			std::cout << "\n=============\n" << d;
			std::cout << "\n-------------\n";
			for (size_t j = 1; j < 4; j++){
				std::cout << d.to_string_piece(j) << std::endl;
			}
			d.draw(draw_file, std::to_string(i), 500, 500, false, false, i==0);
		}
		draw_file.close();

		return 0;
}