#include "polynomials.hh"

#include <stdexcept>
#include <iostream>
#include <chrono>


Eigen::MatrixXd find_coefficients (
    size_t index,
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
){
    switch(index){
        case 1:
            return find_coefficients_01(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 2:
            return find_coefficients_02(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 3:
            return find_coefficients_03(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 4:
            return find_coefficients_04(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 5:
            return find_coefficients_05(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 6:
            return find_coefficients_06(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 7:
            return find_coefficients_07(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 8:
            return find_coefficients_08(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 9:
            return find_coefficients_09(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 10:
            return find_coefficients_10(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 11:
            return find_coefficients_11(xi, yi, thi, xm, ym, xf, yf, thf, r);
        // case 12:
        //     return find_coefficients_12(xi, yi, thi, xm, ym, xf, yf, thf, r);
        // case 13:
        //     return find_coefficients_13(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 14:
            return find_coefficients_14(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 15:
            return find_coefficients_15(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 16:
            return find_coefficients_16(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 17:
            return find_coefficients_17(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 18:
            return find_coefficients_18(xi, yi, thi, xm, ym, xf, yf, thf, r);
        default:
            throw std::invalid_argument("Invalid index for find_coefficients");
    }
}


std::pair<Dubins, Dubins>
find_solution_three_points(
	Configuration2& pi,
	Configuration2& pm,
	Configuration2& pf,
	K_T kmax,
	std::string man_name
){
    const auto& tuple_dubins = P3DP_DICT.at(man_name);
    auto type_dub1 = std::get<1>(tuple_dubins);
    auto type_dub2 = std::get<2>(tuple_dubins);

    Dubins dub1 = Dubins(pi, pm, {kmax}, type_dub1);
    Dubins dub2 = Dubins(pm, pf, {kmax}, type_dub2);
    return std::make_pair(dub1, dub2);
}

std::pair<Dubins, Dubins>
find_solution_three_points(
	Configuration2& pi,
	Configuration2& pm,
	Configuration2& pf,
	K_T kmax,
	size_t man_id
){
  std::string man_name = "";
  // Loop through the keys and values of P3DP_DICT
  for (const auto& [word, data] : P3DP_DICT) {
    if (std::get<0>(data) == man_id) {
      man_name = word;
      break;
    }
  }

 return find_solution_three_points(pi, pm, pf, kmax, man_name);
}


std::pair<double, double> 
find_best_angle(
  Configuration2 ci, 
  Configuration2 cm, 
  Configuration2 cf, 
  double kmax,
  std::vector<double> angles,
  int man_id
){
  double best_th = 0.0;
  double best_len = 1e8;
  
  for (size_t i = 0; i<angles.size(); i++){
    cm.th(angles[i]);

    try{
        auto dubs = find_solution_three_points(ci, cm, cf, kmax, man_id);
        double len = dubs.first.l() + dubs.second.l();

        if (len < best_len){
            best_len = len;
            best_th = angles[i];
        }
    } catch (const std::runtime_error & e){
        std::cerr << "Could not compute man for angle: " << angles[i] << " " << e.what() << std::endl;
    }
  }

  return std::make_pair(best_th, best_len);
}


int main() {
    double xi = 1;
    double yi = 0;
    double r = 1/2.0;
    double kmax = 2.0;
    double thi = 2.0/3.0*M_PI;
    double thf = -5.0/9.0*M_PI;
    double xm = cos(5.0/9.0*M_PI);
    double ym = sin(5.0/9.0*M_PI);
    double xf = cos(2.0/3.0*M_PI);
    double yf = sin(2.0/3.0*M_PI);

    // Print the date and the time
    auto today = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(today);
    std::cout << "Computation started at " << std::ctime(&now_c) << std::endl;

    Configuration2 ci = Configuration2(xi, yi, thi);
    Configuration2 cm = Configuration2(xm, ym, 0.0);
    Configuration2 cf = Configuration2(xf, yf, thf);

    const size_t N_TESTS = 1000;
    
    // (1(6), 2*7(8), 2*7(8), 7(8), 10(8), 11(6), 12(12), 13(12), 14(12))*2

    auto start_01 = std::chrono::high_resolution_clock::now();
    std::vector<double> th_01 = {};
    std::pair<double, double> th_len;
    for (size_t i = 0; i < N_TESTS; i++){
        th_01 = solve_man_01(xi, yi, thi, xm, ym, xf, yf, thf, r);
        th_len = find_best_angle(ci, cm, cf, kmax, th_01, 1);
    }
    auto end_01 = std::chrono::high_resolution_clock::now();
    auto elapsed_01_us = std::chrono::duration_cast<std::chrono::microseconds>(end_01 - start_01);

    std::cout << "Angle for man 1: " << std::endl;
    for (auto th : th_01) {
        std::cout << "\t" << th << std::endl;
    }
    std::cout << "Best angle for man 1 is " << th_len.first << " with len: " << th_len.second << " in " //milliseconds 
                << elapsed_01_us.count() / (double)N_TESTS << "us." << std::endl << std::endl;

    //////////////////////////////////////////////////////////////////////////////////////

    auto start_07 = std::chrono::high_resolution_clock::now();
    std::vector<double> th_07 = {};
    for (size_t i = 0; i < N_TESTS; i++){
        th_07 = solve_man_07(xi, yi, thi, xm, ym, xf, yf, thf, r);
        th_len = find_best_angle(ci, cm, cf, kmax, th_07, 7);
    }
    auto end_07 = std::chrono::high_resolution_clock::now();
    auto elapsed_07_us = std::chrono::duration_cast<std::chrono::microseconds>(end_07 - start_07);

    std::cout << "Angle for man 7: " << std::endl;
    for (auto th : th_07) {
        std::cout << "\t" << th << std::endl;
    }
    std::cout << "Best angle for man 7 is " << th_len.first << " with len: " << th_len.second << " in "
              << elapsed_07_us.count() / (double)N_TESTS << "us." << std::endl << std::endl;

    //////////////////////////////////////////////////////////////////////////////////////

    auto start_10 = std::chrono::high_resolution_clock::now();
    std::vector<double> th_10 = {};
    for (size_t i = 0; i < N_TESTS; i++){
        th_10 = solve_man_10(xi, yi, thi, xm, ym, xf, yf, thf, r);
        th_len = find_best_angle(ci, cm, cf, kmax, th_10, 10);
    }
    auto end_10 = std::chrono::high_resolution_clock::now();
    auto elapsed_10_us = std::chrono::duration_cast<std::chrono::microseconds>(end_10 - start_10);

    std::cout << "Angle for man 10: " << std::endl;
    for (auto th : th_10) {
        std::cout << "\t" << th << std::endl;
    }
    std::cout << "Best angle for man 10 is " << th_len.first << " with len: " << th_len.second << " in "
              << elapsed_10_us.count() / (double)N_TESTS << "us." << std::endl << std::endl;

    //////////////////////////////////////////////////////////////////////////////////////

    auto start_11 = std::chrono::high_resolution_clock::now();
    std::vector<double> th_11 = {};
    for (size_t i = 0; i < N_TESTS; i++){
        th_11 = solve_man_11(xi, yi, thi, xm, ym, xf, yf, thf, r);
        th_len = find_best_angle(ci, cm, cf, kmax, th_11, 11);
    }
    auto end_11 = std::chrono::high_resolution_clock::now();
    auto elapsed_11_us = std::chrono::duration_cast<std::chrono::microseconds>(end_11 - start_11);

    std::cout << "Angle for man 11: " << std::endl;
    for (auto th : th_11) {
        std::cout << "\t" << th << std::endl;
    }
    std::cout << "Best angle for man 11 is " << th_len.first << " with len: " << th_len.second << " in "
              << elapsed_11_us.count() / (double)N_TESTS << "us." << std::endl << std::endl;

    //////////////////////////////////////////////////////////////////////////////////////

    auto start_12 = std::chrono::high_resolution_clock::now();
    std::vector<double> th_12 = {};
    for (size_t i = 0; i < N_TESTS; i++){
        th_12 = solve_man_12(xi, yi, thi, xm, ym, xf, yf, thf, r);
        th_len = find_best_angle(ci, cm, cf, kmax, th_12, 12);
    }
    auto end_12 = std::chrono::high_resolution_clock::now();
    auto elapsed_12_us = std::chrono::duration_cast<std::chrono::microseconds>(end_12 - start_12);

    std::cout << "Angle for man 12: " << std::endl;
    for (auto th : th_12) {
        std::cout << "\t" << th << std::endl;
    }
    std::cout << "Best angle for man 12 is " << th_len.first << " with len: " << th_len.second << " in "
              << elapsed_12_us.count() / (double)N_TESTS << "us." << std::endl << std::endl;

    //////////////////////////////////////////////////////////////////////////////////////

    auto start_13 = std::chrono::high_resolution_clock::now();
    std::vector<double> th_13 = {};
    for (size_t i = 0; i < N_TESTS; i++){
        th_13 = solve_man_13(xi, yi, thi, xm, ym, xf, yf, thf, r);
        th_len = find_best_angle(ci, cm, cf, kmax, th_13, 13);
    }
    auto end_13 = std::chrono::high_resolution_clock::now();
    auto elapsed_13_us = std::chrono::duration_cast<std::chrono::microseconds>(end_13 - start_13);

    std::cout << "Angle for man 13: " << std::endl;
    for (auto th : th_13) {
        std::cout << "\t" << th << std::endl;
    }
    std::cout << "Best angle for man 13 is " << th_len.first << " with len: " << th_len.second << " in "
              << elapsed_13_us.count() / (double)N_TESTS << "us." << std::endl << std::endl;

    //////////////////////////////////////////////////////////////////////////////////////

    auto start_14 = std::chrono::high_resolution_clock::now();
    std::vector<double> th_14 = {};
    for (size_t i = 0; i < N_TESTS; i++){
        th_14 = solve_man_14(xi, yi, thi, xm, ym, xf, yf, thf, r);
        th_len = find_best_angle(ci, cm, cf, kmax, th_14, 14);
    }
    auto end_14 = std::chrono::high_resolution_clock::now();
    auto elapsed_14_us = std::chrono::duration_cast<std::chrono::microseconds>(end_14 - start_14);

    std::cout << "Angle for man 14: " << std::endl;
    for (auto th : th_14) {
        std::cout << "\t" << th << std::endl;
    }
    std::cout << "Best angle for man 14 is " << th_len.first << " with len: " << th_len.second << " in "
              << elapsed_14_us.count() / (double)N_TESTS << "us." << std::endl << std::endl;


    double total_avg_time = elapsed_01_us.count() / (double)N_TESTS +
                            5 * (elapsed_07_us.count() / (double)N_TESTS) +
                            elapsed_10_us.count() / (double)N_TESTS +
                            elapsed_11_us.count() / (double)N_TESTS +
                            elapsed_12_us.count() / (double)N_TESTS +
                            elapsed_13_us.count() / (double)N_TESTS +
                            elapsed_14_us.count() / (double)N_TESTS;
    total_avg_time *= 2.0;
    std::cout << "Total average time for all manouvres: " << total_avg_time << "us." <<  std::endl;

    return 0;
}