#pragma once
#include <iostream>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "dubins.hh"


#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#elif __has_include("Eigen/Dense")
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#else
#error "Eigen library is required to solve the coefficient system."
#endif

using namespace mpdp;
using namespace mpdp::cpu;


const std::map<std::string, std::tuple<int, Dubins::D_TYPE, Dubins::D_TYPE>> P3DP_DICT = {
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

Eigen::MatrixXd find_coefficients_01 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_02 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_03 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_04 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_05 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_06 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_07 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_08 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_09 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_10 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_11 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_12 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_13 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_14 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_15 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_16 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_17 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_18 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);


Eigen::MatrixXd find_coefficients (
    size_t index,
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);


template<int N_COEFF, typename CoeffFunc>
std::vector<double> solve_man(
    CoeffFunc&& coeff_fn,
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance
) {
    static_assert(N_COEFF >= 2, "solve_man requires at least two coefficients.");

    Eigen::MatrixXd coefficients =
        std::forward<CoeffFunc>(coeff_fn)(xi, yi, thi, xm, ym, xf, yf, thf, r);
    const Eigen::Index available_coeffs = coefficients.size();
    if (available_coeffs < N_COEFF) {
        std::cerr << "Requested " << N_COEFF << " coefficients, but only "
                  << available_coeffs << " were provided." << std::endl;
        return {};
    }

    constexpr int degree = N_COEFF - 1;
    const double leading_coefficient = coefficients(degree);
    if (std::abs(leading_coefficient) < 1e-12) {
        std::cerr << "Leading coefficient is too small, cannot build companion matrix." << std::endl;
        return {};
    }

    using CompanionMatrix = Eigen::Matrix<double, degree, degree>;
    CompanionMatrix companion = CompanionMatrix::Zero();
    for (int i = 1; i < degree; ++i) {
        companion(i, i - 1) = 1.0;
    }
    for (int i = 0; i < degree; ++i) {
        companion(0, i) = -coefficients(degree - 1 - i) / leading_coefficient;
    }

    Eigen::EigenSolver<CompanionMatrix> eigen_solver(companion);
    const auto roots = eigen_solver.eigenvalues();

    std::vector<double> solutions;
    solutions.reserve(static_cast<std::size_t>(degree));

    const auto root_count = roots.size();
    for (Eigen::Index i = 0; i < root_count; ++i) {
        const auto& root = roots(i);
        if (std::abs(root.imag()) < imaginary_tolerance) {
            solutions.push_back(2.0 * std::atan(root.real()));
        }
    }
    
    return solutions;
}

std::vector<double> solve_man_01(
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance = 1e-8
);

std::vector<double> solve_man_07(
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance = 1e-8
);

std::vector<double> solve_man_10(
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance = 1e-8
);

std::vector<double> solve_man_11(
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance = 1e-8
);

std::vector<double> solve_man_12(
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance = 1e-8
);

std::vector<double> solve_man_13(
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance = 1e-8
);

std::vector<double> solve_man_14(
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance = 1e-8
);