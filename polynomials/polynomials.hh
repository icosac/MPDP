#pragma once
#include <iostream>

#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#elif __has_include("Eigen/Dense")
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#else
#error "Eigen library is required to solve the coefficient system."
#endif

Eigen::MatrixXd find_coefficients_man_01 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_02 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_03 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_04 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_05 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_06 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_07 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_08 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_09 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_10 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_11 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_12 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_13 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_14 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_15 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_16 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_17 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);

Eigen::MatrixXd find_coefficients_man_18 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);


Eigen::MatrixXd find_coefficients_man (
    size_t index,
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
);