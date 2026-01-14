#include <polynomials.hh>

#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>

#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#elif __has_include("Eigen/Dense")
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#else
#error "Eigen library is required to solve the coefficient system."
#endif

Eigen::MatrixXd find_coefficients_11 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
){    
    double dxi = -r*sin(thi)-xi+xm;
    double dyi = r*cos(thi)-yi+ym;
    double dxf = r*sin(thf)+xf-xm; 
    double dyf = -r*cos(thf)+yf-ym;

    double t2 = dxf * dyi;
    double t3 = r * dxf;
    double t4 = dxi * dyf;
    double t5 = dxi * r;
    double t6 = t2 - t3 + t4 + t5;
    double t7 = t2 - t3 - t4 - t5;
    double t9 = dyf * r;
    double t10 = dyi * r;
    double t12 = -2 * t9 - 2 * t10;
    double t15 = 4 * dxf * dxi;
    double t17 = 4 * dyf * dyi;
    double t18 = 2 * t9;
    double t19 = 2 * t10;
    double t20 = -t15 + t17 - t18 + t19;
    double t22 = t2 + t3 - t4 + t5;
    double t26 = -6 * t2 - 6 * t4;
    double t30 = t15 - t17 - t18 + t19;
    double t34 = t2 + t3 + t4 - t5;

    Eigen::Matrix<double, 7, 1> coefficients;
    coefficients << -t6 * t7,
                    -t12 * t6 - t20 * t7,
                    -t12 * t20 - t22 * t6 - t26 * t7,
                    -t12 * t26 - t20 * t22 - t30 * t7,
                    -t12 * t30 - t22 * t26 - t34 * t7,
                    -t12 * t34 - t22 * t30,
                    -t22 * t34;

    return coefficients;
}


std::vector<double> solve_man_11(
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance
) {
    return solve_man<7>(
        find_coefficients_11,
        xi, yi, thi,
        xm, ym,
        xf, yf, thf,
        r, imaginary_tolerance
    );
}
