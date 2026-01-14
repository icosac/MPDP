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


Eigen::MatrixXd find_coefficients_man_14_p4 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
){
    double dxi = r*sin(thi)-xi+xm;
    double dyi = -r*cos(thi)-yi+ym;
    double dxf = xf-r*sin(thf)-xm;
    double dyf = r*cos(thf)+yf-ym;

    double t2 = dxf * dxf;
    double t3 = dyi * dyi;
    double t4 = t3 * t2;
    double t7 = 2 * r * dyi * t2;
    double t8 = r * r;
    double t9 = t2 * t8;
    double t10 = 3 * t9;
    double t12 = dxf * dxi;
    double t13 = t12 * dyf * dyi;
    double t14 = 2 * t13;
    double t15 = dyf * r;
    double t17 = 2 * t15 * t12;
    double t20 = 2 * dyi * r * t12;
    double t22 = dxi * dxf * t8;
    double t23 = 6 * t22;
    double t24 = dyf * dyf;
    double t25 = dxi * dxi;
    double t26 = t24 * t25;
    double t29 = 2 * r * dyf * t25;
    double t30 = t25 * t8;
    double t31 = 3 * t30;
    double t33 = 4 * t24 * t8;
    double t34 = dyf * t8;
    double t36 = 8 * dyi * t34;
    double t38 = 4 * t8 * t3;
    double t39 = dxf * dyi;
    double t40 = t39 * t15;
    double t42 = t8 * dyf * dxf;
    double t44 = r * t3 * dxf;
    double t45 = t8 * t39;
    double t47 = t24 * dxi * r;
    double t49 = dyi * dxi * t15;
    double t50 = dxi * t34;
    double t52 = dxi * dyi * t8;

    Eigen::Matrix<double, 5, 1> coefficients;
    coefficients << t4 - t7 - t10 - t14 + t17 - t20 - t23 + t26 + t29 - t31 - t33 - t36 - t38,
                    -4 * t40 + 4 * t42 - 4 * t44 + 4 * t45 + 4 * t47 + 4 * t49 + 4 * t50 + 4 * t52,
                    2 * t4 - 10 * t9 - 4 * t13 - 20 * t22 + 2 * t26 - 10 * t30 - t33 - t36 - t38,
                    -4 * t40 - 4 * t42 - 4 * t44 - 4 * t45 + 4 * t47 + 4 * t49 - 4 * t50 - 4 * t52,
                    t4 + t7 - t10 - t14 - t17 + t20 - t23 + t26 - t29 - t31 - t33 - t36 - t38;

    return coefficients;
}


Eigen::MatrixXd find_coefficients_man_14_p8 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
){
    double dxi = r*sin(thi)-xi+xm;
    double dyi = -r*cos(thi)-yi+ym;
    double dxf = xf-r*sin(thf)-xm;
    double dyf = r*cos(thf)+yf-ym;

    double t2 = dxf * dxf;
    double t3 = dyi * dyi;
    double t4 = t3 * t2;
    double t5 = dyi * t2;
    double t6 = r * t5;
    double t7 = 2 * t6;
    double t8 = r * r;
    double t9 = t2 * t8;
    double t10 = 3 * t9;
    double t12 = dxf * dxi;
    double t13 = t12 * dyf * dyi;
    double t14 = 2 * t13;
    double t15 = dyf * r;
    double t16 = t15 * t12;
    double t17 = 2 * t16;
    double t19 = dyi * r * t12;
    double t20 = 2 * t19;
    double t22 = dxi * dxf * t8;
    double t23 = 6 * t22;
    double t24 = dyf * dyf;
    double t25 = dxi * dxi;
    double t26 = t25 * t24;
    double t28 = r * dyf * t25;
    double t29 = 2 * t28;
    double t30 = t25 * t8;
    double t31 = 3 * t30;
    double t32 = t24 * t8;
    double t33 = 4 * t32;
    double t34 = dyf * t8;
    double t35 = dyi * t34;
    double t36 = 8 * t35;
    double t37 = t8 * t3;
    double t38 = 4 * t37;
    double t39 = dxi * t5;
    double t40 = 8 * t39;
    double t43 = 8 * dxi * t2 * r;
    double t44 = t25 * dxf;
    double t45 = dyf * t44;
    double t46 = 8 * t45;
    double t48 = 8 * r * t44;
    double t49 = dyf * dxf;
    double t50 = t3 * t49;
    double t51 = 8 * t50;
    double t52 = dxf * dyi;
    double t53 = t52 * t15;
    double t54 = 12 * t53;
    double t56 = 4 * t8 * t49;
    double t58 = r * t3 * dxf;
    double t59 = 4 * t58;
    double t61 = 28 * t8 * t52;
    double t63 = dxi * dyi * t24;
    double t64 = 8 * t63;
    double t66 = t24 * dxi * r;
    double t67 = 4 * t66;
    double t69 = dyi * dxi * t15;
    double t70 = 12 * t69;
    double t72 = 28 * dxi * t34;
    double t75 = 4 * dxi * dyi * t8;
    double t76 = t25 * t2;
    double t77 = 16 * t76;
    double t78 = 12 * t4;
    double t79 = 12 * t6;
    double t80 = 16 * t9;
    double t81 = 56 * t13;
    double t82 = 28 * t16;
    double t83 = 28 * t19;
    double t84 = 32 * t22;
    double t85 = 12 * t26;
    double t86 = 12 * t28;
    double t87 = 16 * t30;
    double t88 = t3 * t24;
    double t89 = 16 * t88;
    double t92 = 16 * dyi * t24 * r;
    double t93 = 12 * t32;
    double t95 = 16 * t3 * t15;
    double t96 = 24 * t35;
    double t97 = 12 * t37;
    double t98 = 56 * t39;
    double t99 = 56 * t45;
    double t100 = 56 * t50;
    double t101 = 28 * t53;
    double t102 = 20 * t58;
    double t103 = 56 * t63;
    double t104 = 20 * t66;
    double t105 = 28 * t69;

    Eigen::Matrix<double, 9, 1> coefficients;
    coefficients << t4 - t7 - t10 + t14 - t17 + t20 + t23 + t26 + t29 - t31 - t33 - t36 - t38,
                    -t40 + t43 - t46 - t48 + t51 - t54 + t56 + t59 + t61 + t64 - t67 + t70 + t72 + t75,
                    t77 - t78 + t79 - t80 - t81 + t82 - t83 - t84 - t85 - t86 - t87 + t89 - t92 - t93 + t95 + t96 - t97,
                    t98 - t43 + t99 + t48 - t100 + t101 + t56 - t102 + t61 - t103 + t104 - t105 + t72 + t75,
                    -32 * t76 + 38 * t4 - 26 * t9 + 140 * t13 - 76 * t22 + 38 * t26 - 26 * t30 - 32 * t88 - 16 * t32 + 64 * t35 - 16 * t37,
                    -t98 - t43 - t99 + t48 + t100 + t101 - t56 - t102 - t61 + t103 + t104 - t105 - t72 - t75,
                    t77 - t78 - t79 - t80 - t81 - t82 + t83 - t84 - t85 + t86 - t87 + t89 + t92 - t93 - t95 + t96 - t97,
                    t40 + t43 + t46 - t48 - t51 - t54 - t56 + t59 - t61 - t64 - t67 + t70 - t72 - t75,
                    t4 + t7 - t10 + t14 + t17 - t20 + t23 + t26 - t29 - t31 - t33 - t36 - t38;

    return coefficients;
}


std::vector<double> solve_man_14(
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance
){
    auto th_14_4 = solve_man<5>(
        find_coefficients_man_14_p4,
        xi, yi, thi,
        xm, ym,
        xf, yf, thf,
        r, imaginary_tolerance
    );

    auto th_14_8 = solve_man<9>(
        find_coefficients_man_14_p8,
        xi, yi, thi,
        xm, ym,
        xf, yf, thf,
        r, imaginary_tolerance
    );

    std::vector<double> result;
    result.reserve(th_14_4.size() + th_14_8.size());
    result.insert(result.end(), th_14_4.begin(), th_14_4.end());
    result.insert(result.end(), th_14_8.begin(), th_14_8.end());

    return result;
}