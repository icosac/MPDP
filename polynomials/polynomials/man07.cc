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

Eigen::MatrixXd find_coefficients_07 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
){
    double dxi = xm - xi - r*sin(thi);
    double dyi = ym - yi + r*cos(thi);
    double dxf = xf + r*sin(thf) - xm;
    double dyf = yf - r*cos(thf) - ym;

    double t2 = dxf * dxf;
    double t3 = t2 * t2;
    double t4 = dxi * dxi;
    double t5 = t4 * t3;
    double t6 = t4 * t2;
    double t7 = dyf * dyf;
    double t8 = t7 * t6;
    double t9 = 2 * t8;
    double t11 = dyf * r * t6;
    double t12 = 4 * t11;
    double t13 = r * r;
    double t14 = t13 * t6;
    double t15 = 2 * t14;
    double t16 = dyi * dyi;
    double t18 = t13 * t16 * t2;
    double t19 = 16 * t18;
    double t21 = t13 * r;
    double t22 = t21 * dyi * t2;
    double t23 = 32 * t22;
    double t24 = t13 * t13;
    double t25 = t24 * t2;
    double t26 = 16 * t25;
    double t27 = t7 * t7;
    double t28 = t27 * t4;
    double t29 = t7 * dyf;
    double t31 = r * t29 * t4;
    double t32 = 4 * t31;
    double t34 = t13 * t7 * t4;
    double t35 = 10 * t34;
    double t37 = t21 * dyf * t4;
    double t38 = 28 * t37;
    double t39 = t24 * t4;
    double t40 = 15 * t39;
    double t43 = 4 * dyi * dxi * t3;
    double t44 = t2 * dxf;
    double t47 = 8 * r * t4 * t44;
    double t48 = dxi * t2;
    double t49 = dyi * t7;
    double t51 = 8 * t49 * t48;
    double t52 = dyf * dyi;
    double t55 = 16 * r * t52 * t48;
    double t56 = t13 * dyi;
    double t57 = t56 * t48;
    double t58 = 56 * t57;
    double t60 = 64 * t21 * t48;
    double t61 = t4 * dxf;
    double t64 = 8 * r * t7 * t61;
    double t66 = t13 * dyf * t61;
    double t67 = 48 * t66;
    double t69 = 56 * t21 * t61;
    double t70 = dxf * dyf;
    double t73 = 64 * t13 * t16 * t70;
    double t74 = t21 * dyi;
    double t76 = 128 * t74 * t70;
    double t78 = 64 * t24 * t70;
    double t81 = 4 * dyi * t27 * dxi;
    double t83 = dyi * r;
    double t85 = 16 * t83 * t29 * dxi;
    double t86 = t7 * dxi;
    double t87 = t56 * t86;
    double t88 = 40 * t87;
    double t91 = 112 * t52 * dxi * t21;
    double t94 = 60 * t24 * dyi * dxi;
    double t95 = t16 * t3;
    double t96 = 4 * t95;
    double t99 = 32 * t83 * dxi * t44;
    double t100 = 8 * t11;
    double t101 = 16 * t14;
    double t103 = t16 * t7 * t2;
    double t104 = 8 * t103;
    double t106 = r * t16;
    double t108 = 16 * t106 * dyf * t2;
    double t109 = 56 * t18;
    double t110 = 64 * t22;
    double t111 = dxf * dxi;
    double t114 = 32 * r * t49 * t111;
    double t116 = t13 * t52 * t111;
    double t117 = 64 * t116;
    double t120 = 256 * t21 * dyf * t111;
    double t121 = t96 - t99 - t100 + t101 + t104 + t108 - t109 + t110 - t114 - t117 + t120;
    double t123 = 224 * t74 * t111;
    double t124 = 8 * t31;
    double t125 = 48 * t34;
    double t126 = 56 * t37;
    double t127 = t16 * t27;
    double t128 = 4 * t127;
    double t131 = 16 * r * t16 * t29;
    double t133 = t13 * t16 * t7;
    double t134 = 24 * t133;
    double t136 = 128 * t21 * t49;
    double t137 = t24 * t7;
    double t138 = 64 * t137;
    double t141 = 112 * t21 * t16 * dyf;
    double t142 = t24 * t16;
    double t143 = 60 * t142;
    double t144 = t123 - t124 + t125 + t126 + t128 + t131 + t134 - t136 + t138 - t141 - t143;
    double t147 = 32 * r * t16 * t44;
    double t148 = 8 * t57;
    double t149 = 112 * t66;
    double t152 = 32 * t106 * t7 * dxf;
    double t155 = 224 * t21 * t16 * dxf;
    double t156 = 104 * t87;
    double t158 = 256 * t21 * t86;
    double t174 = t96 + t99 + t100 + t101 + t104 - t108 - t109 - t110 + t114 - t117 - t120;
    double t175 = -t123 + t124 + t125 - t126 + t128 - t131 + t134 + t136 + t138 + t141 - t143;

    Eigen::Matrix<double, 9, 1> coefficients;
    coefficients << t5 + t9 + t12 + t15 + t19 - t23 + t26 + t28 + t32 - t35 - t38 - t40,
                    t43 - t47 + t51 + t55 - t58 + t60 - t64 + t67 + t69 + t73 - t76 + t78 + t81 + t85 - t88 - t91 - t94,
                    t121 + t144,
                    t43 + t47 - t147 + t51 - t55 + t148 - t60 + t64 + t149 - t69 - t152 + t76 + t78 + t155 + t81 - t85 - t156 + t158 + t91 - t94,
                    -2 * t5 + 8 * t95 - 4 * t8 - 36 * t14 + 16 * t103 - 80 * t18 - 32 * t25 + 128 * t116 - 2 * t28 + 180 * t34 + 30 * t39 + 8 * t127 - 16 * t133 + 128 * t137 - 120 * t142,
                    -t43 + t47 - t147 - t51 - t55 - t148 - t60 + t64 - t149 - t69 - t152 + t76 - t78 + t155 - t81 - t85 + t156 + t158 + t91 + t94,
                    t174 + t175,
                    -t43 - t47 - t51 + t55 + t58 + t60 - t64 - t67 + t69 - t73 - t76 - t78 - t81 + t85 + t88 - t91 + t94,
                    t5 + t9 - t12 + t15 + t19 + t23 + t26 + t28 - t32 - t35 + t38 - t40;

    return coefficients;
}

std::vector<double> solve_man_07(
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance
) {
    return solve_man<9>(
        find_coefficients_07,
        xi, yi, thi,
        xm, ym,
        xf, yf, thf,
        r, imaginary_tolerance
    );
}
