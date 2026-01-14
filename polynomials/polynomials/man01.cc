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

Eigen::MatrixXd find_coefficients_01 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
){
    double dxf = xf+r*sin(thf)-xm; 
    double dxi = -r*sin(thi)-xi+xm; 
    double dyf = -r*cos(thf)+yf-ym;
    double dyi = r*cos(thi)-yi+ym;

    double t2 = dxf * dxf;
    double t3 = t2 * t2;
    double t4 = dxi * dxi;
    double t5 = t4 * t3;
    double t6 = t4 * t4;
    double t7 = t6 * t2;
    double t8 = t4 * t2;
    double t9 = dyf * dyf;
    double t11 = 2 * t9 * t8;
    double t12 = dyf * r;
    double t13 = t12 * t8;
    double t14 = 4 * t13;
    double t15 = dyi * dyi;
    double t17 = 2 * t15 * t8;
    double t18 = dyi * r;
    double t19 = t18 * t8;
    double t20 = 4 * t19;
    double t21 = t15 * t15;
    double t22 = t21 * t2;
    double t23 = t15 * dyi;
    double t25 = r * t23 * t2;
    double t26 = 4 * t25;
    double t28 = r * r;
    double t29 = t28 * t15 * t2;
    double t30 = 10 * t29;
    double t32 = t28 * r;
    double t33 = t32 * dyi * t2;
    double t34 = 28 * t33;
    double t35 = t28 * t28;
    double t37 = 15 * t35 * t2;
    double t38 = t9 * t9;
    double t39 = t38 * t4;
    double t40 = t9 * dyf;
    double t42 = r * t40 * t4;
    double t43 = 4 * t42;
    double t44 = t9 * t4;
    double t45 = t28 * t44;
    double t46 = 10 * t45;
    double t48 = t32 * dyf * t4;
    double t49 = 28 * t48;
    double t51 = 15 * t35 * t4;
    double t54 = 4 * dyi * dxi * t3;
    double t55 = t2 * dxf;
    double t57 = r * t4 * t55;
    double t58 = 8 * t57;
    double t59 = t4 * dxi;
    double t61 = r * t59 * t2;
    double t62 = 8 * t61;
    double t63 = dxi * t2;
    double t64 = dyi * t9;
    double t66 = 8 * t64 * t63;
    double t68 = dyf * dyi * r;
    double t69 = t68 * t63;
    double t70 = 16 * t69;
    double t71 = r * t15;
    double t72 = t71 * t63;
    double t73 = 8 * t72;
    double t74 = t28 * dyi;
    double t76 = 40 * t74 * t63;
    double t77 = t32 * t63;
    double t78 = 56 * t77;
    double t81 = 4 * dyf * t6 * dxf;
    double t82 = t4 * dxf;
    double t84 = r * t9 * t82;
    double t85 = 8 * t84;
    double t86 = t15 * dyf;
    double t88 = 8 * t86 * t82;
    double t89 = t68 * t82;
    double t90 = 16 * t89;
    double t91 = -t54 + t58 + t62 - t66 - t70 + t73 + t76 - t78 + t81 + t85 + t88 - t90;
    double t94 = 40 * t28 * dyf * t82;
    double t95 = t32 * t82;
    double t96 = 56 * t95;
    double t97 = dyf * dxf;
    double t99 = 4 * t21 * t97;
    double t101 = r * t23 * t97;
    double t102 = 16 * t101;
    double t105 = 40 * t28 * t15 * t97;
    double t106 = t32 * dyi;
    double t107 = t106 * t97;
    double t108 = 112 * t107;
    double t110 = 60 * t35 * t97;
    double t113 = 4 * dyi * t38 * dxi;
    double t115 = t18 * t40 * dxi;
    double t116 = 16 * t115;
    double t117 = t9 * dxi;
    double t119 = 40 * t74 * t117;
    double t121 = t106 * dxi * dyf;
    double t122 = 112 * t121;
    double t125 = 60 * t35 * dyi * dxi;
    double t126 = -t94 - t96 + t99 - t102 - t105 + t108 - t110 - t113 - t116 + t119 + t122 + t125;
    double t128 = 4 * t15 * t3;
    double t131 = 32 * t18 * dxi * t55;
    double t132 = 12 * t13;
    double t133 = 12 * t19;
    double t136 = 8 * t15 * t9 * t2;
    double t139 = 16 * t71 * dyf * t2;
    double t140 = 12 * t25;
    double t141 = 50 * t29;
    double t142 = 84 * t33;
    double t145 = 32 * t12 * t59 * dxf;
    double t146 = dxf * dxi;
    double t149 = 32 * r * t64 * t146;
    double t152 = 32 * r * t86 * t146;
    double t153 = t5 - t128 + t131 - t7 + t11 + t132 - t17 + t133 - t136 - t139 - t22 + t140 + t141 - t142 + t37 + t145 + t149 + t152;
    double t156 = 224 * t32 * dyf * t146;
    double t158 = 224 * t106 * t146;
    double t160 = 4 * t9 * t6;
    double t161 = 12 * t42;
    double t163 = 8 * t15 * t44;
    double t165 = 16 * t18 * t44;
    double t166 = 50 * t45;
    double t167 = 84 * t48;
    double t169 = 4 * t15 * t38;
    double t172 = 16 * r * t15 * t40;
    double t174 = 4 * t21 * t9;
    double t177 = 16 * r * t23 * t9;
    double t179 = 112 * t32 * t64;
    double t181 = 60 * t35 * t9;
    double t183 = 112 * t32 * t86;
    double t185 = 60 * t35 * t15;
    double t186 = -t156 - t158 + t160 + t39 + t161 + t163 - t165 - t166 - t167 - t51 - t169 - t172 + t174 - t177 + t179 - t181 + t183 + t185;
    double t215 = t5 - t128 - t131 - t7 + t11 - t132 - t17 - t133 - t136 + t139 - t22 - t140 + t141 + t142 + t37 - t145 - t149 - t152;
    double t216 = t156 + t158 + t160 + t39 - t161 + t163 + t165 - t166 + t167 - t51 - t169 + t172 + t174 + t177 - t179 - t181 - t183 + t185;
    double t217 = t54 + t58 + t62 + t66 - t70 + t73 - t76 - t78 - t81 + t85 - t88 - t90;
    double t218 = t94 - t96 - t99 - t102 + t105 + t108 + t110 + t113 - t116 - t119 + t122 - t125;

    Eigen::Matrix<double, 7, 1> coefficients;

    coefficients << -t5 + t7 - t11 - t14 + t17 - t20 + t22 - t26 - t30 + t34 - t37 - t39 - t43 + t46 + t49 + t51,
                    t126 + t91,
                    t153 + t186,
                    -224 * dxf * t15 * t32 + 32 * dxf * t71 * t9 + 32 * r * t15 * t55 + 32 * r * t59 * t9 - 224 * t117 * t32 + 32 * t117 * t71 + 32 * t101 - 224 * t107 + 32 * t115 - 224 * t121 - 16 * t57 - 16 * t61 + 32 * t69 - 16 * t72 + 112 * t77 - 16 * t84 + 32 * t89 + 112 * t95,
                    t215 + t216,
                    t217 + t218,
                    -t5 + t7 - t11 + t14 + t17 + t20 + t22 + t26 - t30 - t34 - t37 - t39 + t43 + t46 - t49 + t51;

    return coefficients;
}


std::vector<double> solve_man_01(
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance
) {
    return solve_man<7>(
        find_coefficients_01,
        xi, yi, thi,
        xm, ym,
        xf, yf, thf,
        r, imaginary_tolerance
    );
}
