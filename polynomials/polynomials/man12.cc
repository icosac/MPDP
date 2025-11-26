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

#if __has_include(<unsupported/Eigen/Polynomials>) || __has_include("unsupported/Eigen/Polynomials")
#  if __has_include(<unsupported/Eigen/Polynomials>)
#    include <unsupported/Eigen/Polynomials>
#  else
#    include "unsupported/Eigen/Polynomials"
#  endif
#else
#  error "Eigen's unsupported Polynomial module is required."
#endif


Eigen::MatrixXd find_coefficients_man_12_p4 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
){
    double dxi = xm-xi+r*sin(thi);
    double dyi = ym-yi-r*cos(thi);
    double dxf = xf+r*sin(thf)-xm;
    double dyf = yf-r*cos(thf)-ym;

    double t2 = dxf * dxf;
    double t3 = dyi * dyi;
    double t4 = t2 * t3;
    double t7 = 2 * r * dyi * t2;
    double t8 = r * r;
    double t9 = t2 * t8;
    double t10 = 3 * t9;
    double t11 = dxf * dxi;
    double t12 = dyf * dyi;
    double t13 = t12 * t11;
    double t14 = 2 * t13;
    double t17 = 2 * dyf * r * t11;
    double t18 = dyi * r;
    double t20 = 2 * t18 * t11;
    double t22 = dxi * dxf * t8;
    double t23 = 2 * t22;
    double t24 = dxi * dxi;
    double t25 = dyf * dyf;
    double t26 = t25 * t24;
    double t29 = 2 * dyf * r * t24;
    double t30 = t24 * t8;
    double t32 = 4 * t25 * t8;
    double t33 = t8 * r;
    double t35 = 8 * dyf * t33;
    double t36 = t8 * t8;
    double t37 = 4 * t36;
    double t38 = dxf * dyf;
    double t40 = 4 * t18 * t38;
    double t42 = 4 * t8 * t38;
    double t45 = 4 * r * t3 * dxf;
    double t48 = 4 * t8 * dyi * dxf;
    double t50 = 16 * dxf * t33;
    double t53 = 4 * r * t25 * dxi;
    double t54 = dxi * dyf;
    double t56 = 4 * t18 * t54;
    double t58 = 4 * t8 * t54;
    double t61 = 4 * t8 * dyi * dxi;
    
    Eigen::Matrix<double, 5, 1> coefficients;
    coefficients << t4 - t7 - t10 - t14 + t17 - t20 + t23 + t26 + t29 + t30 - t32 - t35 - t37,
                    -t40 + t42 - t45 + t48 + t50 + t53 + t56 + t58 + t61,
                    8 * t12 * t8 + 4 * t3 * t8 - 4 * t13 - 4 * t22 + 2 * t26 - 2 * t30 - t32 - 8 * t36 + 2 * t4 - 10 * t9,
                    -t40 - t42 - t45 - t48 + t50 + t53 + t56 - t58 - t61,
                    t4 + t7 - t10 - t14 - t17 + t20 + t23 + t26 - t29 + t30 - t32 + t35 - t37;
 
    return coefficients;
}


Eigen::MatrixXd find_coefficients_man_12_p8 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
){

    double dxi = xm-xi+r*sin(thi);
    double dyi = ym-yi-r*cos(thi);
    double dxf = xf+r*sin(thf)-xm;
    double dyf = yf-r*cos(thf)-ym;

    double t2 = dxf * dxf;
    double t3 = dyi * dyi;
    double t4 = t3 * t2;
    double t6 = r * dyi * t2;
    double t7 = 2 * t6;
    double t8 = r * r;
    double t9 = t2 * t8;
    double t10 = 3 * t9;
    double t11 = dxf * dxi;
    double t12 = dyf * dyi;
    double t13 = t12 * t11;
    double t14 = 2 * t13;
    double t16 = dyf * r * t11;
    double t17 = 2 * t16;
    double t18 = dyi * r;
    double t19 = t18 * t11;
    double t20 = 2 * t19;
    double t22 = dxi * dxf * t8;
    double t23 = 2 * t22;
    double t24 = dxi * dxi;
    double t25 = dyf * dyf;
    double t26 = t25 * t24;
    double t28 = r * dyf * t24;
    double t29 = 2 * t28;
    double t30 = t24 * t8;
    double t31 = t25 * t8;
    double t32 = 4 * t31;
    double t33 = t8 * r;
    double t34 = dyf * t33;
    double t35 = 8 * t34;
    double t36 = t8 * t8;
    double t37 = 4 * t36;
    double t38 = dxi * t2;
    double t39 = dyi * t38;
    double t40 = 8 * t39;
    double t42 = 8 * r * t38;
    double t43 = t24 * dxf;
    double t44 = dyf * t43;
    double t45 = 8 * t44;
    double t47 = 8 * r * t43;
    double t48 = dxf * dyf;
    double t49 = t3 * t48;
    double t50 = 8 * t49;
    double t51 = t18 * t48;
    double t52 = 12 * t51;
    double t54 = 4 * t8 * t48;
    double t56 = r * t3 * dxf;
    double t57 = 4 * t56;
    double t60 = 4 * t8 * dyi * dxf;
    double t61 = dxf * t33;
    double t62 = 16 * t61;
    double t63 = t25 * dxi;
    double t64 = dyi * t63;
    double t65 = 8 * t64;
    double t66 = r * t63;
    double t67 = 4 * t66;
    double t68 = dxi * dyf;
    double t69 = t18 * t68;
    double t70 = 12 * t69;
    double t72 = 4 * t8 * t68;
    double t75 = 4 * t8 * dyi * dxi;
    double t76 = t24 * t2;
    double t77 = 16 * t76;
    double t78 = 12 * t4;
    double t79 = 12 * t6;
    double t80 = 16 * t9;
    double t81 = 56 * t13;
    double t82 = 28 * t16;
    double t83 = 28 * t19;
    double t84 = 12 * t26;
    double t85 = 12 * t28;
    double t86 = t25 * t3;
    double t87 = 16 * t86;
    double t90 = 16 * r * dyi * t25;
    double t91 = 12 * t31;
    double t94 = 16 * r * t3 * dyf;
    double t95 = t8 * t12;
    double t96 = 8 * t95;
    double t97 = 16 * t34;
    double t98 = t8 * t3;
    double t99 = 4 * t98;
    double t100 = 16 * t36;
    double t101 = 56 * t39;
    double t102 = 56 * t44;
    double t103 = 56 * t49;
    double t104 = 28 * t51;
    double t105 = 20 * t56;
    double t106 = 48 * t61;
    double t107 = 56 * t64;
    double t108 = 20 * t66;
    double t109 = 28 * t69;

    Eigen::Matrix<double, 9, 1> coefficients;
    coefficients << t4 - t7 - t10 + t14 - t17 + t20 - t23 + t26 + t29 + t30 - t32 - t35 - t37,
                    -t40 + t42 - t45 - t47 + t50 - t52 + t54 + t57 - t60 + t62 + t65 - t67 + t70 - t72 + t75,
                    t77 - t78 + t79 - t80 - t81 + t82 - t83 - t84 - t85 + t87 - t90 - t91 + t94 - t96 - t97 + t99 - t100,
                    t101 - t42 + t102 + t47 - t103 + t104 + t54 - t105 - t60 + t106 - t107 + t108 - t109 - t72 + t75,
                    -32 * t76 + 38 * t4 - 26 * t9 + 140 * t13 + 4 * t22 + 38 * t26 - 2 * t30 - 32 * t86 - 16 * t31 - 16 * t95 + 8 * t98 - 24 * t36,
                    -t101 - t42 - t102 + t47 + t103 + t104 - t54 - t105 + t60 + t106 + t107 + t108 - t109 + t72 - t75,
                    t77 - t78 - t79 - t80 - t81 - t82 + t83 - t84 + t85 + t87 + t90 - t91 - t94 - t96 + t97 + t99 - t100,
                    t40 + t42 + t45 - t47 - t50 - t52 - t54 + t57 + t60 + t62 - t65 - t67 + t70 + t72 - t75,
                    t4 + t7 - t10 + t14 + t17 - t20 - t23 + t26 - t29 + t30 - t32 + t35 - t37;

    return coefficients;
}


std::vector<double> solve_man_12(
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance
){

    auto th_12_4 = solve_man<5>(
        find_coefficients_man_12_p4,
        xi, yi, thi, 
        xm, ym, 
        xf, yf, thf,
        r, imaginary_tolerance
    );

    auto th_12_8 = solve_man<9>(
        find_coefficients_man_12_p8,
        xi, yi, thi, 
        xm, ym, 
        xf, yf, thf,
        r, imaginary_tolerance
    );

    std::vector<double> result;
    result.insert(result.end(), th_12_4.begin(), th_12_4.end());
    result.insert(result.end(), th_12_8.begin(), th_12_8.end());

    return result;
}