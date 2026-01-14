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


Eigen::MatrixXd find_coefficients_man_13_p4 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
){
    double dxf = xf-r*sin(thf)-xm; 
    double dxi = -r*sin(thi)-xi+xm; 
    double dyf = r*cos(thf)+yf-ym; 
    double dyi = r*cos(thi)-yi+ym;

    double t2 = dyi * dyi;
    double t3 = dxf * dxf;
    double t4 = t2 * t3;
    double t7 = 2 * r * dyi * t3;
    double t8 = r * r;
    double t9 = t8 * t3;
    double t10 = dxf * dxi;
    double t11 = dyf * dyi;
    double t12 = t11 * t10;
    double t13 = 2 * t12;
    double t16 = 2 * dyf * r * t10;
    double t17 = dyi * r;
    double t19 = 2 * t17 * t10;
    double t21 = dxi * dxf * t8;
    double t22 = 2 * t21;
    double t23 = dxi * dxi;
    double t24 = dyf * dyf;
    double t25 = t23 * t24;
    double t28 = 2 * r * dyf * t23;
    double t29 = t8 * t23;
    double t30 = 3 * t29;
    double t32 = 4 * t8 * t2;
    double t33 = r * t8;
    double t35 = 8 * t33 * dyi;
    double t36 = t8 * t8;
    double t37 = 4 * t36;
    double t38 = dxf * dyf;
    double t40 = 4 * t17 * t38;
    double t42 = 4 * t8 * t38;
    double t45 = 4 * r * t2 * dxf;
    double t48 = 4 * t8 * dyi * dxf;
    double t51 = 4 * r * t24 * dxi;
    double t52 = dxi * dyf;
    double t54 = 4 * t17 * t52;
    double t56 = 4 * t8 * t52;
    double t59 = 4 * t8 * dyi * dxi;
    double t61 = 16 * t33 * dxi;

    Eigen::MatrixXd coefficients(5,1);
    coefficients << t4 - t7 + t9 - t13 + t16 - t19 + t22 + t25 + t28 - t30 - t32 + t35 - t37,
                    -t40 + t42 - t45 + t48 + t51 + t54 + t56 + t59 - t61,
                    8 * t11 * t8 + 4 * t24 * t8 - 4 * t12 - 4 * t21 + 2 * t25 - 10 * t29 - t32 - 8 * t36 + 2 * t4 - 2 * t9,
                    -t40 - t42 - t45 - t48 + t51 + t54 - t56 - t59 - t61,
                    t4 + t7 + t9 - t13 - t16 + t19 + t22 + t25 - t28 - t30 - t32 - t35 - t37;

    return coefficients;
}


Eigen::MatrixXd find_coefficients_man_13_p8 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
){    
    double dxf = xf-r*sin(thf)-xm; 
    double dxi = -r*sin(thi)-xi+xm; 
    double dyf = r*cos(thf)+yf-ym; 
    double dyi = r*cos(thi)-yi+ym;

    double t2 = r * r;
    double t3 = t2 * t2;
    double t4 = 4 * t3;
    double t5 = r * t2;
    double t6 = t5 * dyi;
    double t7 = 8 * t6;
    double t8 = dxf * dxf;
    double t9 = dxf * dxi;
    double t11 = dxi * dxi;
    double t13 = dyi * dyi;
    double t16 = t2 * (t8 - 2 * t9 - 3 * t11 - 4 * t13);
    double t17 = dxf - dxi;
    double t18 = dxf * dyi;
    double t20 = dxi * dyf + t18;
    double t23 = 2 * r * t20 * t17;
    double t24 = t20 * t20;
    double t25 = t5 * dxi;
    double t26 = 16 * t25;
    double t27 = dyf - dyi;
    double t30 = 4 * t2 * t17 * t27;
    double t32 = 8 * t11 * dxf;
    double t33 = 8 * t8;
    double t34 = dyf * dyf;
    double t36 = dyf * dyi;
    double t45 = r * (-t32 + dxi * (t33 - 4 * t34 + 12 * t36) - 12 * (-dyi / 3 + dyf) * t18);
    double t47 = t20 * (t9 - t36);
    double t48 = 8 * t47;
    double t49 = 16 * t3;
    double t50 = 16 * t6;
    double t58 = t2 * (-16 * t11 + 4 * (dyf + dyi) * (dyf - 3 * dyi));
    double t60 = 12 * t11 * dyf;
    double t63 = 28 * dxf * t27 * dxi;
    double t65 = 16 * t13 * dyf;
    double t68 = 12 * t8 - 16 * t34;
    double t75 = t11 * (16 * t8 - 12 * t34);
    double t76 = t36 * t9;
    double t77 = 56 * t76;
    double t78 = -t68;
    double t79 = t13 * t78;
    double t80 = 48 * t25;
    double t91 = r * (t32 + dxi * (-t33 + 20 * t34 - 28 * t36) + 28 * (dyf - 0.5e1 / 0.7e1 * dyi) * dyi * dxf);
    double t92 = 56 * t47;

    Eigen::MatrixXd coefficients(9,1);
    coefficients << -t4 + t7 + t16 - t23 + t24,
                    -t26 + t30 + t45 - t48,
                    -t49 + t50 + t58 + r * (dyi * t68 - t60 + t63 + t65) + t75 - t77 + t79,
                    -t80 + t30 + t91 + t92,
                    -24 * t3 + t2 * (-2 * t8 + 4 * t9 - 26 * t11 + 8 * t34 - 16 * t36 - 16 * t13) + t11 * (-32 * t8 + 38 * t34) + 140 * t76 + t13 * (38 * t8 - 32 * t34),
                    -t80 - t30 + t91 - t92,
                    -t49 - t50 + t58 + r * (dyi * t78 + t60 - t63 - t65) + t75 - t77 + t79,
                    -t26 - t30 + t45 + t48,
                    -t4 - t7 + t16 + t23 + t24;

    return coefficients;
}


std::vector<double> solve_man_13(
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r, double imaginary_tolerance
){

    auto th_13_4 = solve_man<5>(
        find_coefficients_man_13_p4,
        xi, yi, thi, 
        xm, ym, 
        xf, yf, thf,
        r, 1e-6
    );

    auto th_13_8 = solve_man<9>(
        find_coefficients_man_13_p8,
        xi, yi, thi, 
        xm, ym, 
        xf, yf, thf,
        r, 1e-6
    );

    std::vector<double> result;
    result.reserve(th_13_4.size() + th_13_8.size());
    result.insert(result.end(), th_13_4.begin(), th_13_4.end());
    result.insert(result.end(), th_13_8.begin(), th_13_8.end());

    return result;
}