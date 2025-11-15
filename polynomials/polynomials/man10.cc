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

Eigen::MatrixXd find_coefficients_man_10 (
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
){
    double dxi = xm - xi + r * sin(thi);
    double dyi = ym - yi - r * cos(thi);
    double dxf = xf - r * sin(thf) - xm;
    double dyf = yf + r * cos(thf) - ym;

    double t2 = dxf * dxf;
    double t3 = t2 * t2;
    double t4 = dxi * dxi;
    double t5 = t4 * t3;
    double t6 = t4 * t2;
    double t7 = dyf * dyf;
    double t8 = t7 * t6;
    double t9 = 2 * t8;
    double t11 = dyf * t4;
    double t12 = t11 * t2 * r;
    double t13 = 4 * t12;
    double t14 = r * r;
    double t15 = t14 * t6;
    double t16 = 2 * t15;
    double t17 = dyi * dyi;
    double t19 = t14 * t17 * t2;
    double t20 = 16 * t19;
    double t22 = t14 * r;
    double t23 = t22 * dyi * t2;
    double t24 = 32 * t23;
    double t25 = t14 * t14;
    double t26 = t25 * t2;
    double t27 = 16 * t26;
    double t28 = t7 * t7;
    double t29 = t28 * t4;
    double t31 = t7 * dyf;
    double t32 = t31 * t4 * r;
    double t33 = 4 * t32;
    double t35 = t14 * t7 * t4;
    double t36 = 10 * t35;
    double t37 = t22 * t11;
    double t38 = 28 * t37;
    double t39 = t25 * t4;
    double t40 = 15 * t39;
    double t43 = 4 * dyi * dxi * t3;
    double t44 = t2 * dxf;
    double t47 = 8 * r * t4 * t44;
    double t48 = dxi * t2;
    double t49 = dyi * t7;
    double t51 = 8 * t49 * t48;
    double t52 = dyf * dyi;
    double t55 = 16 * r * t52 * t48;
    double t56 = t14 * dyi;
    double t57 = t56 * t48;
    double t58 = 56 * t57;
    double t60 = 64 * t22 * t48;
    double t61 = t4 * dxf;
    double t64 = 8 * r * t7 * t61;
    double t66 = t14 * dyf * t61;
    double t67 = 48 * t66;
    double t69 = 56 * t22 * t61;
    double t70 = dxf * dyf;
    double t73 = 64 * t14 * t17 * t70;
    double t74 = t22 * dyi;
    double t76 = 128 * t74 * t70;
    double t78 = 64 * t25 * t70;
    double t81 = 4 * dyi * t28 * dxi;
    double t83 = dyi * r;
    double t85 = 16 * t83 * t31 * dxi;
    double t86 = t7 * dxi;
    double t87 = t56 * t86;
    double t88 = 40 * t87;
    double t91 = 112 * t74 * dxi * dyf;
    double t94 = 60 * t25 * dyi * dxi;
    double t95 = t17 * t3;
    double t96 = 4 * t95;
    double t99 = 32 * t83 * dxi * t44;
    double t100 = 8 * t12;
    double t101 = 16 * t15;
    double t103 = t17 * t7 * t2;
    double t104 = 8 * t103;
    double t106 = r * t17;
    double t108 = 16 * t106 * dyf * t2;
    double t109 = 56 * t19;
    double t110 = 64 * t23;
    double t111 = dxf * dxi;
    double t114 = 32 * r * t49 * t111;
    double t116 = t14 * t52 * t111;
    double t117 = 64 * t116;
    double t120 = 256 * t22 * dyf * t111;
    double t121 = t96 + t99 + t100 + t101 + t104 - t108 - t109 - t110 + t114 - t117 - t120;
    double t123 = 224 * t74 * t111;
    double t124 = 8 * t32;
    double t125 = 48 * t35;
    double t126 = 56 * t37;
    double t127 = t17 * t28;
    double t128 = 4 * t127;
    double t131 = 16 * r * t17 * t31;
    double t133 = t14 * t17 * t7;
    double t134 = 24 * t133;
    double t136 = 128 * t22 * t49;
    double t137 = t25 * t7;
    double t138 = 64 * t137;
    double t141 = 112 * t22 * t17 * dyf;
    double t142 = t25 * t17;
    double t143 = 60 * t142;
    double t144 = -t123 + t124 + t125 - t126 + t128 - t131 + t134 + t136 + t138 + t141 - t143;
    double t147 = 32 * r * t17 * t44;
    double t148 = 8 * t57;
    double t149 = 112 * t66;
    double t152 = 32 * t106 * t7 * dxf;
    double t155 = 224 * t22 * t17 * dxf;
    double t156 = 104 * t87;
    double t158 = 256 * t22 * t86;
    double t174 = t96 - t99 - t100 + t101 + t104 + t108 - t109 + t110 - t114 - t117 + t120;
    double t175 = t123 - t124 + t125 + t126 + t128 + t131 + t134 - t136 + t138 - t141 - t143;

    Eigen::Matrix<double, 9, 1> coefficients;
    coefficients << t5 + t9 - t13 + t16 + t20 + t24 + t27 + t29 - t33 - t36 + t38 - t40,
                    t43 + t47 + t51 - t55 - t58 - t60 + t64 + t67 - t69 + t73 + t76 + t78 + t81 - t85 - t88 + t91 - t94,
                    t121 + t144,
                    t43 - t47 + t147 + t51 + t55 + t148 + t60 - t64 + t149 + t69 + t152 - t76 + t78 - t155 + t81 + t85 - t156 - t158 - t91 - t94,
                    -2 * t5 + 8 * t95 - 4 * t8 - 36 * t15 + 16 * t103 - 80 * t19 - 32 * t26 + 128 * t116 - 2 * t29 + 180 * t35 + 30 * t39 + 8 * t127 - 16 * t133 + 128 * t137 - 120 * t142,
                    -t43 - t47 + t147 - t51 + t55 - t148 + t60 - t64 - t149 + t69 + t152 - t76 - t78 - t155 - t81 + t85 + t156 - t158 - t91 + t94,
                    t174 + t175,
                    -t43 + t47 - t51 - t55 + t58 - t60 + t64 - t67 - t69 - t73 + t76 - t78 - t81 - t85 + t88 + t91 + t94,
                    t5 + t9 + t13 + t16 + t20 - t24 + t27 + t29 + t33 - t36 - t38 - t40;
    return coefficients;
}

int main() {
    double xi = 1;
    double yi = 0;
    double r = 1/1.1;
    double thi = 2.0/3.0*M_PI;
    double thf = -5.0/9.0*M_PI;
    double xm = cos(5.0/9.0*M_PI);
    double ym = sin(5.0/9.0*M_PI);
    double xf = cos(2.0/3.0*M_PI);
    double yf = sin(2.0/3.0*M_PI);

    auto start = std::chrono::high_resolution_clock::now();

    Eigen::MatrixXd coefficients = find_coefficients_man_10(xi, yi, thi, xf, yf, thf, xm, ym, r);

    double coeff0 = coefficients(0);
    double coeff1 = coefficients(1);
    double coeff2 = coefficients(2);
    double coeff3 = coefficients(3);
    double coeff4 = coefficients(4);
    double coeff5 = coefficients(5);
    double coeff6 = coefficients(6);
    double coeff7 = coefficients(7);
    double coeff8 = coefficients(8);

    // double coeff0 = t5 + t9 - t13 + t16 + t20 + t24 + t27 + t29 - t33 - t36 + t38 - t40;
    // double coeff1 = t43 + t47 + t51 - t55 - t58 - t60 + t64 + t67 - t69 + t73 + t76 + t78 + t81 - t85 - t88 + t91 - t94;
    // double coeff2 = t121 + t144;
    // double coeff3 = t43 - t47 + t147 + t51 + t55 + t148 + t60 - t64 + t149 + t69 + t152 - t76 + t78 - t155 + t81 + t85 - t156 - t158 - t91 - t94;
    // double coeff4 = -2 * t5 + 8 * t95 - 4 * t8 - 36 * t15 + 16 * t103 - 80 * t19 - 32 * t26 + 128 * t116 - 2 * t29 + 180 * t35 + 30 * t39 + 8 * t127 - 16 * t133 + 128 * t137 - 120 * t142;
    // double coeff5 = -t43 - t47 + t147 - t51 + t55 - t148 + t60 - t64 - t149 + t69 + t152 - t76 - t78 - t155 - t81 + t85 + t156 - t158 - t91 + t94;
    // double coeff6 = t174 + t175;
    // double coeff7 = -t43 + t47 - t51 - t55 + t58 - t60 + t64 - t67 - t69 - t73 + t76 - t78 - t81 - t85 + t88 + t91 + t94;
    // double coeff8 = t5 + t9 + t13 + t16 + t20 - t24 + t27 + t29 + t33 - t36 - t38 - t40;

    // Eigen::Matrix<double, 9, 1> coefficients;
    // coefficients << coeff0,
    //                 coeff1,
    //                 coeff2,
    //                 coeff3,
    //                 coeff4,
    //                 coeff5,
    //                 coeff6,
    //                 coeff7,
    //                 coeff8;

    double coeff[8] = {coeff0, coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7};

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Computing coefficients took: " << elapsed.count() << " seconds" << std::endl;
    start = std::chrono::high_resolution_clock::now();

    if (std::abs(coeff8) < 1e-12) {
        std::cerr << "Leading coefficient is too small, cannot build companion matrix." << std::endl;
        return 1;
    }

    Eigen::Matrix<double, 8, 8> companion = Eigen::Matrix<double, 8, 8>::Zero();
    for (int i = 1; i < 8; ++i) {
        companion(i, i - 1) = 1.0;
    }
    for (int i = 0; i < 8; ++i) {
        companion(0, i) = -coefficients(7 - i) / coeff8;
        // companion(0, i) = -coeff[7 - i] / coeff8;
    }

    Eigen::EigenSolver<Eigen::Matrix<double, 8, 8>> eigen_solver(companion);
    Eigen::VectorXcd roots = eigen_solver.eigenvalues();

    end = std::chrono::high_resolution_clock::now();
    elapsed = end - start;
    std::cout << "Computation Time: " << elapsed.count() << " seconds" << std::endl;
    
    std::cout << "Coeff 0: " << coeff0 << std::endl;
    std::cout << "Coeff 1: " << coeff1 << std::endl;
    std::cout << "Coeff 2: " << coeff2 << std::endl;
    std::cout << "Coeff 3: " << coeff3 << std::endl;
    std::cout << "Coeff 4: " << coeff4 << std::endl;
    std::cout << "Coeff 5: " << coeff5 << std::endl;
    std::cout << "Coeff 6: " << coeff6 << std::endl;
    std::cout << "Coeff 7: " << coeff7 << std::endl;
    std::cout << "Coeff 8: " << coeff8 << std::endl;

    std::cout << std::endl << "Roots of the polynomial:" << std::endl;
    std::vector<double> real_roots;
    const double imag_tolerance = 1e-8;

    for (int i = 0; i < roots.size(); ++i) {
        const auto& root = roots(i);
        std::cout << "Root " << i << ": " << root << std::endl;
        if (std::abs(root.imag()) < imag_tolerance) {
            real_roots.push_back(root.real());
        }
    }

    if (!real_roots.empty()) {
        std::cout << std::endl << "Real roots (|Im| < " << imag_tolerance << "):" << std::endl;
        for (const double real_root : real_roots) {
            std::cout << real_root << " " << 2.0*atan(real_root) << std::endl;
        }
    } else {
        std::cout << std::endl << "No roots with negligible imaginary part found." << std::endl;
    }

    // -2.565758111, -0.9374882144, 0.6208885655, 3.083874499

    return 0;
}
