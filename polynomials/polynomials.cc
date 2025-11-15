#include "polynomials.hh"

Eigen::MatrixXd find_coefficients_man (
    size_t index,
    double xi, double yi, double thi,
    double xm, double ym,
    double xf, double yf, double thf,
    double r
){
    switch(index){
        case 1:
            return find_coefficients_man_01(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 2:
            return find_coefficients_man_02(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 3:
            return find_coefficients_man_03(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 4:
            return find_coefficients_man_04(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 5:
            return find_coefficients_man_05(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 6:
            return find_coefficients_man_06(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 7:
            return find_coefficients_man_07(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 8:
            return find_coefficients_man_08(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 9:
            return find_coefficients_man_09(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 10:
            return find_coefficients_man_10(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 11:
            return find_coefficients_man_11(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 12:
            return find_coefficients_man_12(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 13:
            return find_coefficients_man_13(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 14:
            return find_coefficients_man_14(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 15:
            return find_coefficients_man_15(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 16:
            return find_coefficients_man_16(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 17:
            return find_coefficients_man_17(xi, yi, thi, xm, ym, xf, yf, thf, r);
        case 18:
            return find_coefficients_man_18(xi, yi, thi, xm, ym, xf, yf, thf, r);
        default:
            throw std::invalid_argument("Invalid index for find_coefficients_man");
    }
}