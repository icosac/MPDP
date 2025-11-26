#include <polynomials.hh>

#include <stdexcept>
#include <string>

namespace {
[[noreturn]] void throw_not_implemented(const char* label) {
    throw std::runtime_error(std::string("Polynomial case ") + label + " is not implemented yet.");
}
}  // namespace

Eigen::MatrixXd find_coefficients_02(
    double, double, double,
    double, double,
    double, double, double,
    double) {
    throw_not_implemented("02");
}

Eigen::MatrixXd find_coefficients_03(
    double, double, double,
    double, double,
    double, double, double,
    double) {
    throw_not_implemented("03");
}

Eigen::MatrixXd find_coefficients_04(
    double, double, double,
    double, double,
    double, double, double,
    double) {
    throw_not_implemented("04");
}

Eigen::MatrixXd find_coefficients_05(
    double, double, double,
    double, double,
    double, double, double,
    double) {
    throw_not_implemented("05");
}

Eigen::MatrixXd find_coefficients_06(
    double, double, double,
    double, double,
    double, double, double,
    double) {
    throw_not_implemented("06");
}

Eigen::MatrixXd find_coefficients_08(
    double, double, double,
    double, double,
    double, double, double,
    double) {
    throw_not_implemented("08");
}

Eigen::MatrixXd find_coefficients_09(
    double, double, double,
    double, double,
    double, double, double,
    double) {
    throw_not_implemented("09");
}

Eigen::MatrixXd find_coefficients_14(
    double, double, double,
    double, double,
    double, double, double,
    double) {
    throw_not_implemented("14");
}

Eigen::MatrixXd find_coefficients_15(
    double, double, double,
    double, double,
    double, double, double,
    double) {
    throw_not_implemented("15");
}

Eigen::MatrixXd find_coefficients_16(
    double, double, double,
    double, double,
    double, double, double,
    double) {
    throw_not_implemented("16");
}

Eigen::MatrixXd find_coefficients_17(
    double, double, double,
    double, double,
    double, double, double,
    double) {
    throw_not_implemented("17");
}

Eigen::MatrixXd find_coefficients_18(
    double, double, double,
    double, double,
    double, double, double,
    double) {
    throw_not_implemented("18");
}
