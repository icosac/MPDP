/**
 * @file 3PMD.hh
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License Agero 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Header file for the 3 Points Markov-Dubins problem.
 */

#ifndef INC_3PMD_HH
#define INC_3PMD_HH

#include <configuration.hh>
#include <dp.hh>
#include <dubins.hh>
#include <timeperf.hh>

#include <map>
#include <tuple>
#include <set>

constexpr std::map<std::string, std::tuple<int, Dubins::D_TYPE, Dubins::D_TYPE>> P3DP_DICT = {
 {"RLRRLR", {1, Dubins::D_TYPE::RLR, Dubins::D_TYPE::RLR}},
 {"LRLLRL", {2, Dubins::D_TYPE::LRL, Dubins::D_TYPE::LRL}},
 {"RLRRSR", {3, Dubins::D_TYPE::RLR, Dubins::D_TYPE::RSR}},
 {"RLRRSL", {4, Dubins::D_TYPE::RLR, Dubins::D_TYPE::RSL}},
 {"LRLLSL", {5, Dubins::D_TYPE::LRL, Dubins::D_TYPE::LSL}},
 {"LRLLSR", {6, Dubins::D_TYPE::LRL, Dubins::D_TYPE::LSR}},
 {"RSRRLR", {7, Dubins::D_TYPE::RSR, Dubins::D_TYPE::RLR}},
 {"LSRRLR", {8, Dubins::D_TYPE::LSR, Dubins::D_TYPE::RLR}},
 {"RSLLRL", {9, Dubins::D_TYPE::RSL, Dubins::D_TYPE::LRL}},
 {"LSLLRL", {10, Dubins::D_TYPE::LSL, Dubins::D_TYPE::LRL}},
 {"RSRRSR", {11, Dubins::D_TYPE::RSR, Dubins::D_TYPE::RSR}},
 {"LSRRSR", {12, Dubins::D_TYPE::LSR, Dubins::D_TYPE::RSR}},
 {"RSRRSL", {13, Dubins::D_TYPE::RSR, Dubins::D_TYPE::RSL}},
 {"LSRRSL", {14, Dubins::D_TYPE::LSR, Dubins::D_TYPE::RSL}},
 {"LSLLSL", {15, Dubins::D_TYPE::LSL, Dubins::D_TYPE::LSL}},
 {"RSLLSL", {16, Dubins::D_TYPE::RSL, Dubins::D_TYPE::LSL}},
 {"LSLLSR", {17, Dubins::D_TYPE::LSL, Dubins::D_TYPE::LSR}},
 {"RSLLSR", {18, Dubins::D_TYPE::RSL, Dubins::D_TYPE::LSR}}
};

int main3PMDBruteForce();
int main3PMDBruteForceWithPlot();
void main3PDP();

void generateDataset3PDPCircle(int argc, char** argv);
void generateDataset3PDPRect(int argc, char** argv);

#endif //INC_3PMD_HH
