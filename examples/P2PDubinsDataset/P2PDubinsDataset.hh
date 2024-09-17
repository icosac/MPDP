// System includes
#include <random>
#include <fstream>
#include <iomanip>
#include <iostream>

// Library includes
#include <dubins.hh>
#include <timeperf.hh>

/**
 * @brief Function to generate a random dataset of Dubins for testing.
 * @param fix_initial_pos Whether the initial position should be fixed to (0,0,0) or not
 * @return
 */
int genDSDubinsP2P(bool fix_initial_pos = false);

int testDubins();