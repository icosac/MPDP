#ifndef DP_HH
#define DP_HH

#include <settings.hh>
#include <utils.hh>
#include <typedefs.hh>
#include <configuration.hh>
#include <clothoidG1.hh>
#include <dubins.hh>

#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <omp.h>

#ifndef CURVE
#define CURVE Dubins
#endif

namespace DP {
  /*!
   * The wrapper to call to solve the dynamic programming problem for multi-point path finding.
   * @param points A vector of points the path should go through.
   * @param discr The number of sampling to consider for each point.
   * @param fixedAngles A vector stating which angles should not be changed.
   * @param nRefinements The number of times the algorithm should be called back in order to narrow the sampling intervals finding more precise final values.
   * @param params A list of additional parameters to pass to the curve constructor. Default is nullptr.
   * @return A vector containing the total length of the path and the best angles in the rest of the positions.
   */
  std::vector<real_type> solveDP (std::vector<Configuration2> points, int discr, const std::vector<bool>& fixedAngles, int nRefinements, real_type* params=nullptr);
  
} //namespace DP



#endif //DP_HH

















