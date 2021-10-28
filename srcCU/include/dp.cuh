#ifndef DP_CUH
#define DP_CUH

#ifndef CURVE
#define CURVE Dubins
#endif

#include <utils.cuh>
#include <typedefs.hh>
#include <configuration.cuh>
#include <dubins.cuh>
#include <constants.cuh>

#include <iostream>
#include <set>
#include <cmath>
#include <vector>
#include <sstream>
#include <algorithm>

namespace DP{
  /*!
 * The wrapper to call to solve the dynamic programming problem for multi-point path finding.
 * @param points A vector of points the path should go through.
 * @param discr The number of sampling to consider for each point.
 * @param fixedAngles A vector stating which angles should not be changed.
 * @param nRefinements The number of times the algorithm should be called back in order to narrow the sampling intervals finding more precise final values.
 * @param params A list of additional parameters to pass to the curve constructor. Default is nullptr.
 * @return A vector containing the total length of the path and the best angles in the rest of the positions.
 */
  std::vector<real_type> solveDP(std::vector<Configuration2>& points, int discr, const std::vector<bool> fixedAngles, std::vector<real_type> params, short type=2, bool guessInitialAnglesVal=false, uint nIter=1, uint threads=128, Angle _fullAngle=m_2pi);
} //namespace DP


#endif //DP_CUH

















