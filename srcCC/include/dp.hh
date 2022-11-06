#ifndef DP_HH
#define DP_HH

#ifndef CURVE
#define CURVE Dubins
#endif

#include <utils.hh>
#include <typedefs.hh>
#include <configuration.hh>
#include <clothoidG1.hh>
#include <dubins.hh>
#include <rs.hh>

#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <omp.h>


namespace DP {
  /*!
   * The wrapper to call to solve the dynamic programming problem for multi-point path finding.
   * @param points A vector of points the path should go through.
   * @param discr The number of sampling to consider for each point.
   * @param fixedAngles A vector stating which angles should not be changed.
   * @param params A list of refinements, that is the number of parameters to pass to the curve constructor.
   * @param nRefinements The number of times the algorithm should be called back in order to narrow the sampling intervals finding more precise final values. Default is 1.
   * @param saveAngles If set to `true`, then the points angles are changed to the best angles found at the end of the algorithm, otherwise, they are only returned. Default is true.
   * @return A pair where the first element is the computed length of the curve and the second one is a vector containing the best angles.
   */
  std::pair<LEN_T, std::vector<Angle> >
  solveDP(CURVE_TYPE curveType, std::vector<Configuration2>& points, const std::vector<bool>& fixedAngles,
          std::vector<real_type> params, int discr, int nRefinements, bool saveAngles=true);
  
} //namespace DP



#endif //DP_HH

















