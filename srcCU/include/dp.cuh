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
 * @param params A list of parameters to pass to the curve constructor.
 * @param nRef The number of refinements, that is the number of time the algorithm should be called back in order to narrow the sampling intervals finding more precise final values. Default is 1.
 * @param saveAngles  If set to `true`, then the points angles are changed to the best angles found at the end of the algorithm, otherwise, they are only returned. Default is true.
 * @param type The type of the function. Default is 2.
 *  - 1: A GPU accelerated version of the CPU provided solution (best with smaller GPU memory);
 *  - 2: A GPU accelerated version of the CPU provided solution, plus it uses a pipeline to improve performances (best with larger GPU memory).
 * @param threads The number of threads to use to parallelize. The blocks are computed considering the number of samples. Default set to 128.
 * @param _fullAngle The initial period in which to consider the samples. Default is 2*\pi.
 * @return A vector containing the total length of the path and the best angles in the rest of the positions.
 */
  std::vector<real_type> solveDP(std::vector<Configuration2>& points, int discr, const std::vector<bool> fixedAngles,
                                 std::vector<real_type> params, uint nRefs=1, bool saveAngles=false,
                                 short type=2, uint threads=128, Angle _fullAngle=m_2pi);
} //namespace DP


#endif //DP_CUH

















