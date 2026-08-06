/**
 * @file curve.cuh
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief What the GPU dynamic programming solver expects of a curve family.
 */

#ifndef MPDP_CURVE_CUH
#define MPDP_CURVE_CUH

namespace mpdp {
namespace gpu {

//! The point-to-point curve families the GPU DP can be instantiated with.
enum class CurveKind {
	DUBINS = 0,			///< Markov-Dubins path, forward motion only.
	REEDS_SHEPP = 1		///< Reeds-Shepp path, forward and backward motion.
};

}	 // namespace gpu
}	 // namespace mpdp

#endif	// MPDP_CURVE_CUH
