/**
 * @file rs.cuh
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Device-side Reeds-Shepp solver - skeleton, not implemented yet.
 */

#ifndef MPDP_RS_CUH
#define MPDP_RS_CUH

#include <curve.cuh>
#include <math_utils.cuh>

namespace mpdp {
namespace gpu {

/*!
 * Reeds-Shepp solver - not implemented yet.
 */
template <typename T>
struct ReedsShepp {
	using S = Scalar<T>;

	static constexpr int kNumParams		 = 1;
	static constexpr bool kImplemented = false;
	static constexpr CurveKind kKind	 = CurveKind::REEDS_SHEPP;

	//! Mirrors `Dubins::Solution` so the DP can stay generic.
	struct Solution {
		T len		 = T (0);
		T s1		 = T (0);
		T s2		 = T (0);
		T s3		 = T (0);
		T k1		 = T (0);
		T k2		 = T (0);
		T k3		 = T (0);
		int type = 0;
	};

	//! Placeholder that never reports a usable path.
	__host__ __device__ static __forceinline__ Solution
	solve (T, T, T, T, T, T, const T*)
	{
		Solution out;
		out.len = S::huge();
		return out;
	}

	//! Placeholder that never reports a usable path.
	__host__ __device__ static __forceinline__ T
	length (T x0, T y0, T th0, T x1, T y1, T th1, const T* params)
	{
		return solve (x0, y0, th0, x1, y1, th1, params).len;
	}
};

}	 // namespace gpu
}	 // namespace mpdp

#endif	// MPDP_RS_CUH
