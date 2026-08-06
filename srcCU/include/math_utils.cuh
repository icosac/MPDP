/**
 * @file math_utils.cuh
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Precision-parametric math helpers for the GPU solvers.
 */

#ifndef MPDP_MATH_UTILS_CUH
#define MPDP_MATH_UTILS_CUH

#include <cuda_runtime.h>

#include <cmath>

namespace mpdp {
namespace gpu {

/*!
 * Scalar traits: math functions and constants for one precision.
 */
template <typename T>
struct Scalar;

//! Double precision specialisation.
template <>
struct Scalar<double> {
	using type = double;

	__host__ __device__ static __forceinline__ double sin (double x) { return ::sin (x); }
	__host__ __device__ static __forceinline__ double cos (double x) { return ::cos (x); }
	__host__ __device__ static __forceinline__ double acos (double x) { return ::acos (x); }
	__host__ __device__ static __forceinline__ double sqrt (double x) { return ::sqrt (x); }
	__host__ __device__ static __forceinline__ double abs (double x) { return ::fabs (x); }
	__host__ __device__ static __forceinline__ double
	atan2 (double y, double x)
	{
		return ::atan2 (y, x);
	}
	__host__ __device__ static __forceinline__ double
	hypot (double x, double y)
	{
		return ::hypot (x, y);
	}

	__host__ __device__ static __forceinline__ double
	pi()
	{
		return 3.14159265358979323846264338328;
	}
	__host__ __device__ static __forceinline__ double
	twoPi()
	{
		return 6.28318530717958647692528676656;
	}
	//! A large finite value used as "no path yet"; stays finite under addition.
	__host__ __device__ static __forceinline__ double
	huge()
	{
		return 1.0e300;
	}
};

//! Single precision specialisation.
template <>
struct Scalar<float> {
	using type = float;

	__host__ __device__ static __forceinline__ float sin (float x) { return ::sinf (x); }
	__host__ __device__ static __forceinline__ float cos (float x) { return ::cosf (x); }
	__host__ __device__ static __forceinline__ float acos (float x) { return ::acosf (x); }
	__host__ __device__ static __forceinline__ float sqrt (float x) { return ::sqrtf (x); }
	__host__ __device__ static __forceinline__ float abs (float x) { return ::fabsf (x); }
	__host__ __device__ static __forceinline__ float
	atan2 (float y, float x)
	{
		return ::atan2f (y, x);
	}
	__host__ __device__ static __forceinline__ float
	hypot (float x, float y)
	{
		return ::hypotf (x, y);
	}

	__host__ __device__ static __forceinline__ float
	pi()
	{
		return 3.14159265358979323846264338328f;
	}
	__host__ __device__ static __forceinline__ float
	twoPi()
	{
		return 6.28318530717958647692528676656f;
	}
	__host__ __device__ static __forceinline__ float
	huge()
	{
		return 1.0e30f;
	}
};

/*!
 * Standardises an angle to [0, 2*pi).
 *
 * @param ang The angle to standardise.
 * @return The standardised angle.
 */
template <typename T>
__host__ __device__ __forceinline__ T
mod2pi (T ang)
{
	const T two_pi = Scalar<T>::twoPi();
	while (ang < T (0)) { ang += two_pi; }
	while (ang >= two_pi) { ang -= two_pi; }
	return ang;
}

}	 // namespace gpu
}	 // namespace mpdp

#endif	// MPDP_MATH_UTILS_CUH
