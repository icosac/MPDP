/**
 * @file dubins.cuh
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Device-side Dubins solver used by the GPU dynamic programming.
 */

#ifndef MPDP_DUBINS_CUH
#define MPDP_DUBINS_CUH

#include <curve.cuh>
#include <math_utils.cuh>

namespace mpdp {
namespace gpu {

template <typename T>
struct Dubins {
	using S = Scalar<T>;

	//! Number of entries `params` must hold: the maximum curvature.
	static constexpr int kNumParams = 1;
	//! Whether this solver is implemented.
	static constexpr bool kImplemented = true;
	//! The family this solver belongs to.
	static constexpr CurveKind kKind = CurveKind::DUBINS;

	//! The six Dubins words, plus the "no solution" marker.
	enum Type { INVALID = 0, LRL, RLR, LSL, LSR, RSL, RSR };

	//! Full solution of one point-to-point problem.
	struct Solution {
		T len	 = T (0);	 ///< Total length.
		T s1	 = T (0);	 ///< Length of the first arc.
		T s2	 = T (0);	 ///< Length of the middle arc.
		T s3	 = T (0);	 ///< Length of the last arc.
		T k1	 = T (0);	 ///< Curvature of the first arc.
		T k2	 = T (0);	 ///< Curvature of the middle arc.
		T k3	 = T (0);	 ///< Curvature of the last arc.
		int type = INVALID;	 ///< The chosen word.
	};

	/*!
	 * Solves the point-to-point Dubins problem.
	 * @param x0 Initial abscissa.
	 * @param y0 Initial ordinate.
	 * @param th0 Initial heading.
	 * @param x1 Final abscissa.
	 * @param y1 Final ordinate.
	 * @param th1 Final heading.
	 * @param params `params[0]` is the maximum curvature.
	 * @return The shortest of the six maneuvers.
	 */
	__host__ __device__ static __forceinline__ Solution
	solve (T x0, T y0, T th0, T x1, T y1, T th1, const T* params)
	{
		// ---- scale to standard ------------------------------------------------
		const T dx		 = x1 - x0;
		const T dy		 = y1 - y0;
		const T phi		 = S::atan2 (dy, dx);
		const T lambda = S::hypot (dx, dy) * T (0.5);
		const T sKmax	 = params[0] * lambda;
		const T sth0	 = mod2pi (th0 - phi);
		const T sth1	 = mod2pi (th1 - phi);

		const T invK	= T (1) / sKmax;
		const T sin_0 = S::sin (sth0);
		const T cos_0 = S::cos (sth0);
		const T sin_1 = S::sin (sth1);
		const T cos_1 = S::cos (sth1);

		const T Ksq	 = sKmax * sKmax;
		const T dcos = S::cos (sth0 - sth1);
		const T dcos2 = cos_0 - cos_1;
		const T dsin = sin_0 - sin_1;
		const T scos = cos_0 + cos_1;
		const T ssin = sin_0 + sin_1;
		const T dth	 = sth0 - sth1;

		const T two_pi = S::twoPi();

		T len = S::huge();
		T ss1 = T (0), ss2 = T (0), ss3 = T (0);
		T sk1 = T (0), sk2 = T (0), sk3 = T (0);
		int type = INVALID;

		T C, Sc, temp1, temp2, temp3, t1, t2, t3, lc;

		// ---- LRL --------------------------------------------------------------
		C		 = -dcos2;
		Sc	 = T (2) * sKmax + dsin;
		temp1 = S::atan2 (C, Sc);
		temp2 = T (0.125) * (T (6) - T (4) * Ksq + T (2) * dcos - T (4) * sKmax * dsin);
		if (S::abs (temp2) <= T (1))
		{
			t2 = invK * mod2pi (two_pi - S::acos (temp2));
			t1 = invK * mod2pi (-sth0 + temp1 + T (0.5) * t2 * sKmax);
			t3 = invK * mod2pi (-dth + (t2 - t1) * sKmax);
			lc = t1 + t2 + t3;
			if (lc < len)
			{
				len = lc;
				ss1 = t1; ss2 = t2; ss3 = t3;
				sk1 = T (1); sk2 = T (-1); sk3 = T (1);
				type = LRL;
			}
		}

		// ---- RLR --------------------------------------------------------------
		C		 = dcos2;
		Sc	 = T (2) * sKmax - dsin;
		temp1 = S::atan2 (C, Sc);
		temp2 = T (0.125) * (T (6) - T (4) * Ksq + T (2) * dcos + T (4) * sKmax * dsin);
		if (S::abs (temp2) <= T (1))
		{
			t2 = invK * mod2pi (two_pi - S::acos (temp2));
			t1 = invK * mod2pi (sth0 - temp1 + T (0.5) * t2 * sKmax);
			t3 = invK * mod2pi (dth + (t2 - t1) * sKmax);
			lc = t1 + t2 + t3;
			if (lc < len)
			{
				len = lc;
				ss1 = t1; ss2 = t2; ss3 = t3;
				sk1 = T (-1); sk2 = T (1); sk3 = T (-1);
				type = RLR;
			}
		}

		// ---- LSL --------------------------------------------------------------
		C		 = cos_1 - cos_0;
		Sc	 = T (2) * sKmax + dsin;
		temp1 = S::atan2 (C, Sc);
		temp2 = T (2) + T (4) * Ksq - T (2) * dcos + T (4) * sKmax * dsin;
		if (temp2 >= T (0))
		{
			temp3 = invK * S::sqrt (temp2);
			t1		= invK * mod2pi (temp1 - sth0);
			t2		= temp3;
			t3		= invK * mod2pi (sth1 - temp1);
			lc		= t1 + t2 + t3;
			if (lc < len)
			{
				len = lc;
				ss1 = t1; ss2 = t2; ss3 = t3;
				sk1 = T (1); sk2 = T (0); sk3 = T (1);
				type = LSL;
			}
		}

		// ---- LSR --------------------------------------------------------------
		C		 = scos;
		Sc	 = T (2) * sKmax + ssin;
		temp1 = S::atan2 (-C, Sc);
		temp2 = T (-2) + T (4) * Ksq + T (2) * dcos + T (4) * sKmax * ssin;
		if (temp2 >= T (0))
		{
			t2		= invK * S::sqrt (temp2);
			temp3 = -S::atan2 (T (-2), t2 * sKmax);
			t1		= invK * mod2pi (-sth0 + temp1 + temp3);
			t3		= invK * mod2pi (-sth1 + temp1 + temp3);
			lc		= t1 + t2 + t3;
			if (lc < len)
			{
				len = lc;
				ss1 = t1; ss2 = t2; ss3 = t3;
				sk1 = T (1); sk2 = T (0); sk3 = T (-1);
				type = LSR;
			}
		}

		// ---- RSL --------------------------------------------------------------
		C		 = scos;
		Sc	 = T (2) * sKmax - ssin;
		temp1 = S::atan2 (C, Sc);
		temp2 = T (-2) + T (4) * Ksq + T (2) * dcos - T (4) * sKmax * ssin;
		if (temp2 >= T (0))
		{
			t2		= invK * S::sqrt (temp2);
			temp3 = S::atan2 (T (2), t2 * sKmax);
			t1		= invK * mod2pi (sth0 - temp1 + temp3);
			t3		= invK * mod2pi (sth1 - temp1 + temp3);
			lc		= t1 + t2 + t3;
			if (lc < len)
			{
				len = lc;
				ss1 = t1; ss2 = t2; ss3 = t3;
				sk1 = T (-1); sk2 = T (0); sk3 = T (1);
				type = RSL;
			}
		}

		// ---- RSR --------------------------------------------------------------
		C		 = cos_0 - cos_1;
		Sc	 = T (2) * sKmax - dsin;
		temp1 = S::atan2 (C, Sc);
		temp2 = T (2) + T (4) * Ksq - T (2) * dcos - T (4) * sKmax * dsin;
		if (temp2 >= T (0))
		{
			temp3 = invK * S::sqrt (temp2);
			t1		= invK * mod2pi (sth0 - temp1);
			t2		= temp3;
			t3		= invK * mod2pi (temp1 - sth1);
			lc		= t1 + t2 + t3;
			if (lc < len)
			{
				len = lc;
				ss1 = t1; ss2 = t2; ss3 = t3;
				sk1 = T (-1); sk2 = T (0); sk3 = T (-1);
				type = RSR;
			}
		}

		// ---- scale back -------------------------------------------------------
		Solution out;
		out.s1	 = ss1 * lambda;
		out.s2	 = ss2 * lambda;
		out.s3	 = ss3 * lambda;
		out.k1	 = sk1 * params[0];
		out.k2	 = sk2 * params[0];
		out.k3	 = sk3 * params[0];
		out.len	 = out.s1 + out.s2 + out.s3;
		out.type = type;
		return out;
	}

	/*!
	 * Length-only entry point, which is all DP on GPU needs.
	 * @param x0 Initial abscissa.
	 * @param y0 Initial ordinate.
	 * @param th0 Initial heading.
	 * @param x1 Final abscissa.
	 * @param y1 Final ordinate.
	 * @param th1 Final heading.
	 * @param params `params[0]` is the maximum curvature.
	 * @return The length of the shortest maneuver.
	 */
	__host__ __device__ static __forceinline__ T
	length (T x0, T y0, T th0, T x1, T y1, T th1, const T* params)
	{
		return solve (x0, y0, th0, x1, y1, th1, params).len;
	}
};

}	 // namespace gpu
}	 // namespace mpdp

#endif	// MPDP_DUBINS_CUH
