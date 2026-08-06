/**
 * @file rs.cuh
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Device-side Reeds-Shepp solver used by the GPU dynamic programming.
 *
 * A port of `RS::reeds_shepp` from srcCC/rs.cc: twelve base words, each tried
 * with four sign variants of the standardised problem, giving the classic 48
 * Reeds-Shepp maneuvers. The dispatch order and the maneuver numbering are
 * identical to the CPU, so a length difference can never be blamed on a
 * different search order.
 */

#ifndef MPDP_RS_CUH
#define MPDP_RS_CUH

#include <curve.cuh>
#include <math_utils.cuh>

namespace mpdp {
namespace gpu {

/*!
 * Reeds-Shepp path solver.
 */
template <typename T>
struct ReedsShepp {
	using S = Scalar<T>;

	//! Number of entries `params` must hold: the maximum curvature.
	static constexpr int kNumParams = 1;
	//! Whether this solver is implemented.
	static constexpr bool kImplemented = true;
	//! The family this solver belongs to.
	static constexpr CurveKind kKind = CurveKind::REEDS_SHEPP;

	/*!
	 * The chosen maneuver.
	 *
	 * `t`, `u` and `v` are the three maneuver parameters already multiplied by
	 * the turning radius, matching `RS::_L1`, `_L2` and `_L3` on the CPU. Note
	 * that for the words containing a straight segment `u` is a length rather
	 * than an angle, and srcCC scales it by the radius anyway; this port keeps
	 * that behaviour so the two agree field by field.
	 */
	struct Solution {
		T len		 = T (0);	 ///< Total length of the path.
		T t			 = T (0);	 ///< First maneuver parameter, times the turning radius.
		T u			 = T (0);	 ///< Second maneuver parameter, times the turning radius.
		T v			 = T (0);	 ///< Third maneuver parameter, times the turning radius.
		int man	 = 0;			 ///< Maneuver number, 1 to 48; 0 means none was found.
	};

 private:
	//! Tolerance below which two coordinates count as coincident, as in srcCC.
	__host__ __device__ static __forceinline__ T eps3() { return T (1.0e-14); }
	//! Pi, matching the literal used by srcCC/rs.cc.
	__host__ __device__ static __forceinline__ T mpi() { return T (3.1415926535897932385); }
	//! Pi/2, matching the literal used by srcCC/rs.cc.
	__host__ __device__ static __forceinline__ T mpiDiv2() { return T (1.5707963267948966192); }
	//! 2*pi, matching the literal used by srcCC/rs.cc.
	__host__ __device__ static __forceinline__ T mpiMul2() { return T (6.2831853071795864770); }

	/*!
	 * The atan2 variant srcCC/rs.cc uses, which returns an angle in [0, 2*pi)
	 * and is built on `atan(y/x)` rather than on the library `atan2`. It is not
	 * interchangeable with `atan2`: for y == 0 and x > 0 it returns 2*pi, not 0.
	 * @param y Ordinate.
	 * @param x Abscissa.
	 * @return The angle.
	 */
	__host__ __device__ static __forceinline__ T
	myAtan2 (T y, T x)
	{
		if (x == T (0) && y == T (0)) { return T (0); }
		if (x == T (0)) { return y > T (0) ? mpiDiv2() : -mpiDiv2(); }
		const T a = S::atan (y / x);
		if (a > T (0)) { return x > T (0) ? a : a + mpi(); }
		return x > T (0) ? a + mpiMul2() : a + mpi();
	}

	//! The best maneuver found so far.
	struct Best {
		T len;
		T t, u, v;
		int man;
	};

	/*!
	 * Keeps a candidate maneuver when it is shorter than the incumbent.
	 * @param b The incumbent.
	 * @param len Length of the candidate.
	 * @param t First parameter of the candidate.
	 * @param u Second parameter of the candidate.
	 * @param v Third parameter of the candidate.
	 * @param man Maneuver number of the candidate.
	 */
	__host__ __device__ static __forceinline__ void
	keep (Best& b, T len, T t, T u, T v, int man)
	{
		if (len < b.len)
		{
			b.len = len;
			b.t		= t;
			b.u		= u;
			b.v		= v;
			b.man = man;
		}
	}

	// ---- the twelve base words ------------------------------------------------
	// Each returns the length of the maneuver, or `huge()` when the geometry
	// makes it impossible, and writes its three parameters into t, u and v.

	//! C | C | C
	__host__ __device__ static __forceinline__ T
	c_c_c (T radcurv, T x, T y, T phi, T rs, T rc, T& t, T& u, T& v)
	{
		const T r4 = T (4) * radcurv;
		const T a = x - rs;
		const T b = y + rc;
		if (S::abs (a) < eps3() && S::abs (b) < eps3()) { return S::huge(); }
		const T u1 = S::sqrt (a * a + b * b);
		if (u1 > r4) { return S::huge(); }
		const T theta = myAtan2 (b, a);
		const T alpha = S::acos (u1 / r4);
		t = mod2pi (mpiDiv2() + alpha + theta);
		u = mod2pi (mpi() - T (2) * alpha);
		v = mod2pi (phi - t - u);
		return radcurv * (t + u + v);
	}

	//! C | C C
	__host__ __device__ static __forceinline__ T
	c_cc (T radcurv, T x, T y, T phi, T rs, T rc, T& t, T& u, T& v)
	{
		const T r4 = T (4) * radcurv;
		const T a = x - rs;
		const T b = y + rc;
		if (S::abs (a) < eps3() && S::abs (b) < eps3()) { return S::huge(); }
		const T u1 = S::sqrt (a * a + b * b);
		if (u1 > r4) { return S::huge(); }
		const T theta = myAtan2 (b, a);
		const T alpha = S::acos (u1 / r4);
		t = mod2pi (mpiDiv2() + alpha + theta);
		u = mod2pi (mpi() - T (2) * alpha);
		v = mod2pi (t + u - phi);
		return radcurv * (t + u + v);
	}

	//! C S C, first form
	__host__ __device__ static __forceinline__ T
	csca (T radcurv, T x, T y, T phi, T rs, T rc, T& t, T& u, T& v)
	{
		const T a = x - rs;
		const T b = y + rc;
		t = mod2pi (myAtan2 (b, a));
		u = S::sqrt (a * a + b * b);
		v = mod2pi (phi - t);
		return radcurv * (t + v) + u;
	}

	//! C S C, second form
	__host__ __device__ static __forceinline__ T
	cscb (T radcurv, T x, T y, T phi, T rs, T rc, T& t, T& u, T& v)
	{
		const T r2 = T (2) * radcurv;
		const T sq2 = T (4) * radcurv * radcurv;
		const T a = x + rs;
		const T b = y - rc;
		const T u1 = S::sqrt (a * a + b * b);
		if (u1 < r2) { return S::huge(); }
		const T theta = myAtan2 (b, a);
		u = S::sqrt (u1 * u1 - sq2);
		const T alpha = myAtan2 (r2, u);
		t = mod2pi (theta + alpha);
		v = mod2pi (t - phi);
		return radcurv * (t + v) + u;
	}

	//! C Cu | Cu C
	__host__ __device__ static __forceinline__ T
	ccu_cuc (T radcurv, T x, T y, T phi, T rs, T rc, T& t, T& u, T& v)
	{
		const T r2 = T (2) * radcurv, r4 = T (4) * radcurv;
		const T a = x + rs;
		const T b = y - rc;
		if (S::abs (a) < eps3() && S::abs (b) < eps3()) { return S::huge(); }
		const T u1 = S::sqrt (a * a + b * b);
		if (u1 > r4) { return S::huge(); }
		const T theta = myAtan2 (b, a);
		T alpha;
		if (u1 > r2)
		{
			alpha = S::acos ((u1 / T (2) - radcurv) / r2);
			t			= mod2pi (mpiDiv2() + theta - alpha);
			u			= mod2pi (mpi() - alpha);
			v			= mod2pi (phi - t + T (2) * u);
		}
		else
		{
			alpha = S::acos ((u1 / T (2) + radcurv) / r2);
			t			= mod2pi (mpiDiv2() + theta + alpha);
			u			= mod2pi (alpha);
			v			= mod2pi (phi - t + T (2) * u);
		}
		return radcurv * (T (2) * u + t + v);
	}

	//! C | Cu Cu | C
	__host__ __device__ static __forceinline__ T
	c_cucu_c (T radcurv, T x, T y, T phi, T rs, T rc, T& t, T& u, T& v)
	{
		const T r2 = T (2) * radcurv;
		const T sqr = radcurv * radcurv;
		const T sqr2 = T (4) * radcurv * radcurv;
		const T a = x + rs;
		const T b = y - rc;
		if (S::abs (a) < eps3() && S::abs (b) < eps3()) { return S::huge(); }
		const T u1 = S::sqrt (a * a + b * b);
		if (u1 > T (6) * radcurv) { return S::huge(); }
		const T theta = myAtan2 (b, a);
		const T va1 = (T (5) * sqr - u1 * u1 / T (4)) / sqr2;
		if (va1 < T (0) || va1 > T (1)) { return S::huge(); }
		u = S::acos (va1);
		const T va2 = S::sin (u);
		const T alpha = S::asin (r2 * va2 / u1);
		t = mod2pi (mpiDiv2() + theta + alpha);
		v = mod2pi (t - phi);
		return radcurv * (T (2) * u + t + v);
	}

	//! C | C2 S C, first form
	__host__ __device__ static __forceinline__ T
	c_c2sca (T radcurv, T x, T y, T phi, T rs, T rc, T& t, T& u, T& v)
	{
		const T r2 = T (2) * radcurv;
		const T sq2 = T (4) * radcurv * radcurv;
		const T a = x - rs;
		const T b = y + rc;
		const T u1 = S::sqrt (a * a + b * b);
		if (u1 < r2) { return S::huge(); }
		const T theta = myAtan2 (b, a);
		u = S::sqrt (u1 * u1 - sq2) - r2;
		if (u < T (0)) { return S::huge(); }
		const T alpha = myAtan2 (r2, u + r2);
		t = mod2pi (mpiDiv2() + theta + alpha);
		v = mod2pi (t + mpiDiv2() - phi);
		return radcurv * (t + mpiDiv2() + v) + u;
	}

	//! C | C2 S C, second form
	__host__ __device__ static __forceinline__ T
	c_c2scb (T radcurv, T x, T y, T phi, T rs, T rc, T& t, T& u, T& v)
	{
		const T r2 = T (2) * radcurv;
		const T a = x + rs;
		const T b = y - rc;
		const T u1 = S::sqrt (a * a + b * b);
		if (u1 < r2) { return S::huge(); }
		const T theta = myAtan2 (b, a);
		t = mod2pi (mpiDiv2() + theta);
		u = u1 - r2;
		v = mod2pi (phi - t - mpiDiv2());
		return radcurv * (t + mpiDiv2() + v) + u;
	}

	//! C | C2 S C2 | C
	__host__ __device__ static __forceinline__ T
	c_c2sc2_c (T radcurv, T x, T y, T phi, T rs, T rc, T& t, T& u, T& v)
	{
		const T r2 = T (2) * radcurv, r4 = T (4) * radcurv;
		const T sq2 = T (4) * radcurv * radcurv;
		const T a = x + rs;
		const T b = y - rc;
		const T u1 = S::sqrt (a * a + b * b);
		if (u1 < r4) { return S::huge(); }
		const T theta = myAtan2 (b, a);
		u = S::sqrt (u1 * u1 - sq2) - r4;
		if (u < T (0)) { return S::huge(); }
		const T alpha = myAtan2 (r2, u + r4);
		t = mod2pi (mpiDiv2() + theta + alpha);
		v = mod2pi (t - phi);
		return radcurv * (t + mpi() + v) + u;
	}

	//! C C | C
	__host__ __device__ static __forceinline__ T
	cc_c (T radcurv, T x, T y, T phi, T rs, T rc, T& t, T& u, T& v)
	{
		const T r2 = T (2) * radcurv, r4 = T (4) * radcurv;
		const T sqr = radcurv * radcurv;
		const T a = x - rs;
		const T b = y + rc;
		if (S::abs (a) < eps3() && S::abs (b) < eps3()) { return S::huge(); }
		const T u1 = S::sqrt (a * a + b * b);
		if (u1 > r4) { return S::huge(); }
		const T theta = myAtan2 (b, a);
		u = S::acos ((T (8) * sqr - u1 * u1) / (T (8) * sqr));
		T va = S::sin (u);
		if (S::abs (va) < T (0.001)) { va = T (0); }
		if (S::abs (va) < T (0.001) && S::abs (u1) < T (0.001)) { return S::huge(); }
		const T alpha = S::asin (r2 * va / u1);
		t = mod2pi (mpiDiv2() - alpha + theta);
		v = mod2pi (t - u - phi);
		return radcurv * (t + u + v);
	}

	//! C S C2 | C, first form
	__host__ __device__ static __forceinline__ T
	csc2_ca (T radcurv, T x, T y, T phi, T rs, T rc, T& t, T& u, T& v)
	{
		const T r2 = T (2) * radcurv;
		const T sq2 = T (4) * radcurv * radcurv;
		const T a = x - rs;
		const T b = y + rc;
		const T u1 = S::sqrt (a * a + b * b);
		if (u1 < r2) { return S::huge(); }
		const T theta = myAtan2 (b, a);
		u = S::sqrt (u1 * u1 - sq2) - r2;
		if (u < T (0)) { return S::huge(); }
		const T alpha = myAtan2 (u + r2, r2);
		t = mod2pi (mpiDiv2() + theta - alpha);
		v = mod2pi (t - mpiDiv2() - phi);
		return radcurv * (t + mpiDiv2() + v) + u;
	}

	//! C S C2 | C, second form
	__host__ __device__ static __forceinline__ T
	csc2_cb (T radcurv, T x, T y, T phi, T rs, T rc, T& t, T& u, T& v)
	{
		const T r2 = T (2) * radcurv;
		const T a = x + rs;
		const T b = y - rc;
		const T u1 = S::sqrt (a * a + b * b);
		if (u1 < r2) { return S::huge(); }
		const T theta = myAtan2 (b, a);
		t = mod2pi (theta);
		u = u1 - r2;
		v = mod2pi (-t - mpiDiv2() + phi);
		return radcurv * (t + mpiDiv2() + v) + u;
	}

 public:
	/*!
	 * Solves the point-to-point Reeds-Shepp problem.
	 * @param x0 Initial abscissa.
	 * @param y0 Initial ordinate.
	 * @param th0 Initial heading.
	 * @param x1 Final abscissa.
	 * @param y1 Final ordinate.
	 * @param th1 Final heading.
	 * @param params `params[0]` is the maximum curvature.
	 * @return The shortest of the 48 maneuvers.
	 */
	__host__ __device__ static __forceinline__ Solution
	solve (T x0, T y0, T th0, T x1, T y1, T th1, const T* params)
	{
		const T radcurv = T (1) / params[0];

		// ---- change of coordinates --------------------------------------------
		const T dx		= x1 - x0;
		const T dy		= y1 - y0;
		const T theta = myAtan2 (dy, dx);
		const T alpha = theta - th0;
		const T vard	= S::sqrt (dx * dx + dy * dy);
		const T x			= S::cos (alpha) * vard;
		const T y			= S::sin (alpha) * vard;
		const T phi		= th1 - th0;

		const T sphi = S::sin (phi);
		const T cphi = S::cos (phi);

		const T ap = radcurv * sphi;
		const T am = -radcurv * sphi;
		const T b1 = radcurv * (cphi - T (1));
		const T b2 = radcurv * (cphi + T (1));

		Best best;
		best.len = S::huge();
		best.t = best.u = best.v = T (0);
		best.man = 0;

		T var, tn, un, vn;

		// ---- C | C | C (maneuvers 1-4) ---------------------
		// srcCC assigns maneuver 1 unconditionally, even when it is invalid.
		var = c_c_c (radcurv, x, y, phi, ap, b1, tn, un, vn);
		best.len = var; best.t = tn; best.u = un; best.v = vn; best.man = 1;
		// Maneuver 2 is switched off in srcCC/rs.cc by `if (var < length && false)`.
		// Kept here, disabled, so the GPU matches the CPU exactly.
		// var = c_c_c (radcurv, -x, y, -phi, am, b1, tn, un, vn);
		// keep (best, var, tn, un, vn, 2);
		var = c_c_c (radcurv, x, -y, -phi, am, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 3);
		// Maneuver 4 is switched off in srcCC/rs.cc by `if (var < length && false)`.
		// Kept here, disabled, so the GPU matches the CPU exactly.
		// var = c_c_c (radcurv, -x, -y, phi, ap, b1, tn, un, vn);
		// keep (best, var, tn, un, vn, 4);

		// ---- C | C C (maneuvers 5-8) -----------------------
		var = c_cc (radcurv, x, y, phi, ap, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 5);
		var = c_cc (radcurv, -x, y, -phi, am, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 6);
		var = c_cc (radcurv, x, -y, -phi, am, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 7);
		var = c_cc (radcurv, -x, -y, phi, ap, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 8);

		// ---- C S C (a) (maneuvers 9-12) ---------------------
		var = csca (radcurv, x, y, phi, ap, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 9);
		var = csca (radcurv, x, -y, -phi, am, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 10);
		var = csca (radcurv, -x, y, -phi, am, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 11);
		var = csca (radcurv, -x, -y, phi, ap, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 12);

		// ---- C S C (b) (maneuvers 13-16) ---------------------
		var = cscb (radcurv, x, y, phi, ap, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 13);
		var = cscb (radcurv, x, -y, -phi, am, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 14);
		var = cscb (radcurv, -x, y, -phi, am, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 15);
		var = cscb (radcurv, -x, -y, phi, ap, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 16);

		// ---- C Cu | Cu C (maneuvers 17-20) -------------------
		var = ccu_cuc (radcurv, x, y, phi, ap, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 17);
		var = ccu_cuc (radcurv, x, -y, -phi, am, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 18);
		var = ccu_cuc (radcurv, -x, y, -phi, am, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 19);
		var = ccu_cuc (radcurv, -x, -y, phi, ap, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 20);

		// ---- C | Cu Cu | C (maneuvers 21-24) -----------------
		var = c_cucu_c (radcurv, x, y, phi, ap, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 21);
		var = c_cucu_c (radcurv, x, -y, -phi, am, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 22);
		var = c_cucu_c (radcurv, -x, y, -phi, am, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 23);
		var = c_cucu_c (radcurv, -x, -y, phi, ap, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 24);

		// ---- C | C2 S C (a) (maneuvers 25-28) ----------------
		var = c_c2sca (radcurv, x, y, phi, ap, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 25);
		var = c_c2sca (radcurv, x, -y, -phi, am, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 26);
		var = c_c2sca (radcurv, -x, y, -phi, am, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 27);
		var = c_c2sca (radcurv, -x, -y, phi, ap, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 28);

		// ---- C | C2 S C (b) (maneuvers 29-32) ----------------
		var = c_c2scb (radcurv, x, y, phi, ap, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 29);
		var = c_c2scb (radcurv, x, -y, -phi, am, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 30);
		var = c_c2scb (radcurv, -x, y, -phi, am, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 31);
		var = c_c2scb (radcurv, -x, -y, phi, ap, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 32);

		// ---- C | C2 S C2 | C (maneuvers 33-36) ---------------
		var = c_c2sc2_c (radcurv, x, y, phi, ap, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 33);
		var = c_c2sc2_c (radcurv, x, -y, -phi, am, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 34);
		var = c_c2sc2_c (radcurv, -x, y, -phi, am, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 35);
		var = c_c2sc2_c (radcurv, -x, -y, phi, ap, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 36);

		// ---- C C | C (maneuvers 37-40) -----------------------
		var = cc_c (radcurv, x, y, phi, ap, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 37);
		var = cc_c (radcurv, x, -y, -phi, am, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 38);
		var = cc_c (radcurv, -x, y, -phi, am, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 39);
		var = cc_c (radcurv, -x, -y, phi, ap, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 40);

		// ---- C S C2 | C (a) (maneuvers 41-44) ----------------
		var = csc2_ca (radcurv, x, y, phi, ap, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 41);
		var = csc2_ca (radcurv, x, -y, -phi, am, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 42);
		var = csc2_ca (radcurv, -x, y, -phi, am, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 43);
		var = csc2_ca (radcurv, -x, -y, phi, ap, b1, tn, un, vn);
		keep (best, var, tn, un, vn, 44);

		// ---- C S C2 | C (b) (maneuvers 45-48) ----------------
		var = csc2_cb (radcurv, x, y, phi, ap, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 45);
		var = csc2_cb (radcurv, x, -y, -phi, am, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 46);
		var = csc2_cb (radcurv, -x, y, -phi, am, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 47);
		var = csc2_cb (radcurv, -x, -y, phi, ap, b2, tn, un, vn);
		keep (best, var, tn, un, vn, 48);
		Solution out;
		out.len = best.len;
		out.t		= best.t * radcurv;
		out.u		= best.u * radcurv;
		out.v		= best.v * radcurv;
		out.man = best.man;
		return out;
	}

	/*!
	 * Length-only entry point, which is all the DP needs.
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

#endif	// MPDP_RS_CUH
