/**
 * @file dp.hh
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License Agero 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief This file contains the function for the dynamic programming algorithm.
 */

#ifndef DP_HH
#define DP_HH

#ifndef CUDA_ON

// System includes
#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <omp.h>

// Library includes
#include <utils.hh>
#include <typedefs.hh>
#include <configuration.hh>
#include <dubins.hh>
#include <rs.hh>

#ifdef DEBUG
#define printMatrix(type)                                                \
	if (type == 0) { std::cout << "Angle table: " << std::endl; }          \
	else if (type == 1) { std::cout << "Length table: " << std::endl; }    \
	for (int i = 0; i < MATRIX.size(); i++)                                \
	{                                                                      \
		for (int j = 0; j < MATRIX[i].size(); j++)                           \
		{                                                                    \
			if (type == 0)                                                     \
			{                                                                  \
				std::cout << std::setprecision (4) << MATRIX[i][j].th() << "\t"; \
			}                                                                  \
			else if (type == 1)                                                \
			{                                                                  \
				std::cout << std::setprecision (4) << MATRIX[i][j].l() << "\t";  \
			}                                                                  \
		}                                                                    \
		std::cout << std::endl;                                              \
	}                                                                      \
	std::cout << std::endl << std::endl;
#else
#define printMatrix(type)
#endif

static K_T Kmax = DUBINS_DEFAULT_KMAX;
#define MATRIX this->matrix

/*!
 * This function returns (up to) two circles through two points, given the radius. Credit
 * to Marco Frego & Paolo Bevilacqua.
 * @param x1 The initial abscissa coordinate.
 * @param y1 The initial ordinate coordinate.
 * @param x2 The final abscissa coordinate.
 * @param y2 The final ordinate coordinate.
 * @param r The radius of the possible circles.
 * @param XC The abscissas of points to use to calculate the tangents.
 * @param YC The abscissas of points to use to calculate the tangents.
 */
void
circles (
		double x1,
		double y1,
		double x2,
		double y2,
		double r,
		std::vector<double>& XC,
		std::vector<double>& YC);

//////////////////////////////////////////////////////////////////////////////////////////

class DP {
private:
	class Cell {
	private:
		Angle _th;	///< Angle of the final point of the point to point curve.
		LEN_T _l;		///< Length of the point to point curve.
		int _next;	///< Id to the next cell for dynamic programming.

	public:
		/*!
		 * Default void constructor which returns a cell initialized with ANGLE::FREE, max
		 * length and -1 as next cell.
		 */
		Cell() : _th (ANGLE::FREE), _l (std::numeric_limits<LEN_T>::max()), _next (-1) {}

		/*!
		 * Constructor that takes in input an angle, a length and the next id and returns a
		 * DP::Cell.
		 * @param th The initial angle of the point to point curve.
		 * @param l The length of the point to point curve. Default is `LEN_T` max value.
		 * @param next The next id of the cell. Default is -1.
		 */
		explicit Cell (Angle th, LEN_T l = std::numeric_limits<LEN_T>::max(), int next = -1)
				: _th (th), _l (l), _next (next)
		{}

		/*!
		 * Returns the angle.
		 * @return the angle.
		 */
		Angle
		th() const
		{
			return this->_th;
		}
		/*!
		 * Returns the length.
		 * @return the length.
		 */
		LEN_T
		l() const { return this->_l; }
		/*!
		 * Returns the next it.
		 * @return the next it.
		 */
		int
		next() const
		{
			return this->_next;
		}

		/*!
		 * Sets the new angle.
		 * @param th The new angle to be set.
		 * @return the new set angle.
		 */
		Angle
		th (Angle th)
		{
			this->_th = th;
			return this->th();
		}
		/*!
		 * Sets the new length.
		 * @param th The new length to be set.
		 * @return the new set length.
		 */
		LEN_T
		l (LEN_T l)
		{
			this->_l = l;
			return this->l();
		}
		/*!
		 * Sets the new next id.
		 * @param th The new next id to be set.
		 * @return the new set next id.
		 */
		int
		next (int next)
		{
			this->_next = next;
			return this->next();
		}

		/*!
		 * Creates a deep copy of a cell to `this`.
		 * @param d The cell to copy from.
		 * @return `*this`.
		 */
		Cell&
		copy (const Cell& d)
		{
			this->th (d.th());
			this->l (d.l());
			this->next (d.next());

			return *this;
		}

		/*!
		 * Overrides the assign operator (=) to make a deep copy of a cell to `this`.
		 * @param d The cell to copy from.
		 * @return `*this`.
		 */
		Cell&
		operator= (const Cell& d)
		{
			this->copy (d);
			return *this;
		}

		/*!
		 * Function to print the most essential info about `DP::Cell`.
		 * @param pretty An additional truth value to print a prettier version. Default is
		 * `false`.
		 * @return A `std::stringstream` object containing the data of `DP::Cell`.
		 */
		std::stringstream
		to_string (bool pretty = false) const
		{
			std::stringstream out;
			if (pretty) { out << "th: " << this->th() << " l: " << this->l(); }
			else
			{
				out << "<" << (Angle)(this->th() * 1.0) << ", " << (LEN_T)(this->l()) << ">";
			}
			return out;
		}

		/*! This function overrides the << operator so to print with `std::cout` the most
			 essential info about the `DP::Cell`.
				\param[in] out The out stream.
				\param[in] data The `DP::Cell` to print.
				\returns An output stream to be printed.
		*/
		friend std::ostream&
		operator<< (std::ostream& out, const Cell& data)
		{
			out << data.to_string().str();
			return out;
		}
	};

	/// The matrix used for the dynamic programming algorithm.
	std::vector<std::vector<Cell>> matrix;

	/*!
	 * This function sets the possible angles to consider for each point in the matrix. It
	 * is specialized for the first iteration.
	 * @param discr The number of samples to consider for each angle.
	 * @param fixedAngles A vector of bools stating which angles should not be changed.
	 * @param points A vector containing all the points through the multi-point curve should
	 * go.
	 */
	void
	setSamplingAngles (
			int discr,
			const std::vector<bool>& fixedAngles,
			const std::vector<Configuration2>& points);

	/*!
	 * This function sets the possible angles to consider for each point in the matrix. It
	 * is specialized for the first iteration.
	 * @param points A vector containing all the points through the multi-point curve should
	 * go.
	 * @param fixedAngles A vector of `bool` stating which angles should not be changed.
	 * @param hrange The value to move from the previous computed angle to find new sampling
	 * angles.
	 * @param hn The number to cycle to add new samples (set to half the value of new
	 * samplings since the added samples are mirrored w.r.t. the previously found angle).
	 */
	void
	setSamplingAngles (
			const std::vector<Configuration2>& points,
			const std::vector<bool>& fixedAngles,
			double hrange,
			int hn);

	/*!
	 * Function to compute the backward phase once all the matrix is filled.
	 * @param points The list of points.
	 * @return Returns the best angles to go through to compute the shortest multi-point
	 * path.
	 */
	std::pair<LEN_T, std::vector<Angle>>
	bestAngles (std::vector<Configuration2>* points = nullptr);

	/*!
	 * Given two points, this function returns some particular angles, for example the angle
	 * that would allow to connect the two points with a straight line. Credit to Marco
	 * Frego & Paolo Bevilacqua.
	 * @param i The ith point to consider.
	 * @param thPrev Where to save the angle for the (i-1)th point.
	 * @param thCur Where to save the angle for the ith point.
	 * @param points All the points. (Simpler and cheaper than passing two `Configuration2`)
	 */
	void
	guessInitialAngles (
			const uint i,
			std::vector<Cell>& thPrev,
			std::vector<Cell>& thCur,
			const std::vector<Configuration2>& points);

	/*!
	 * Function that computes the point to point curves and finds the best values for each
	 * initial angle considering the rest of the curve.
	 * @param points The list of points which the multi-point path should go through.
	 * @param params Additional parameters to pass to the curve constructor.
	 * @return A pair where the first element is the computed length of the curve and the
	 * second one is a vector containing the best angles.
	 */
	template <class CurveT = Dubins>
	std::pair<LEN_T, std::vector<Angle>>
	solveDPInner (std::vector<Configuration2>& points, std::vector<real_type>& params)
	{
		static_assert (
				std::is_base_of<Curve, CurveT>::value, "CurveT must be a subclass of Curve");

		printMatrix (0);

		for (uint idx = (points.size() - 1); idx > 0; --idx)
		{	 // Cycle between all points starting from the last one
			Configuration2* c0 = &points[idx - 1];
			Configuration2* c1 = &points[idx];

			#pragma omp parallel for
			for (int i = 0; i < (int)(MATRIX[idx - 1].size()); ++i)
			{	 // Consider the previous angle
				int bestJ		 = -1;
				double bestL = std::numeric_limits<LEN_T>::max();

				for (uint j = 0; j < MATRIX[idx].size(); ++j)
				{	 // Consider the next angle
					// Compute Dubins
					CurveT curve = CurveT (
							c0->x(), c0->y(), MATRIX[idx - 1][i].th(), c1->x(), c1->y(),
							MATRIX[idx][j].th(), params);
					LEN_T curL = curve.l();
					if (idx == (points.size() - 1) && i == 1)
					{
						COUT (MATRIX[idx - 1][i])
						COUT (MATRIX[idx][j])
						COUT (dub)
						COUT (curL)
					}

					if (idx < (points.size() - 1)) { curL += MATRIX[idx][j].l(); }

					if (curL < bestL)
					{
						bestL = curL;
						bestJ = j;
					}
				}
				MATRIX[idx - 1][i].l (bestL);
				MATRIX[idx - 1][i].next (bestJ);
			}
		}
		printMatrix (1) return bestAngles (&points);
	}

public:
	DP() {}	 ///< Default constructor.

	/*!
	 * The wrapper to call to solve the dynamic programming problem for multi-point path
	 * finding.
	 * @param points A vector of points the path should go through.
	 * @param discr The number of sampling to consider for each point.
	 * @param fixedAngles A vector stating which angles should not be changed.
	 * @param nRef The number of times the algorithm should be called back in order to
	 * narrow the sampling intervals finding more precise final values.
	 * @param params A list of additional parameters to pass to the curve constructor.
	 * @param saveAngles If set to `true`, then the points angles are changed to the best
	 * angles found at the end of the algorithm, otherwise, they are only returned. Default
	 * is true
	 * @return A pair where the first element is the computed length of the curve and the
	 * second one is a vector containing the best angles.
	 */
	template <class CurveT = Dubins>
	std::pair<LEN_T, std::vector<Angle>>
	solveDP (
			std::vector<Configuration2>& points,
			const std::vector<bool>& fixedAngles,
			std::vector<real_type> params,
			int discr,
			int nRef,
			bool saveAngles = true)
	{
		static_assert (
				std::is_base_of<Curve, CurveT>::value, "CurveT must be a subclass of Curve");
		MATRIX.clear();

		// TODO this should be independent of the params argument
		Kmax = params[0];

		std::pair<LEN_T, std::vector<Angle>> ret;
		std::vector<Angle> bestA;

		std::vector<Configuration2> compPoints;
		for (auto p : points)
		{
			compPoints.emplace_back (Configuration2 (p.x(), p.y(), p.th()));
		}

		// First round
		setSamplingAngles (discr, fixedAngles, compPoints);
		ret		= solveDPInner (compPoints, params);
		bestA = ret.second;

		// Other refinements
		double hrange = 2.0 * m_pi;
		for (int ref = 0; ref < nRef; ++ref)
		{
			COUT (ref)
			hrange = hrange / discr * 1.5;
			setSamplingAngles (compPoints, fixedAngles, hrange, discr / 2);
			printV (bestA)

					ret = solveDPInner<CurveT> (compPoints, params);
			bestA		= ret.second;
		}

		if (saveAngles)
		{
			for (uint i = 0; i < points.size(); i++) { points[i].th (bestA[i]); }
		}

		return ret;
	}
};

#endif
#endif	// DP_HH
