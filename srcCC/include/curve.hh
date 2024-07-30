/**
 * @file curve.hh
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License Agero 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief This file contains the definition of the Curve class, which is the base abstract
 * class for all the curves.
 * @todo The points should be made as smart pointers.
 */

#ifndef CURVE_HH
#define CURVE_HH

// System includes
#include <vector>

// Library includes
#include <configuration.hh>
#ifdef MPDP_DRAW
// #include <clothoids/clothoidAsyPlot.hh>
#endif

enum class CURVE_TYPE { INVALID, DUBINS, DUBINS_ARC, RS };	///< Possible types of CURVE

class Curve {
 private:
	Configuration2 _ci;							 ///< Initial `Configuration`
	Configuration2 _cf;							 ///< Final `Configuration`
	CURVE_TYPE _type;								 ///< Type of curve
	std::vector<real_type> _params;	 ///< Parameters of curve

 public:
	/*!
	 * @brief Void constructor.
	 */
	Curve() : _ci(), _cf(), _type (CURVE_TYPE::INVALID), _params ({}) {}
	/*!
	 * @brief Constructor to only set the type of the curve.
	 */
	Curve (CURVE_TYPE type = CURVE_TYPE::INVALID) : _ci(), _cf(), _type (type), _params ({})
	{}

	/*!
	 * @brief Constructor that takes two `Configuration2` and the type of the curve.
	 * @param ci Initial configuration.
	 * @param cf Final configuration.
	 * @param type Type of the curve.
	 * @param params The parameters of the curve, such as the curvature.
	 */
	Curve (
			Configuration2 ci,
			Configuration2 cf,
			CURVE_TYPE type								= CURVE_TYPE::INVALID,
			std::vector<real_type> params = {})
			: _ci (ci), _cf (cf), _type (type), _params (params)
	{}

	/**
	 * @brief Returns a pointer to the initial `Configuration2`.
	 * @return A pointer to the initial `Configuration2`.
	 */
	Configuration2*
	ci()
	{
		return &(this->_ci);
	}

	/**
	 * @brief Returns a pointer to the final `Configuration2`.
	 * @return A pointer to the final `Configuration2`.
	 */
	Configuration2*
	cf()
	{
		return &(this->_cf);
	}

	/**
	 * @brief Sets the initial `Configuration2`.
	 * @param ci The initial `Configuration2`.
	 */
	void
	ci (Configuration2 ci)
	{
		this->_ci = ci;
	}

	/**
	 * @brief Sets the final `Configuration2`.
	 * @param cf The final `Configuration2`.
	 */
	void
	cf (Configuration2 cf)
	{
		this->_cf = cf;
	}

	/**
	 * @brief Returns the type of the curve.
	 * @return The type of the curve.
	 */
	CURVE_TYPE
	type() const { return this->_type; }

	/**
	 * @brief Returns the parameters of the curve
	 * @return A vector of `real_type` containing the parameters of the curve.
	 */
	std::vector<real_type>
	params() const
	{
		return this->_params;
	}

	virtual LEN_T
	l() const = 0;	///< Returns the length of the curve.

	virtual void
	solve() = 0;	///< Solves the curve depending on the type.

	virtual std::vector<std::vector<double>>
	split_wise() = 0;	 ///< Splits the curve.

	/**
	 * @brief Splits the curve into `num_split` parts.
	 * @param num_split The number of parts to split the curve into.
	 * @return A vector of `Configuration2` containing the split points.
	 */
	//  virtual std::vector<Configuration2> split(int num_split) = 0;

#ifdef MPDP_DRAW
//  virtual void draw() = 0;                              ///< Draws the curve.
#endif
};

#endif	// CURVE_HH
