/**
 * @file dubins.hh
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief This file contains the function for solving the PP Dubins.
 */

#ifndef DUBINS_HH
#define DUBINS_HH

// Library includes
#include <curve.hh>
#include <utils.hh>

#define DUBINS_DEFAULT_KMAX 0.01

// System includes
#include <cmath>
#include <limits>
#include <fstream>

class Dubins : public Curve {
 public:
	// From LaValle's book the order is LRL RLR LSL LSR RSL RSR
	//! The possible types of Dubins.
	enum D_TYPE {
		INVALID,	///< If a Dubins is not valid, possibly not initialized or no solution.
		LRL,			///< Dubins is of type left-right-left.
		RLR,			///< Dubins is of type right-left-right.
		LSL,			///< Dubins is of type left-straight-left.
		LSR,			///< Dubins is of type left-straight-right.
		RSL,			///< Dubins is of type right-straight-left.
		RSR				///< Dubins is of type right-straight-right.
	};

	//! The possible types of Dubins as strings.
	std::string D_TYPE_STR[7] = {"INVALID", "LRL", "RLR", "LSL", "LSR", "RSL", "RSR"};

 private:
	//! The type of the Dubins path.
	D_TYPE _type;
	//! The maximum curvature and the curvature for each part of the Dubins.
	K_T _kmax = 0.0, _k1 = 0.0, _k2 = 0.0, _k3 = 0.0;
	//! The lengths of each part of the Dubins.
	LEN_T _s1 = 0.0, _s2 = 0.0, _s3 = 0.0;

	/*!
	 * Function to standardize the components. Credit to Marco Frego & Paolo Bevilacqua
	 */
	void
	scaleToStandard (Angle& phi, real_type& lambda, Angle& sth0, Angle& sth1, K_T& sKmax);

	/*!
	 * Given the standardized version, compute the best word. Credit to Marco Frego & Paolo
	 * Bevilacqua.
	 * @param th0 The initial standardized angle.
	 * @param th1 The final standardized angle.
	 * @param lambda A multiplier.
	 * @param sKmax The standardized curvature.
	 */
	void
	computeBest (Angle th0, Angle th1, real_type lambda, K_T& sKmax);

	/*!
	 * Function to solve the Dubins curve.
	 */
	void
	solve() override
	{
		real_type lambda;
		K_T sKmax;
		Angle phi, sth0, sth1;
		this->scaleToStandard (phi, lambda, sth0, sth1, sKmax);
		this->computeBest (sth0, sth1, lambda, sKmax);
	}

	/**
	 * Struct to store the input data for the Dubins computation of the maneuvers.
	 */
	static struct Man_input_data {
		real_type th0;
		real_type th1;
		K_T sKmax;
		real_type invK;
		real_type cos_0;
		real_type cos_1;
		real_type Ksq;
		real_type dcos;
		real_type dcos2;
		real_type dsin;
		real_type scos;
		real_type ssin;
		real_type dth;
	} Man_input_data;

	/**
	 * Struct to store the output data for the Dubins computation of the maneuvers.
	 */
	static struct Man_output_data {
		real_type len;
		real_type ss1;
		real_type ss2;
		real_type ss3;
		real_type sk1;
		real_type sk2;
		real_type sk3;
		D_TYPE type;
	} Man_output_data;

 public:
	/*!
	 * Void constructor to initialize a Dubins object.
	 */
	Dubins() : Curve (CURVE_TYPE::DUBINS), _type (D_TYPE::INVALID), _kmax (0) {}

	/*!
	 * Constructor to initialize a Dubins object with an initial and a final
	 * `Configuration2` and additional possible parameters.
	 * @param ci The initial `Configuration2`.
	 * @param cf The final `Configuration2`
	 * @param params Additional parameters to pass. Default is `nullptr`, in such case
	 * DUBINS_DEFAULT_KMAX==0.01 is used.
	 */
	Dubins (Configuration2 ci, Configuration2 cf, std::vector<real_type> params = {})
			: Curve (ci, cf, CURVE_TYPE::DUBINS, params), _type (D_TYPE::INVALID)
	{
		if (params.empty()) { this->_kmax = DUBINS_DEFAULT_KMAX; }
		else { this->_kmax = params[0]; }
		solve();
	}

	/*!
	 * Constructor to initialize a Dubins object with an initial and a final
	 * `Configuration2` and additional possible parameters.
	 * @param x0 The initial abscissa.
	 * @param y0 The final ordinate.
	 * @param th0 The initial angle.
	 * @param x1 The initial abscissa.
	 * @param y1 The final ordinate.
	 * @param th1 The final angle.
	 * @param params Additional parameters to pass. Default is `nullptr`, in such case
	 * DUBINS_DEFAULT_KMAX==0.01 is used.
	 */
	Dubins (
			real_type x0,
			real_type y0,
			Angle th0,
			real_type x1,
			real_type y1,
			real_type th1,
			std::vector<real_type> params = {})
			: Curve (
						Configuration2 (x0, y0, th0),
						Configuration2 (x1, y1, th1),
						CURVE_TYPE::DUBINS,
						params),
				_type (D_TYPE::INVALID)
	{
		if (params.empty()) { this->_kmax = DUBINS_DEFAULT_KMAX; }
		else { this->_kmax = params[0]; }
		solve();
	}

	/*!
	 * Constructor to initialize a Dubins object with an initial and a final
	 * `Configuration2` and additional possible parameters.
	 * @param ci The initial `Configuration2`.
	 * @param cf The final `Configuration2`
	 * @param kmax The curvature of the Dubins parts.
	 */
	Dubins (Configuration2 ci, Configuration2 cf, real_type kmax)
			: Curve (ci, cf, CURVE_TYPE::DUBINS), _type (D_TYPE::INVALID), _kmax (kmax)
	{
		solve();
	}

	/*!
	 * Constructor to initialize a Dubins object with an initial and a final
	 * `Configuration2` and additional possible parameters.
	 * @param ci The initial `Configuration2`.
	 * @param cf The final `Configuration2`
	 * @param kmax The curvature of the Dubins parts.
	 */
	Dubins (
			Configuration2 ci, Configuration2 cf, std::vector<real_type> params, D_TYPE type)
			: Curve (ci, cf, CURVE_TYPE::DUBINS, params), _type (D_TYPE::INVALID),
				_kmax (params[0])
	{
		this->comp_man (type);
	}

	K_T
	kmax() const
	{
		return this->_kmax;
	}	 ///< Returns the maximum curvature.
	inline K_T
	k (int id) const
	{
		switch (id)
		{
			case 1:
				return this->k1();
				break;
			case 2:
				return this->k2();
				break;
			case 3:
				return this->k3();
				break;
			default:
				return std::numeric_limits<K_T>::quiet_NaN();
		}
	}

	/**
	 * @brief Returns the curvature of the first part of the Dubins.
	 * @return The curvature of the first part of the Dubins.
	 */
	inline K_T
	k1() const
	{
		return this->_k1;
	}

	/**
	 * @brief Returns the curvature of the middle part of the Dubins.
	 * @return The curvature of the middle part of the Dubins.
	 */
	inline K_T
	k2() const
	{
		return this->_k2;
	}

	/**
	 * @brief Returns the curvature of the final part of the Dubins.
	 * @return The curvature of the final part of the Dubins.
	 */
	inline K_T
	k3() const
	{
		return this->_k3;
	}

	/**
	 * @brief Returns the length of the first part of the Dubins.
	 * @return The length of the first part of the Dubins.
	 */
	inline LEN_T
	s1() const
	{
		return this->_s1;
	}

	/**
	 * @brief Returns the length of the middle part of the Dubins.
	 * @return The length of the middle part of the Dubins.
	 */
	inline LEN_T
	s2() const
	{
		return this->_s2;
	}

	/**
	 * @brief Returns the length of the final part of the Dubins.
	 * @return The length of the final part of the Dubins.
	 */
	inline LEN_T
	s3() const
	{
		return this->_s3;
	}

	/**
	 * @brief Returns the length of the Dubins part.
	 * @param id 1,2,3 to select the part of the Dubins.
	 * @return The length of the selected part of the Dubins. If the id is not valid, it
	 * returns `std::numeric_limits<LEN_T>::quiet_NaN()`.
	 */
	inline LEN_T
	L (int id) const
	{
		switch (id)
		{
			case 1:
				return this->s1();
				break;
			case 2:
				return this->s2();
				break;
			case 3:
				return this->s3();
				break;
			default:
				return std::numeric_limits<LEN_T>::quiet_NaN();
		}
	}

	/**
	 * @brief Returns the length of the whole Dubins maneuver.
	 * @return The length of the Dubins.
	 */
	inline LEN_T
	l() const override
	{
		return (this->s1() + this->s2() + this->s3());
	}

	/**
	 * @brief Returns the type of the Dubins.
	 * @return The type of the Dubins.
	 */
	D_TYPE
	dtype() const { return this->_type; }

	/**
	 * @brief Sets the maximum curvature and returns the new set value.
	 * @param kmax The new maximum curvature.
	 * @return The new maximum curvature.
	 */
	K_T
	kmax (K_T kmax)
	{
		this->_kmax = kmax;
		return this->kmax();
	}

	/**
	 * @brief Sets the curvature of the first part of the Dubins and returns the new set
	 * value
	 * @param k1 The new curvature of the first part of the Dubins.
	 * @return The new curvature of the first part of the Dubins.
	 */
	K_T
	k1 (K_T k1)
	{
		this->_k1 = k1;
		return this->k1();
	}

	/**
	 * @brief Sets the curvature of the middle part of the Dubins and returns the new set
	 * value
	 * @param k2 The new curvature of the middle part of the Dubins.
	 * @return The new curvature of the middle part of the Dubins.
	 */
	K_T
	k2 (K_T k2)
	{
		this->_k2 = k2;
		return this->k2();
	}

	/**
	 * @brief Sets the curvature of the final part of the Dubins and returns the new set
	 * value
	 * @param k3 The new curvature of the final part of the Dubins.
	 * @return The new curvature of the final part of the Dubins.
	 */
	K_T
	k3 (K_T k3)
	{
		this->_k3 = k3;
		return this->k3();
	}

	/**
	 * @brief Sets the length of the first part of the Dubins and returns the new set value
	 * @param s1 The new length of the first part of the Dubins.
	 * @return The new length of the first part of the Dubins.
	 */
	LEN_T
	s1 (LEN_T s1)
	{
		this->_s1 = s1;
		return this->s1();
	}

	/**
	 * @brief Sets the length of the middle part of the Dubins and returns the new set value
	 * @param s2 The new length of the middle part of the Dubins.
	 * @return The new length of the middle part of the Dubins.
	 */
	LEN_T
	s2 (LEN_T s2)
	{
		this->_s2 = s2;
		return this->s2();
	}
	/**
	 * @brief Sets the length of the final part of the Dubins and returns the new set value
	 * @param s2 The new length of the final part of the Dubins.
	 * @return The new length of the final part of the Dubins.
	 */
	LEN_T
	s3 (LEN_T s3)
	{
		this->_s3 = s3;
		return this->s3();
	}

	/**
	 * @brief Sets the word of the Dubins and returns the new set value
	 * @param type The new word of the Dubins.
	 * @return The new word of the Dubins.
	 */
	D_TYPE
	dtype (D_TYPE type)
	{
		this->_type = type;
		return this->dtype();
	}

	// LRL RLR LSL LSR RSL RSR
	bool
	comp_LRL (const struct Man_input_data& data_in, struct Man_output_data& data_out) const;
	bool
	comp_RLR (const struct Man_input_data& data_in, struct Man_output_data& data_out) const;
	bool
	comp_LSL (const struct Man_input_data& data_in, struct Man_output_data& data_out) const;
	bool
	comp_LSR (const struct Man_input_data& data_in, struct Man_output_data& data_out) const;
	bool
	comp_RSL (const struct Man_input_data& data_in, struct Man_output_data& data_out) const;
	bool
	comp_RSR (const struct Man_input_data& data_in, struct Man_output_data& data_out) const;

	LEN_T
	comp_man (D_TYPE type);

	std::vector<std::vector<double>>
	split_wise() override;

	/*!
	 * Function to print the stringy word of the Dubins.
	 * @return A string containing the word of the computed Dubins.
	 */
	inline std::string
	man_to_string() const
	{
		return D_TYPE_STR[this->dtype()];
	}

	/*!
	 * Function to print the most essential info about `Dubins`.
	 * @param str An additional string to add at the beginning.
	 * @return A `std::stringstream` object containing the data of `Dubins`.
	 */
	std::stringstream
	to_string (const std::string& str = "")
	{
		std::stringstream out;
		out << (str.empty() ? "" : str + " ") << "c0: " << this->ci()->to_string().str()
				<< "\tc1: " << this->cf()->to_string().str() << "\tk: " << this->kmax()
				<< "\tl: " << this->l() << "\ttype: " << this->man_to_string();
		return out;
	}

	/*! This function overload the << operator so to print with `std::cout` the most
		 essential info about the `Dubins`.
			\param[in] out The out stream.
			\param[in] data The Dubins to print.
			\returns An output stream to be printed.
	*/
	friend std::ostream&
	operator<< (std::ostream& out, Dubins& data)
	{
		out << data.to_string().str();
		return out;
	}

#ifdef MPDP_DRAW
	void
	draw (
			std::ofstream& file,
			std::string label = "",
			size_t width			= 8,
			size_t height			= 8,
			bool solve				= false,
			bool close				= false,
			bool init					= false);
#endif	// MPDP_DRAW
};

#endif	// DUBINS_HH
