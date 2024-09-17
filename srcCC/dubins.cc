/**
 * @file dubins.cc
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief This file contains the source code for some of the functions for PP Dubins.
 */

#ifndef CUDA_ON
#include <dubins.hh>

/**
 * @brief Function to scale to a standard settings the values.
 */
void
Dubins::scaleToStandard (
		Angle& phi, real_type& lambda, Angle& sth0, Angle& sth1, K_T& sKmax)
{
	// Credit to Marco Frego & Paolo Bevilacqua.
	real_type dx = this->cf()->x() - this->ci()->x();
	real_type dy = this->cf()->y() - this->ci()->y();
	phi					 = atan2 (dy, dx);
	lambda			 = hypot (dx, dy) * 0.5;
	sKmax				 = this->kmax() * lambda;
	sth0				 = mod2pi (this->ci()->th() - phi);
	sth1				 = mod2pi (this->cf()->th() - phi);
}

/**
 *
 * @param man_in_data
 * @param man_out_data
 * @return true if a maneuver could be computed, not necessarily the shortest one.
 */
bool
Dubins::comp_LRL (
		const struct Man_input_data& man_in_data, struct Man_output_data& man_out_data) const
{
	//  std::cout << "Computing LRL" << std::endl;
	real_type C			= -man_in_data.dcos2;
	real_type S			= 2 * man_in_data.sKmax + man_in_data.dsin;
	real_type temp1 = std::atan2 (C, S);
	real_type temp2 = 0.125 * (6 - 4 * man_in_data.Ksq + 2 * man_in_data.dcos -
														 4 * man_in_data.sKmax * man_in_data.dsin);
	if (std::abs (temp2) <= 1)
	{
		real_type t2 = man_in_data.invK * mod2pi (2 * m_pi - std::acos (temp2));
		real_type t1 = man_in_data.invK *
									 mod2pi (-man_in_data.th0 + temp1 + 0.5 * t2 * man_in_data.sKmax);
		real_type t3 =
				man_in_data.invK * mod2pi (-man_in_data.dth + (t2 - t1) * man_in_data.sKmax);
		real_type lc = t1 + t2 + t3;
		//    std::cout << std::setprecision(12) << "LRL len: " << lc << std::endl;
		if (lc < man_out_data.len)
		{
			man_out_data.len	= lc;
			man_out_data.ss1	= t1;
			man_out_data.ss2	= t2;
			man_out_data.ss3	= t3;
			man_out_data.sk1	= 1;
			man_out_data.sk2	= -1;
			man_out_data.sk3	= 1;
			man_out_data.type = D_TYPE::LRL;
		}
		return true;
	}
	return false;
}

/**
 *
 * @param man_in_data
 * @param man_out_data
 * @return true if a maneuver could be computed, not necessarily the shortest one.
 */
bool
Dubins::comp_RLR (
		const struct Man_input_data& man_in_data, struct Man_output_data& man_out_data) const
{
	//  std::cout << "Computing RLR" << std::endl;
	real_type C			= man_in_data.dcos2;
	real_type S			= 2 * man_in_data.sKmax - man_in_data.dsin;
	real_type temp1 = std::atan2 (C, S);
	real_type temp2 = 0.125 * (6 - 4 * man_in_data.Ksq + 2 * man_in_data.dcos +
														 4 * man_in_data.sKmax * man_in_data.dsin);
	if (std::abs (temp2) <= 1)
	{
		real_type t2 = man_in_data.invK * mod2pi (2 * m_pi - std::acos (temp2));
		real_type t1 = man_in_data.invK *
									 mod2pi (man_in_data.th0 - temp1 + 0.5 * t2 * man_in_data.sKmax);
		real_type t3 =
				man_in_data.invK * mod2pi (man_in_data.dth + (t2 - t1) * man_in_data.sKmax);
		real_type lc = t1 + t2 + t3;
		//    std::cout << std::setprecision(12) << "RLR len: " << lc << std::endl;
		if (lc < man_out_data.len)
		{
			man_out_data.len	= lc;
			man_out_data.ss1	= t1;
			man_out_data.ss2	= t2;
			man_out_data.ss3	= t3;
			man_out_data.sk1	= -1;
			man_out_data.sk2	= 1;
			man_out_data.sk3	= -1;
			man_out_data.type = D_TYPE::RLR;
		}
		return true;
	}
	return false;
}

/**
 *
 * @param man_in_data
 * @param man_out_data
 * @return true if a maneuver could be computed, not necessarily the shortest one.
 */
bool
Dubins::comp_LSL (
		const struct Man_input_data& man_in_data, struct Man_output_data& man_out_data) const
{
	//  std::cout << "Computing LSL" << std::endl;
	real_type C			= man_in_data.cos_1 - man_in_data.cos_0;
	real_type S			= 2 * man_in_data.sKmax + man_in_data.dsin;
	real_type temp1 = std::atan2 (C, S);
	real_type temp2 = 2 + 4 * man_in_data.Ksq - 2 * man_in_data.dcos +
										4 * man_in_data.sKmax * man_in_data.dsin;
	//  std::cout << C << " " << S << " " << temp1 << " " << temp2 << std::endl;
	if (temp2 >= 0)
	{
		real_type temp3 = man_in_data.invK * std::sqrt (temp2);
		real_type t1		= man_in_data.invK * mod2pi (temp1 - man_in_data.th0);
		real_type t2		= temp3;
		real_type t3		= man_in_data.invK * mod2pi (man_in_data.th1 - temp1);
		real_type lc		= t1 + t2 + t3;
		//    std::cout << std::setprecision(12) << "LSL len: " << lc << std::endl;
		if (lc < man_out_data.len)
		{
			man_out_data.len	= lc;
			man_out_data.ss1	= t1;
			man_out_data.ss2	= t2;
			man_out_data.ss3	= t3;
			man_out_data.sk1	= 1;
			man_out_data.sk2	= 0;
			man_out_data.sk3	= 1;
			man_out_data.type = D_TYPE::LSL;
		}
		return true;
	}
	return false;
}

/**
 *
 * @param man_in_data
 * @param man_out_data
 * @return true if a maneuver could be computed, not necessarily the shortest one.
 */
bool
Dubins::comp_LSR (
		const struct Man_input_data& man_in_data, struct Man_output_data& man_out_data) const
{
	//  std::cout << "Computing LSR" << std::endl;
	real_type C			= man_in_data.scos;
	real_type S			= 2 * man_in_data.sKmax + man_in_data.ssin;
	real_type temp1 = std::atan2 (-C, S);
	real_type temp2 = -2 + 4 * man_in_data.Ksq + 2 * man_in_data.dcos +
										4 * man_in_data.sKmax * man_in_data.ssin;
	if (temp2 >= 0)
	{
		real_type t2		= man_in_data.invK * std::sqrt (temp2);
		real_type temp3 = -std::atan2 (-2, (t2 * man_in_data.sKmax));
		real_type t1		= man_in_data.invK * mod2pi (-man_in_data.th0 + temp1 + temp3);
		real_type t3		= man_in_data.invK * mod2pi (-man_in_data.th1 + temp1 + temp3);
		real_type lc		= t1 + t2 + t3;
		//    std::cout << std::setprecision(12) << "LSR len: " << lc << std::endl;
		if (lc < man_out_data.len)
		{
			man_out_data.len	= lc;
			man_out_data.ss1	= t1;
			man_out_data.ss2	= t2;
			man_out_data.ss3	= t3;
			man_out_data.sk1	= 1;
			man_out_data.sk2	= 0;
			man_out_data.sk3	= -1;
			man_out_data.type = D_TYPE::LSR;
		}
		return true;
	}
	return false;
}

/**
 *
 * @param man_in_data
 * @param man_out_data
 * @return true if a maneuver could be computed, not necessarily the shortest one.
 */
bool
Dubins::comp_RSL (
		const struct Man_input_data& man_in_data, struct Man_output_data& man_out_data) const
{
	//  std::cout << "Computing RSL" << std::endl;
	real_type C			= man_in_data.scos;
	real_type S			= 2 * man_in_data.sKmax - man_in_data.ssin;
	real_type temp1 = std::atan2 (C, S);
	real_type temp2 = -2 + 4 * man_in_data.Ksq + 2 * man_in_data.dcos -
										4 * man_in_data.sKmax * man_in_data.ssin;
	//  std::cout << std::setprecision(12) << "temp2: " << temp2 << std::endl;
	if (temp2 >= 0)
	{
		real_type t2		= man_in_data.invK * std::sqrt (temp2);
		real_type temp3 = std::atan2 (2, (t2 * man_in_data.sKmax));
		real_type t1		= man_in_data.invK * mod2pi (man_in_data.th0 - temp1 + temp3);
		real_type t3		= man_in_data.invK * mod2pi (man_in_data.th1 - temp1 + temp3);
		real_type lc		= t1 + t2 + t3;
		//    std::cout << std::setprecision(12) << "RSL len: " << lc << std::endl;
		if (lc < man_out_data.len)
		{
			man_out_data.len	= lc;
			man_out_data.ss1	= t1;
			man_out_data.ss2	= t2;
			man_out_data.ss3	= t3;
			man_out_data.sk1	= -1;
			man_out_data.sk2	= 0;
			man_out_data.sk3	= 1;
			man_out_data.type = D_TYPE::RSL;
		}
		return true;
	}
	return false;
}

/**
 *
 * @param man_in_data
 * @param man_out_data
 * @return true if a maneuver could be computed, not necessarily the shortest one.
 */
bool
Dubins::comp_RSR (
		const struct Man_input_data& man_in_data, struct Man_output_data& man_out_data) const
{
	//  std::cout << "Computing RSR" << std::endl;
	real_type C			= man_in_data.cos_0 - man_in_data.cos_1;
	real_type S			= 2 * man_in_data.sKmax - man_in_data.dsin;
	real_type temp1 = std::atan2 (C, S);
	real_type temp2 = 2 + 4 * man_in_data.Ksq - 2 * man_in_data.dcos -
										4 * man_in_data.sKmax * man_in_data.dsin;
	if (temp2 >= 0)
	{
		real_type temp3 = man_in_data.invK * std::sqrt (temp2);
		real_type t1		= man_in_data.invK * mod2pi (man_in_data.th0 - temp1);
		real_type t2		= temp3;
		real_type t3		= man_in_data.invK * mod2pi (temp1 - man_in_data.th1);
		real_type lc		= t1 + t2 + t3;
		//    std::cout << std::setprecision(12) << "RSR len: " << lc << std::endl;
		if (lc < man_out_data.len)
		{
			man_out_data.len	= lc;
			man_out_data.ss1	= t1;
			man_out_data.ss2	= t2;
			man_out_data.ss3	= t3;
			man_out_data.sk1	= -1;
			man_out_data.sk2	= 0;
			man_out_data.sk3	= -1;
			man_out_data.type = D_TYPE::RSR;
		}
		return true;
	}
	return false;
}

/**
 * @brief Compute the maneuver given the type.
 * @param type The type of maneuver.
 * @param data_in The input data.
 * @param data_out The output data.
 * @return True if the maneuver could be computed.
 */
LEN_T
Dubins::comp_man (D_TYPE type)
{
	real_type lambda;
	real_type sKmax;
	Angle phi, sth0, sth1;
	this->scaleToStandard (phi, lambda, sth0, sth1, sKmax);

	real_type invK	= real_type (1) / sKmax;
	real_type sin_0 = sin (sth0);
	real_type cos_0 = cos (sth0);
	real_type sin_1 = sin (sth1);
	real_type cos_1 = cos (sth1);

	real_type Ksq		= sKmax * sKmax;
	real_type dcos	= cos (sth0 - sth1);
	real_type dcos2 = cos_0 - cos_1;
	real_type dsin	= sin_0 - sin_1;
	real_type scos	= cos_0 + cos_1;
	real_type ssin	= sin_0 + sin_1;

	real_type dth = sth0 - sth1;

	struct Man_input_data man_input_data = {
			.th0	 = sth0,
			.th1	 = sth1,
			.sKmax = sKmax,
			.invK	 = invK,
			.cos_0 = cos_0,
			.cos_1 = cos_1,
			.Ksq	 = Ksq,
			.dcos	 = dcos,
			.dcos2 = dcos2,
			.dsin	 = dsin,
			.scos	 = scos,
			.ssin	 = ssin,
			.dth	 = dth};

	struct Man_output_data man_output_data = {
			.len	= std::numeric_limits<real_type>::max(),
			.ss1	= 0.0,
			.ss2	= 0.0,
			.ss3	= 0.0,
			.sk1	= 0.0,
			.sk2	= 0.0,
			.sk3	= 0.0,
			.type = D_TYPE::INVALID};

	switch (type)
	{
		case D_TYPE::LRL: {
			//      std::cout << "Choosing LRL" << std::endl;
			this->comp_LRL (man_input_data, man_output_data);
			break;
		}
		case D_TYPE::RLR: {
			//      std::cout << "Choosing RLR" << std::endl;
			this->comp_RLR (man_input_data, man_output_data);
			break;
		}
		case D_TYPE::LSL: {
			//      std::cout << "Choosing LSL" << std::endl;
			this->comp_LSL (man_input_data, man_output_data);
			break;
		}
		case D_TYPE::LSR: {
			//      std::cout << "Choosing LSR" << std::endl;
			this->comp_LSR (man_input_data, man_output_data);
			break;
		}
		case D_TYPE::RSL: {
			//      std::cout << "Choosing RSL" << std::endl;
			this->comp_RSL (man_input_data, man_output_data);
			break;
		}
		case D_TYPE::RSR: {
			//      std::cout << "Choosing RSR" << std::endl;
			this->comp_RSR (man_input_data, man_output_data);
			break;
		}
		default: {
			throw std::runtime_error ("Invalid maneuver type.");
		}
	}

	if (man_output_data.type == D_TYPE::INVALID)
	{
		throw std::runtime_error ("No valid Dubins path found.");
	}

	// ScaleFromStandard
	this->dtype (man_output_data.type);
	this->s1 (man_output_data.ss1 * lambda);
	this->k1 (man_output_data.sk1 * this->kmax());
	this->s2 (man_output_data.ss2 * lambda);
	this->k2 (man_output_data.sk2 * this->kmax());
	this->s3 (man_output_data.ss3 * lambda);
	this->k3 (man_output_data.sk3 * this->kmax());

	return this->l();
}

/*!
 * Given the standardized version, compute the best word. Credit to Marco Frego & Paolo
 * Bevilacqua.
 * @param th0 The initial standardized angle.
 * @param th1 The final standardized angle.
 * @param lambda A multiplier.
 * @param sKmax The standardized curvature.
 */
void
Dubins::computeBest (Angle th0, Angle th1, real_type lambda, K_T& sKmax)
{
	real_type invK	= real_type (1) / sKmax;
	real_type sin_0 = sin (th0);
	real_type cos_0 = cos (th0);
	real_type sin_1 = sin (th1);
	real_type cos_1 = cos (th1);

	real_type Ksq		= sKmax * sKmax;
	real_type dcos	= cos (th0 - th1);
	real_type dcos2 = cos_0 - cos_1;
	real_type dsin	= sin_0 - sin_1;
	real_type scos	= cos_0 + cos_1;
	real_type ssin	= sin_0 + sin_1;

	real_type dth = th0 - th1;

	struct Man_input_data man_input_data = {
			.th0	 = th0,
			.th1	 = th1,
			.sKmax = sKmax,
			.invK	 = invK,
			.cos_0 = cos_0,
			.cos_1 = cos_1,
			.Ksq	 = Ksq,
			.dcos	 = dcos,
			.dcos2 = dcos2,
			.dsin	 = dsin,
			.scos	 = scos,
			.ssin	 = ssin,
			.dth	 = dth};

	//  std::cout << "Maneuver Input Data:" << std::endl
	//            << "th0: " << th0 << std::endl
	//            << "th1: " << th1 << std::endl
	//            << "sKmax: " << sKmax << std::endl
	//            << "invK: " << invK << std::endl
	//            << "Ksq: " << Ksq << std::endl
	//            << "dcos: " << dcos << std::endl
	//            << "dcos2: " << dcos2 << std::endl
	//            << "dsin: " << dsin << std::endl
	//            << "scos: " << scos << std::endl
	//            << "ssin: " << ssin << std::endl
	//            << "dth: " << dth << std::endl;

	struct Man_output_data man_output_data = {
			.len	= std::numeric_limits<real_type>::max(),
			.ss1	= 0.0,
			.ss2	= 0.0,
			.ss3	= 0.0,
			.sk1	= 0.0,
			.sk2	= 0.0,
			.sk3	= 0.0,
			.type = D_TYPE::INVALID};

	if (!this->comp_LRL (man_input_data, man_output_data))
	{
		//    std::cout << "Cannot find valid LRL" << std::endl;
	}
	if (!this->comp_RLR (man_input_data, man_output_data))
	{
		//    std::cout << "Cannot find valid RLR" << std::endl;
	}
	if (!this->comp_LSL (man_input_data, man_output_data))
	{
		//    std::cout << "Cannot find valid LSL" << std::endl;
	}
	if (!this->comp_LSR (man_input_data, man_output_data))
	{
		//    std::cout << "Cannot find valid LSR" << std::endl;
	}
	if (!this->comp_RSL (man_input_data, man_output_data))
	{
		//    std::cout << "Cannot find valid RSL" << std::endl;
	}
	if (!this->comp_RSR (man_input_data, man_output_data))
	{
		//    std::cout << "Cannot find valid RSR" << std::endl;
	}

	if (man_output_data.type == D_TYPE::INVALID)
	{
		//    throw std::runtime_error("No valid Dubins path found.");
	}

	// ScaleFromStandard
	this->dtype (man_output_data.type);
	this->s1 (man_output_data.ss1 * lambda);
	this->k1 (man_output_data.sk1 * this->kmax());
	this->s2 (man_output_data.ss2 * lambda);
	this->k2 (man_output_data.sk2 * this->kmax());
	this->s3 (man_output_data.ss3 * lambda);
	this->k3 (man_output_data.sk3 * this->kmax());
}

std::vector<std::vector<double>>
Dubins::split_wise()
{
	std::vector<std::vector<double>> res;
	uint n_split = 5;
	std::vector<double> tmp;

	Configuration2 curr = *this->ci();
	for (uint seg = 1; seg < 4; seg++)
	{
		// If the segment has length 0, skip it
		if (this->L (seg) == 0) continue;

		// If the segment is a straight line, add only the beginning and the end
		Configuration2 next_goal_c = circleLine (this->L (seg), this->k (seg), curr);
		if (this->k (seg) == 0 && next_goal_c.th() == curr.th())
		{
			tmp = {curr.x(), curr.y(), curr.th(), this->k (seg), this->L (seg)};
			res.push_back (tmp);
			curr = circleLine (this->L (seg), this->k (seg), curr);
		}
		// If the segment is a circle, add n_split points
		else
		{
			// Dynamically change the value of n_split so that the points are at least 0.01
			// apart
			n_split = std::min ((uint)3, (uint)std::ceil (this->L (seg) / 0.01));

			for (uint i = 0; i < n_split; i++)
			{
				tmp = {curr.x(), curr.y(), curr.th(), this->k (seg), (this->L (seg) / n_split)};
				res.push_back (tmp);
				curr = circleLine (this->L (seg) / n_split, this->k (seg), curr);
			}
		}
	}

	return res;
}

#ifdef MPDP_DRAW
void
Dubins::draw (
		std::ofstream& file,
		std::string label,
		size_t width,
		size_t height,
		bool solve,
		bool close,
		bool init)
{
	if (solve) { this->solve(); }

	if (init) { initAsyFile (file); }

	// Initial point
	Configuration2 c = this->ci()[0];
	file << "p = clothoidPoints((" << c.x() << "," << c.y() << "), " << c.th() << ","
			 << this->k1() << ", 0, " << this->s1() << ");" << std::endl;
	file << "draw(p,royalblue);" << std::endl;
	file << "dot((" << c.x() << "," << c.y() << "), red);" << std::endl;
	if (label != "")
	{
		file << "label(\"$" << label << "$\", (" << c.x() << ", " << c.y() << "));"
				 << std::endl;
	}

	// Intermediate point
	c = circleLine (this->s1(), this->k1(), c);
	file << "p = clothoidPoints((" << c.x() << "," << c.y() << "), " << c.th() << ","
			 << this->k2() << ", 0, " << this->s2() << ");" << std::endl;
	file << "draw(p,royalblue);" << std::endl;
	file << "dot((" << c.x() << "," << c.y() << "), red);" << std::endl;

	// Final point
	c = circleLine (this->s2(), this->k2(), c);
	file << "p = clothoidPoints((" << c.x() << "," << c.y() << "), " << c.th() << ","
			 << this->k3() << ", 0, " << this->s3() << ");" << std::endl;
	file << "draw(p,royalblue);" << std::endl;
	file << "dot((" << c.x() << "," << c.y() << "), red);" << std::endl;

	// Plot points
	file << "dot((" << this->ci()->x() << "," << this->ci()->y() << "), black);"
			 << std::endl;
	file << "dot((" << this->cf()->x() << "," << this->cf()->y() << "), purple+3bp);"
			 << std::endl;

	if (close) { file.close(); }
}
#endif	// MPDP_DRAW

#endif	// CUDA_ON