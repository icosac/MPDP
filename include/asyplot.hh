/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi and Marco Frego                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |      email: marco.frego@unitn.it                                         |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef ASYPLOT_H
#define ASYPLOT_H

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
 

using std::string;
using std::ofstream;

class AsyPlot {
public:
	AsyPlot( string filename );
	~AsyPlot();

	void
	writeLine(string const & line) const;


	void dot( double x, double y, string const & penna="black" ) const ;


	void
	drawRect( double x0, double y0,
				double x1, double y1,
				double x2, double y2,
				double x3, double y3,
				string const & penna="black") const ;

	void
	fillRect(double x0, double y0,
	  			double x1, double y1,
				double x2, double y2,
				double x3, double y3,
				string const & pennaBordo = "black",
				string const & pennaRiemp = "white" ) const;

	void
	drawLine( double x0, double y0,
				double x1, double y1,
				std::string const & penna="black" ) const ;

	void
	drawAngle(double s0, double th0,
		double s1, double th1, double dk,
		string const & penna, int npts = 100) const;


	void
	label( string const & text,
			double      x,
			double      y,
			string const & placement = "",
			string const & penna = "black" ) const ;

	void
	draw(string const & text,
			string const & penna) const;

	void
	displayAxes( string const & labX,
					string const & labY,
					double      xmin,
					double      xmax,
					double      ymin,
					double      ymax ) const ;

	void
	xyLimits( double      xmin,
				double      xmax,
				double      ymin,
				double      ymax,
				bool           crop) const;
private:
	mutable ofstream file;
	string  filename;
	bool showAxes;
	bool openFile();
	bool closeFile();
	void initFile();
	void displayAxes() const ;
	void compileFile();
};


#endif
