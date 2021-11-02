#ifndef IOUTILS_HH
#define IOUTILS_HH

#ifndef CUDA_ON 
#include<configuration.hh>
#else
#include<configuration.cuh>
#endif

#include<iostream>
#include<fstream>
#include<vector>
#include<limits>
#include<cstdlib>

#include<asyplot.hh>

void PrintScientific1D(real_type d){
  if (d == 0)
  {
    printf ("%*d", 6, 0);
    return;
  }

  int exponent  = (int)floor(log10( fabs(d)));  // This will round down the exponent
  real_type base   = d * pow(10, -1.0*exponent);

  printf("%1.1lfe%+01d", base, exponent);
}

void PrintScientific2D(real_type d){
  if (d == 0)
  {
    printf ("%*d", 7, 0);
    return;
  }

  int exponent  = (int)floor(log10( fabs(d)));  // This will round down the exponent
  real_type base   = d * pow(10, -1.0*exponent);

  printf("%1.1lfe%+02d", base, exponent);
}

std::vector<std::string> splitString(std::string str){
	size_t pos1=0;
	std::vector<std::string> ret;
	while(pos1!=std::string::npos){
		pos1=str.find(" ");
		//std::cout << "[" << str << "] " << pos1 << " " << std::string::npos << std::endl;
		if(pos1!=std::string::npos){
			ret.emplace_back(str.substr(0, pos1));
			str=str.substr(pos1+1, str.length());
		}
		else{
			std::string str1=str.substr(0, str.length());
			if(str1.length()>0){
				ret.emplace_back(str1);
			}
		}
	} 
	return ret;
}

/*!
 * Function to read from a file the `Configuration2`.
 * @param inputFile A `fstream` to the input file.
 * @param type The type of file: 0 or 1. Default to 0.
 * @param close Whether to close the file after reading or not. Default to true.
 * @return Returns a pair where the first element is vector of `Configuration2` and the second element is a vector of booleans for the fixed angles.
 */
std::pair<std::vector<Configuration2>, std::vector<bool> > 
readConfigurationsFromFile(std::fstream& inputFile, int type=0, bool close=true){
	std::vector<Configuration2> points;
	std::vector<bool> fixedAngles;

	if (type==0){
		real_type x, y;
		unsigned short fixed; 
		std::string th, appString;

		while(getline(inputFile, appString)){
			int i=0; 
			for(auto el : splitString(appString)){
				switch(i){
					case 0: x=atof(el.c_str()); break;
					case 1: y=atof(el.c_str()); break;
					case 2: {
						th=el;
						if (th.find("FREE")!=std::string::npos){ 
							points.emplace_back(Configuration2(x, y, ANGLE::FREE));
							fixedAngles.emplace_back(false);
						}
						else{
							points.emplace_back(Configuration2(x, y, atof(th.c_str())));
						}
						break;
					}
					case 3: {
						if(el=="0")	{ fixedAngles.emplace_back(false); }
						else 				{ fixedAngles.emplace_back(true); }		
					}
				}
				i++;
			}			
		}
	}
	if (type==1){
		uint counter=0;
		size_t dim; 
		real_type value;
		std::string th;
		int fixed;
		
		inputFile >> dim;
		points.resize(dim);

		while (inputFile >> value){
			if((int)(counter/dim)==0){
				points[counter%dim].x(value);
			}
			else if((int)(counter/dim)==1){
				points[counter%dim].y(value);
			}
			else {
				points[counter%dim].th(value);
			}
			counter++;
			if ((int)(counter/dim)==2){ break; }
		}
		counter=0;
		while (inputFile >> th){
			if (th.find("FREE")!=std::string::npos){ 
				points[counter].th(ANGLE::FREE);
			}
			else{
				points[counter].th(atof(th.c_str()));
			}
			counter++;
			if ((int)(counter/dim)==1){ break; }
		}
		while (inputFile >> fixed){
			if(fixed==0){ fixedAngles.emplace_back(false); }
			else 				{ fixedAngles.emplace_back(true); }
		}
	}

	if (close) { inputFile.close(); }

	return std::pair<std::vector<Configuration2>, std::vector<bool> >(points, fixedAngles);
}

/*!
 * Function to read from a file the `Configuration2`.
 * @param inputFile A string containing the path to the input file.
 * @param type The type of file: 0 or 1. Default to 0.
 * @return Returns a pair where the first element is vector of `Configuration2` and the second element is a vector of booleans for the fixed angles.
 */
std::pair<std::vector<Configuration2>, std::vector<bool> >
readConfigurationsFromFile(const char* filename, int type=0){
	std::fstream inputFile;
	inputFile.open(filename, std::fstream::in);
	return readConfigurationsFromFile(inputFile,  type, true);
}

/*!
 * Function to read from a file the points for the `Configuration2`. The first and final angle are manually set
 * @param inputFile A `fstream` to the input file.
 * @param thi The initial angle to be set. Default is ANGLE::FREE.
 * @param thf The initial angle to be set. Default is ANGLE::FREE.
 * @param type The type of file: 0 or 1. Default to 0.
 * @param close Whether to close the file after reading or not. Default to true.
 * @return Returns a pair where the first element is vector of `Configuration2` and the second element is a vector of booleans for the fixed angles.
 */
std::pair<std::vector<Configuration2>, std::vector<bool> >
readPointsFromFile(std::fstream& inputFile, Angle thi=ANGLE::FREE, Angle thf=ANGLE::FREE, int type=0, bool close=true){
	std::vector<Configuration2> points;
	
	if (type==0){
		real_type x, y;
		while(inputFile >> x >> y){
			points.emplace_back(Configuration2(x, y, ANGLE::FREE));
		}	
	}
	
	if (type==1){
		uint counter=0;
		size_t dim; 
		real_type value;
		
		inputFile >> dim;
		points.resize(dim);

		while (inputFile >> value){
			if(counter<dim){
				points[counter%dim].x(value);
			}
			else {
				points[counter%dim].y(value);
			}
			counter++;
		}
	}

	points.front().th(thi);
	points.back().th(thf);
	
	if (close) { inputFile.close(); }
	
	std::vector<bool> fixedAngles(points.size(), false); fixedAngles.front()=true; fixedAngles.back()=true;

	return std::pair<std::vector<Configuration2>, std::vector<bool> >(points, fixedAngles);
}

/*!
 * Function to read from a file the points for the `Configuration2`. The first and final angle are manually set
 * @param inputFile A string containing the path to the input file.
 * @param thi The initial angle to be set. Default is ANGLE::FREE.
 * @param thf The initial angle to be set. Default is ANGLE::FREE.
 * @param type The type of file: 0 or 1. Default to 0.
 * @return Returns a pair where the first element is vector of `Configuration2` and the second element is a vector of booleans for the fixed angles.
 */
std::pair<std::vector<Configuration2>, std::vector<bool> >
readPointsFromFile(const char* filename, Angle thi=ANGLE::FREE, Angle thf=ANGLE::FREE, int type=0){
	std::fstream inputFile;
	inputFile.open(filename, std::fstream::in);
  return readPointsFromFile(inputFile, thi, thf, type, true);
}

/*!
 * This function returns informations regarding the Dubins composing the path.
 * @param points A vector of points with or without the newly computed angles.
 * @param kmax The maximum curvature used to compute the best angles.
 * @param vtheta The best angles. Default is nullptr, instead if set, the function will assign the angles to the point before starting.
 * @param len The length the MPMD should have. Deafult is 0.0, instead if set, it's used to check whether the now computed Dubins are correct.
 */
void getMPMDInfo(std::vector<Configuration2> points, K_T kmax, const std::vector<Angle>* vtheta=NULL, LEN_T len=0.0){
	std::vector<Dubins> dubinss;
	LEN_T Len=0.0;
	for (uint i=0; i<points.size()-1; i++){
		if (vtheta!=NULL){
			if (i==0) { points[i].th((*vtheta)[i]); }
			points[i+1].th((*vtheta)[i+1]);
		}
		dubinss.emplace_back(Dubins(points[i], points[i+1], kmax));
		Len+=dubinss.back().l();
	}
	if (len!=0 && !eq<double>(len, Len, 1e-13)){
		std::cout << "ERROR before: " << std::setprecision(20) << len << " after: " << Len << " difference: ";
    PrintScientific2D(std::abs(len-Len));
    std::cout << std::endl;
	}
	else{
    std::cout << "Total length: " << Len << std::endl;
	}
	for (Dubins d : dubinss){
		std::cout << d << std::endl;
	}
}

/*!
 * Function to compute the sinc of a value.
 * @param x The value.
 * @return sinc(x)
 */
static inline double sinc(double x){
	if (std::abs(x) < 0.002) {
		double xs = x * x;
		return 1 - xs / 6. * (1 - xs / 20.0);
	}
	else
	{
		return std::sin(x) / x;
	}
}

/*!
 * Credit to Marco Frego and Enrico Bertolazzi.
 */
static inline double f(double ell, double k, double th){
	double tmp = k * ell*0.5;
	return ell * sinc(tmp)*std::cos(th + tmp);
}

/*!
 * Credit to Marco Frego and Enrico Bertolazzi.
 */
static inline double g(double ell, double k, double th){
	double tmp = k * ell*0.5;
	return ell * sinc(tmp)*std::sin(th + tmp);
}

/*!
 * Credit to Marco Frego and Enrico Bertolazzi.
 */
void draw(Dubins & dc, AsyPlot& plot, string const & penna)
{
  std::stringstream str;
  
  double x0 = dc.ci()->x();
  double y0 = dc.ci()->y();
  double th0 = dc.ci()->th();
  double k0 = dc.k1();
  double s0 = dc.s1();

  double xM0 = f(s0, k0, th0) + x0;
  double yM0 = g(s0, k0, th0) + y0;
  double thetaM0 = th0 + k0*s0;
  double k1 = dc.k2();
  double s1 = dc.s2();

  double xM1 = f(s1, k1, thetaM0) + xM0;
  double yM1 = g(s1, k1, thetaM0) + yM0;
  double thetaM1 = thetaM0 + k1*s1;
  double k2 = dc.k3();
  double s2 = dc.s3();

  str		<< "path pclot = clothoidPoints(("
      << x0 << ','
      << y0 << "),"
      << th0 << ','
      << k0 << ','
      << 0 << ','
      << s0 << ','
      << "50,0);\n"
      << "pen penna = " << penna << ";\n"
      << "draw(pclot, penna);\n\n";

  str << "path pclot = clothoidPoints(("
    << xM0 << ','
    << yM0 << "),"
    << thetaM0 << ','
    << k1 << ','
    << 0 << ','
    << s1 << ','
    << "50,0);\n"
    << "pen penna = " << penna << ";\n"
    << "draw(pclot, penna);\n\n";

  str << "path pclot = clothoidPoints(("
    << xM1 << ','
    << yM1 << "),"
    << thetaM1 << ','
    << k2 << ','
    << 0 << ','
    << s2 << ','
    << "50,0);\n"
    << "pen penna = " << penna << ";\n"
    << "draw(pclot, penna);\n\n";

  // str << "dot((" << x0 << ','
  //   << y0 << "), red+4bp);\n\n";

  // str << "dot((" << x1 << ','
  //   << y1 << "), red+4bp);\n\n";

  // str << "dot((" << xM0 << ','
  //   << yM0 << "), royalblue+2bp);\n\n";
  // str << "dot((" << xM1 << ','
  //   << yM1 << "), royalblue+2bp);\n\n";

  plot.writeLine(str.str());
  
}

void drawSolution(std::vector<Configuration2> points, double Kmax, const char* filename, const std::vector<Angle>* vtheta=NULL){
  AsyPlot ap(filename);
  for (int i=1; i<points.size(); ++i) {
  	if (vtheta!=NULL){
  		if (i==1) { points[0].th((*vtheta)[0]); }
  		points[i].th((*vtheta)[i]);
  	}
    Dubins d(points[i-1], points[i], Kmax);
    draw(d, ap, "blue");
  }

  for (int i=0; i<points.size(); ++i) {
    ap.dot(points[i].x(), points[i].y(), "red+2bp");
  }
}



#endif //IOUTILS_HH