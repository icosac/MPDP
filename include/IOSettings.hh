#ifndef IOSETTINGS_HH
#define IOSETTINGS_HH

namespace SET{
	//Input
	const bool 	generator = false; //Remember to set nPoints
	const char* filename  = "file.txt";
				char* outfile		= ""; //.asy output file

	//Curve parameters
	const size_t nPoints 	 = 10; //This needs to be set only if SET::generator is true
	const K_T 	 Kmax 	 	 = 0.5;
	const size_t discr   	 = 4;
	const size_t nRef    	 = 4;
	const bool	 saveAngles  = true;
	const short  CUDAFunType = 2;

	//Generator: parameters used if SET::generator is true
	const size_t xSize = 10;
	const size_t ySize = 10;

	//Input by file: SET::filename must be a valid path to an input file
	//Please manually change the function name in the source file to choose between taking points or Configuration2 in input.
	const int 	type = 0;
	const Angle thi  = ANGLE::FREE; //Mind to change these values if taking points in input.
	const Angle thf	 = ANGLE::FREE;

};

#endif //IOSETTINGS_HH