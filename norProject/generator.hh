#ifndef GENERATOR_HH
#define GENERATOR_HH

#include<cmath>
#include<vector>
#include<random>
#include<ctime>

#ifndef CUDA_ON 
#include<configuration.hh>
#include<utils.hh>
#else
#include<configuration.cuh>
#include<utils.cuh>
#endif


std::vector<Configuration2> getPoints(size_t size, size_t xSize=100, size_t ySize=100,
																			std::vector<bool>* fixedAngles=nullptr){
	std::vector<Configuration2> points (size);

	//Create random elements.
	std::uniform_real_distribution<double> unifx(0, xSize);
	std::uniform_real_distribution<double> unify(0, ySize);
	std::uniform_real_distribution<double> unifth(0, m_2pi);
	std::default_random_engine re(time(NULL));

	//Set values.
	for(uint i=0; i<size; i++){
		points[i].x(unifx(re));
		points[i].y(unify(re));
		if(fixedAngles!=nullptr && fixedAngles[0][i]){
			points[i].th(unifth(re));
		}
		else{
			points[i].th(ANGLE::FREE);
		}
	}

	for(auto point : points){
		std::cout << std::setprecision(20) << point << std::endl;
	}

	return points;
}

#endif //GENERATOR_HH
