/**
 * @file main.cu
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Main file for the Dubins and Reed-Shepp paths computation using CUDA.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>

#include <utils.cuh>
#include <dubins.cuh>
#include <dp.cuh>
#include <timeperf.hh>
#include <utilities.hh>

#include <tests.hh>

std::vector<Configuration2> example1 = {
		Configuration2 (0, 0, -2.0 * M_PI / 8.0), Configuration2 (2, 2, ANGLE::FREE),
		Configuration2 (6, -1, ANGLE::FREE), Configuration2 (8, 1, 2.0 * M_PI / 8.0)};

std::vector<std::string> testsNames = {"Kaya Example 1", "Kaya Example 2",
																			 "Kaya Example 3", "Kaya Example 4",
																			 "Omega",					 "Circuit"};

std::vector<std::vector<Configuration2>> Tests = {kaya1, kaya2, kaya3, kaya4, omega, spa};

std::vector<K_T> Ks								= {3.0, 3.0, 5.0, 3.0, 3.0, 3.0};
std::vector<uint> discrs					= {4, 16, 90, 360};
std::vector<uint> refins					= {1, 2, 4, 8, 16};
std::vector<LEN_T> exampleLenghts = {
		3.41557885807514871601142658619,
		6.27803455030931356617429628386,
		11.9162126542854860389297755319,
		7.46756219733842652175326293218,
		41.0725016438839318766440555919,
		6988.66098639942993031581863761};	 // the last length is SPA

std::string
nameTest (std::string name, std::string add = "", std::string conc = " ")
{
	if (add == "") { return name; }
	else { return name + conc + add; }
}

__global__ void
dubinsL (Configuration2 c0, Configuration2 c1, real_type k, real_type* L)
{
	Dubins dubins (c0, c1, k);
	L[0] += dubins.l();
}

int
main (int argc, char* argv[])
{
	cudaFree (0);

	int devicesCount;
	cudaGetDeviceCount (&devicesCount);
	cudaDeviceProp deviceProperties;
	cudaGetDeviceProperties (&deviceProperties, 0);

	std::cout << "Running CUDA" << std::endl;

	if (argc == 1)
	{
		for (int testID = 0; testID < 6; testID++)
		{
			// if (testID!=3){continue;}
			real_type dLen = exampleLenghts[testID];

			std::vector<bool> fixedAngles;
			for (uint i = 0; i < Tests[testID].size(); i++)
			{
				if (i == 0 || i == Tests[testID].size() - 1) { fixedAngles.push_back (true); }
				else { fixedAngles.push_back (false); }
			}
			std::vector<real_type> curveParamV = {Ks[testID], 3};
			real_type* curveParam							 = curveParamV.data();

			for (auto DISCR : discrs)
			{
				if (DISCR != 360) { continue; }
				for (auto r : refins)
				{
					// if (r!=16){continue;}
					TimePerf tp, tp1;
					std::vector<Configuration2> points = Tests[testID];

					tp.start();
					LEN_T Length =
							DP::solveDP (points, fixedAngles, curveParamV, DISCR, r, true, 2).first;
					auto time1 = tp.getTime();

					LEN_T* Length1;
					cudaMallocManaged (&Length1, sizeof (LEN_T));
					for (unsigned int idjijij = points.size() - 1; idjijij > 0; idjijij--)
					{
						dubinsL<<<1, 1>>> (points[idjijij - 1], points[idjijij], Ks[testID], Length1);
						cudaDeviceSynchronize();

						Dubins c (points[idjijij - 1], points[idjijij], Ks[testID]);
						Length += c.l();
					}

					printf ("%3d & %2d & ", DISCR, r);
					PrintScientific2D ((Length - exampleLenghts[testID]) * 1000.0);
					// printf(" & ");
					// PrintScientific2D((Length1[0]-Length)*1000.0);
					// printf(" & ");
					// PrintScientific2D((Length1[0]-exampleLenghts[testID])*1000.0);
					printf (" & ");
					PrintScientific1D (time1);
					// printf("&%.16f", Length);
					// printf("&%.16f\\\\\n", Length1[0]);
					printf ("\\\\\n");

					cudaFree (Length1);
				}
			}
			printf ("\n\n\n\n");
		}
	}
	return 0;
}
