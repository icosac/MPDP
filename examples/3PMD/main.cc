/**
 * @file main.cpp
 * @author Enrico Saccon <enricosaccon96@gmail.com>
 * @license This project is released under the GNU Public License 3.0.
 * @copyright Copyright 2020 Enrico Saccon. All rights reserved.
 * @brief Main file for the Dubins and Reed-Shepp paths computation.
 */

#include <3PMD.hh>

int main(int argc, char** argv){
	// main3PMDBruteForce();
	// main3PMDBruteForceWithPlot();
	// main3PDP();

	// compute3Pman();

	generateDataset3PDPCircle(argc, argv);
	// generateDataset3PDPRect(argc, argv);
	return 0;
}