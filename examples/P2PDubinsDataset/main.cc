#include "P2PDubinsDataset.hh"

/**
 * @brief Main to execute either function to generate a random dataset of P2P Dubins 
 * (genDSDubinsP2P) or to test the Dubins procedure against a dataset (testDubins). 
 */
int main(int argc, char** argv){
  return genDSDubinsP2P(true);
  // return testDubins();
}