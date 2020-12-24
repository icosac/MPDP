#define DEBUG
//#define CLOTHOID CLOTHOID //Need to rewrite CLOTHOID because of typedefs.hh
#define DUBINS DUBINS

#if defined(CLOTHOID)
#include <clothoidG1.hh>
#elif defined(DUBINS)
#include <dubins.hh>
#endif
#include <dp.hh>
#include <timeperf.hh>

#include<iostream>
#include<cmath>
#include<random>
#include<vector>
#include<utility>
#include<iomanip>
#include<algorithm>

using namespace std;

#define SIZE 349

//double X[SIZE] = {2.9265642,2.6734362,2.5109322,1.9078122,1.1859282,1.9249962,2.8265562,0.00468420000000025,-2.826567,-1.9437558,-1.1859438,-1.9062558,-2.501565,-2.6734386,-2.9265642,-2.6187522,-1.1406318,-0.8968758,-1.4562558,-1.9062558,-0.00468780000000013,1.9078122,1.4468682,0.8968722,1.1406282,2.6187522, 2.9265642 };
//double Y[SIZE] = {-1.707808758,-1.707808758,-2.367185958,-2.582810358,-2.582810358,-1.167184758,0.915619242,3.178123242,0.915619242,-1.150000758,-2.582810358,-2.582810358,-2.393750358,-1.707808758,-1.707808758,-3.178123242,-3.178123242,-2.989063158,-0.915616758,0.925003242,2.953123242,0.925003242,-0.915616758,-2.989063158,-3.178123242,-3.178123242, -1.707808758 };

double XY[SIZE*2]= {-475.906, 116.543, -397.488, 193.652, -301.255, 291.346, -246.139, 346.406, -224.85,  367.324, -217.014, 374.77, -213.437, 378.085, -209.305, 380.619, -204.295, 382.007, -200.013, 381.8, -196.202, 380.19, -188.959, 376.596, -181.489, 372.948, -174.962, 369.955, -168.654, 366.995, -159.089, 363.543, -143.082, 359.952, -127.159, 358.404, -114.102, 359.089, -100.154, 360.808, -83.4549, 365.093, -71.9814, 370.041, -62.4792, 374.727, -53.408,  380.282, -45.9236, 386.642, -39.8024, 391.743, -35.8917, 395.707, -31.0709, 401.026, -18.9996, 414.249, -11.0071, 423.137, -4.65741, 430.081, 5.4544, 441.119, 13.4286,  449.776, 23.675, 460.189, 35.4117,  472.497, 43.3489,  480.37, 50.0755,  486.952, 60.1424,  496.321, 67.2874,  503.059, 75.4194,  510.653, 84.4191,  518.893, 92.9769,  526.594, 103.27, 535.662, 111.91, 543.114, 120.782,  550.421, 131.729,  559.417, 140.771,  566.794, 152.423,  576.031, 189.823,  605.257, 210.436,  620.614, 241.594,  640.625, 268.355,  657.42, 297.003,  673.761, 349.202,  702.094, 368.643,  712.255, 380.76, 718.774, 386.161,  720.697, 390.891,  721.054, 395.532,  720.525, 398.992,  718.894, 402.457,  716.573, 405.428,  713.454, 407.407,  710.222, 408.593,  706.687, 408.792,  706.688, 409.53, 702.687, 411.625,  692.004, 412.695,  685.282, 414.098,  677.737, 416.595,  663.032, 419.724,  645.422, 424.11, 621.841, 426.416,  610.201, 428.212,  602.925, 431.363,  595.841, 434.056,  591.592, 437.822,  587.549, 441.597,  584.408, 446.632,  581.579, 452.387,  579.518, 459.109,  578.595, 465.818,  578.655, 476.779,  579.029, 488.601,  579.782, 531.443,  581.181, 564.525,  582.331, 624.118,  584.569, 637.453,  584.82, 644.987,  584.275, 651.711,  583.082, 657.835,  580.863, 662.436,  578.27, 668.032,  573.915, 671.916,  570.344, 676.496,  565.281, 682.583,  557.15, 689.683,  547.306, 710.058,  519.191, 731.613,  489.688, 745.257,  470.631, 752.274,  460.452, 755.795,  455.018, 760.231,  447.534, 764.555,  439.671, 768.628,  432.119, 774.028,  420.742, 778.271,  410.985, 786.587,  388.731, 798.397,  355.42, 805.44, 336.231, 811.206,  321.426, 817.647,  306.695, 824.544,  291.801, 833.031,  274.954, 839.854,  263.3, 852.015,  243.844, 860.742,  229.78, 863.622,  223.592, 865.574,  216.627, 864.742,  210.838, 861.594,  204.657, 856.948,  198.378, 853.849,  194.904, 848.858,  189.947, 842.145,  182.575, 835.807,  175.213, 830.261,  167.932, 825.118,  160.217, 820.456,  152.119, 816.058,  141.91, 812.201,  130.916, 809.731,  120.391, 808.048,  110.759, 807.012,  101.358, 806.406,  93.098, 805.731,  83.4651, 805.05, 75.1089, 803.896,  65.1221, 803.173,  58.9104, 802.12, 53.1602, 800.613,  45.9217, 798.69, 38.5464, 796.473,  31.3983, 794.229,  24.8616, 791.476,  17.7252, 786.846,  7.3409, 781.463,  -3.14673, 773.309,  -16.831, 762.915,  -30.7363, 752.168,  -43.2059, 743.583,  -52.0516, 731.71, -61.998, 719.5,  -70.7443, 706.676,  -78.808, 695.355,  -84.8304, 681.35, -90.4866, 673.649,  -93.5025, 664.082,  -96.6336, 654.893,  -98.9343, 647.948,  -100.514, 641.985,  -101.674, 598.849,  -109.354, 541.932,  -119.438, 521.686,  -123.518, 508.983,  -125.414, 493.083,  -128.32, 482.794,  -130.1, 469.689,  -132.187, 452.535,  -135.03, 441.273,  -137.827, 426.795,  -140.984, 413.243,  -143.962, 406.695,  -145.078, 402.613,  -145.426, 397.199,  -144.56, 392.829,  -141.774, 389.174,  -137.169, 387.386,  -131.565, 386.226,  -125.067, 384.267,  -116.668, 381.553,  -108.947, 379.333,  -102.711, 376.297,  -95.7949, 372.691,  -88.7854, 367.29, -81.7638, 361.563,  -75.4276, 353.267,  -68.665, 346.631,  -64.8707, 336.896,  -61.1643, 328.056,  -57.7864, 317.591,  -55.125, 299.526,  -51.5055, 248.72, -42.3966, 198.651,  -32.946, 178.779,  -29.3033, 165.762,  -27.2629, 154.386,  -25.7433, 145.099,  -24.9739, 130.861,  -23.9411, 118.438,  -23.5809, 108.801,  -23.6276, 93.4582,  -24.0974, 81.7273,  -25.0471, 67.7252,  -26.5797, 53.7597,  -28.5032, 38.2081,  -31.6256, 13.601, -37.6885, -8.39568, -44.6658, -28.9515, -52.5344, -52.8631, -63.8544, -74.8296, -76.0151, -93.5525, -88.0733, -114.187, -103.325, -137.471, -123.587, -159.552, -146.454, -176.848, -166.229, -194.072, -185.896, -209.545, -203.229, -230.651, -227.58, -238.354, -237.777, -244.564, -247.871, -248.841, -256.288, -252.676, -265.697, -255.609, -275.264, -257.617, -284.094, -258.983, -292.728, -260.029, -301.374, -260.022, -312.498, -259.401, -328.297, -258.987, -341.401, -258.767, -350.344, -259.113, -360.154, -259.764, -367.884, -261.406, -378.026, -263.725, -387.421, -265.985, -394.869, -268.773, -402.138, -272.242, -409.944, -276.488, -418.183, -281.408, -425.738, -292.649, -440.021, -304.259, -453.864, -322.914, -475.92, -341.151, -497.417, -363.285, -522.71, -385.12,  -548.083, -401.337, -566.799, -409.86,  -576.907, -417.604, -584.439, -423.995, -590.228, -431.111, -596.782, -440.554, -604.475, -448.782, -610.262, -459.192, -617.404, -473.362, -625.635, -481.058, -629.759, -492.748, -635.25, -504.391, -640.069, -517.472, -644.886, -526.525, -648.198, -533.819, -650.271, -549.915, -655.244, -574.636, -662.299, -606.693, -671.359, -656.283, -685.542, -690.858, -695.223, -719.191, -703.166, -749.379, -711.521, -772.031, -717.207, -783.521, -719.99, -798.887, -723.712, -803.941, -724.099, -809.157, -723.579, -813.861, -721.887, -817.244, -719.819, -820.79,  -717.147, -823.811, -713.722, -826.366, -709.955, -829.081, -704.067, -834.097, -691.843, -843.86,  -667.793, -863.801, -616.461, -876.063, -585.435, -889.437, -551.673, -895.676, -535.754, -900.78,  -522.593, -903.716, -512.996, -905.543, -505.888, -906.2, -499.263, -906.177, -492.609, -905.724, -486.121, -904.339, -478.627, -901.984, -471.497, -899.23,  -465.559, -895.378, -459.031, -891.02,  -453.628, -886.392, -448.606, -881.882, -445.148, -874.257, -440.048, -864.563, -434.441, -850.728, -427.065, -818.262, -409.435, -786.167, -391.852, -774.342, -385.395, -764.007, -379.649, -753.533, -373.669, -743.059, -367.523, -732.805, -361.367, -723.861, -355.801, -720.166, -353.598, -717.217, -351.78, -713.825, -349.348, -709.461, -346.222, -704.921, -342.523, -702.431, -339.341, -701.13,  -336.135, -700.553, -333.151, -700.63,  -330.488, -701.235, -327.91, -702.202, -325.445, -703.248, -323.589, -704.965, -321.509, -707.015, -319.374, -709.313, -317.198, -734.623, -295.13, -751.97,  -280.146, -761.26,  -272.499, -767.273, -267.573, -774.024, -261.212, -777.692, -256.746, -781.674, -251.326, -784.93,  -246.004, -787.651, -240.19, -790.047, -232.723, -791.473, -224.996, -791.66,  -217.165, -790.941, -210.356, -789.128, -202.57, -786.066, -194.68, -781.735, -187.25, -775.521, -179.33, -764.991, -168.615, -741.021, -145.189, -710.565, -114.649, -677.374, -82.0715, -639.464, -44.885, -600.597, -6.61522, -554.566, 38.8255, -475.906, 116.543};

//vector<Configuration2<double> > createPoints(){
//  vector<Configuration2<double> > ret;
//  for (int i=0; i<SIZE; i++){
//    ret.push_back(Configuration2<double>(X[i], Y[i], ANGLE::INVALID));
//  }
//  return ret;
//}

vector<Configuration2<double> > createPoints(){
  vector<Configuration2<double> > ret;
  for (int i=0; i<SIZE*2; i+=2){
    ret.push_back(Configuration2<double>(XY[i], XY[i+1], ANGLE::INVALID));
  }
  return ret;
}

#include<fstream>

#define READ_FROM_FILE()                                                  \
  ifstream input("build_clothoid.txt");                                   \
  ofstream out ("matlab.txt", fstream::out);                              \
    float x0, y0, th0, x1, y1, th1, k, dk, l;                             \
    int i=0;                                                              \
    while (input >> x0 >> y0 >> th0 >> x1 >> y1 >> th1 >> k >> dk >> l){  \
      i++;                                                                \
      Configuration2<float>ci(x0, y0, th0);                               \
      Configuration2<float>cf(x1, y1, th1);                               \
      ClothoidG1<float>c(ci, cf);

#define CLOSE_FILE() } input.close();

vector<Configuration2<double> > kaya1={
        Configuration2<double> (0, 0, -M_PI/3.0),
        Configuration2<double> (-0.1, 0.3, ANGLE::INVALID),
        Configuration2<double> (0.2, 0.8, ANGLE::INVALID),
        Configuration2<double> (1, 1, -M_PI/6.0)
};

vector<Configuration2<double> > kaya2={
        Configuration2<double> (0, 0, -M_PI/3.0),
        Configuration2<double> (-0.1, 0.3, ANGLE::INVALID),
        Configuration2<double> (0.2, 0.8, ANGLE::INVALID),
        Configuration2<double> (1, 1, ANGLE::INVALID),
        Configuration2<double> (0.5, 0.5, ANGLE::INVALID),
        Configuration2<double> (0.5, 0, -M_PI/6.0)
};

vector<Configuration2<double> > kaya4={
       Configuration2<double>(0.5, 1.2, 5*M_PI/6.0),
       Configuration2<double>(0.0, 0.5, ANGLE::INVALID),
       Configuration2<double>(0.5, 0.5, ANGLE::INVALID),
       Configuration2<double>(1.0, 0.5, ANGLE::INVALID),
       Configuration2<double>(1.5, 0.5, ANGLE::INVALID),
       Configuration2<double>(2.0, 0.5, ANGLE::INVALID),
       Configuration2<double>(2.0, 0.0, ANGLE::INVALID),
       Configuration2<double>(1.5, 0.0, ANGLE::INVALID),
       Configuration2<double>(1.0, 0.0, ANGLE::INVALID),
       Configuration2<double>(0.5, 0.0, ANGLE::INVALID),
       Configuration2<double>(0.0, 0.0, ANGLE::INVALID),
       Configuration2<double>(0.0, -0.5, 0)
};

vector<Configuration2<double> > kaya3={
       Configuration2<double>(0.5, 1.2, 5.0*M_PI/6.0),
       Configuration2<double>(0, 0.8, ANGLE::INVALID),
       Configuration2<double>(0, 0.4, ANGLE::INVALID),
       Configuration2<double>(0.1, 0, ANGLE::INVALID),
       Configuration2<double>(0.4, 0.2, ANGLE::INVALID),
       Configuration2<double>(0.5, 0.5, ANGLE::INVALID),
       Configuration2<double>(0.6, 1, ANGLE::INVALID),
       Configuration2<double>(1, 0.8, ANGLE::INVALID),
       Configuration2<double>(1, 0, ANGLE::INVALID),
       Configuration2<double>(1.4, 0.2, ANGLE::INVALID),
       Configuration2<double>(1.2, 1, ANGLE::INVALID),
       Configuration2<double>(1.5, 1.2, ANGLE::INVALID),
       Configuration2<double>(2, 1.5, ANGLE::INVALID),
       Configuration2<double>(1.5, 0.8, ANGLE::INVALID),
       Configuration2<double>(1.5, 0, ANGLE::INVALID),
       Configuration2<double>(1.7, 0.6, ANGLE::INVALID),
       Configuration2<double>(1.9, 1, ANGLE::INVALID),
       Configuration2<double>(2, 0.5, ANGLE::INVALID),
       Configuration2<double>(1.9, 0, ANGLE::INVALID),
       Configuration2<double>(2.5, 0.6, 0),
};

vector<vector<Configuration2<double> > > Tests = {
  kaya1, kaya2, kaya3, kaya4
};

vector<K_T> Ks = {3.0, 3.0, 5.0, 3.0};
vector<uint> discrs = {4, 120, 360, 720, 2000};

#define DISCR 6

int main (){
  cout << "C++" << endl;
#if false
  for (uint discr : discrs){
    cout << "Discr: " << discr << endl;
    for (uint j=0; j<Tests.size(); j++){
      std::vector<bool> fixedAngles;
      vector<Configuration2<double> > v=Tests[j];
      for (int i=0; i<v.size(); i++){
        if (i==0 || i==v.size()-1) {
          fixedAngles.push_back(true);
        }
        else {
          fixedAngles.push_back(false);
        }
      }
      std::vector<real_type> curveParamV={Ks[j]};
      real_type* curveParam=curveParamV.data();

      TimePerf tp;
      tp.start();
      cout << "\t";
      DP::solveDP<Dubins<double> >(v, discr, fixedAngles, curveParam, false);
      auto time=tp.getTime();
      cout << "\tExample " << j << " completed in " << time << " ms" << endl;
    }
  }
#else
#define KAYA kaya1
  std::vector<bool> fixedAngles;
  for (int i=0; i<KAYA.size(); i++){
    if (i==0 || i==KAYA.size()-1) {
      fixedAngles.push_back(true);
    }
    else {
      fixedAngles.push_back(false);
    }
  }
  std::vector<real_type> curveParamV={3.0};
  real_type* curveParam=curveParamV.data();

  DP::solveDP<Dubins<double> >(KAYA, DISCR, fixedAngles, curveParam, false);
#endif
  //READ_FROM_FILE()
  //CLOSE_FILE()
  return 0;
}
