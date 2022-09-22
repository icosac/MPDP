#ifndef TYPEDEFS_HH
#define TYPEDEFS_HH

//General
typedef double real_type; ///<Typedef to describe the real type //I would love to put double here, but then CUDA breaks with atan2 not been defined for doubles... really don't know why
typedef int int_type;     ///<Typedef to describe integers
typedef unsigned int uint;///<Typedef to abbreviate unsigned int

//configuration.hh
typedef double Angle;

//curve.hh
typedef double LEN_T;///<Typedef to describe the length
typedef double K_T;  ///<Typedef to describe the curvature
namespace {
  enum CURVE_TYPE { INVALID, CLOTHOID, DUBINS, DUBINS_ARC }; ///< Possible types of CURVE
}

enum ANGLE { INVALID = 0 };

#ifdef CUDA_ON
#define BOTH __host__ __device__
#else
#define BOTH
#endif

#endif //TYPEDEFS_HH
