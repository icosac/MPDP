#ifndef TYPEDEFS_HH
#define TYPEDEFS_HH

typedef double real_type;  ///<Typedef to describe the real type
typedef int int_type       ///<Typedef to describe integers
typedef unsigned int uint; ///<Typedef to abbreviate unsigned int
typedef real_type Angle;   ///<Typedef if the Angle. This should have been a class.
typedef double 	LEN_T;     ///<Typedef to describe the length
typedef double 	K_T;       ///<Typedef to describe the curvature

enum ANGLE { FREE = 0 };

#ifdef CUDA_ON
#define BOTH __host__ __device__
#else
#define BOTH
#endif

#endif //TYPEDEFS_HH
