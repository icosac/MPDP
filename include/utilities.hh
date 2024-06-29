//
// Created by Enrico Saccon on 22/09/22.
//

#ifndef UTILITIES_HH
#define UTILITIES_HH

#include <fstream>

// inline void PrintScientific1D(real_type d){
//   if (d == 0)
//   {
//     printf ("%*d", 6, 0);
//     return;
//   }

//   int exponent  = (int)floor(log10( fabs(d)));  // This will round down the exponent
//   real_type base   = d * pow(10, -1.0*exponent);

//   printf("%1.1lfe%+01d", base, exponent);
// }

// inline void PrintScientific2D(real_type d){
//   if (d == 0)
//   {
//     printf ("%*d", 7, 0);
//     return;
//   }

//   int exponent  = (int)floor(log10( fabs(d)));  // This will round down the exponent
//   real_type base   = d * pow(10, -1.0*exponent);

//   printf("%1.1lfe%+02d", base, exponent);
// }

inline void initAsyFile(std::ofstream& file){
  file << "import graph;\n"
       << "include \"clothoidLib.asylib\";\n"
       << "size(14cm,7cm);\n"
       << "\n\n\n"
       << "path p;\n";
}


#endif //UTILITIES_HH
