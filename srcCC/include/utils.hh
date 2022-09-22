#ifndef UTILS_HH
#define UTILS_HH

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <limits>

#include <typedefs.hh>

//#define DEBUG

#ifdef DEBUG
#define COUT(x) cout << x;
#else
#define COUT(x)
#endif

#ifndef ASSERT
  #define ASSERT(COND,MSG)         \
    if ( !(COND) ) {                        \
      std::ostringstream ost ;              \
      ost << "On line: " << __LINE__        \
          << " file: " << __FILE__          \
          << MSG << '\n' ;          \
      throw std::runtime_error(ost.str()) ; \
    }
#endif //ASSERT

extern real_type const epsi        ;
extern real_type const m_pi        ; // pi
extern real_type const m_pi_2      ; // pi/2
extern real_type const m_2pi       ; // 2*pi
extern real_type const m_1_pi      ; // 1/pi
extern real_type const m_1_sqrt_pi ; // 1/sqrt(pi)

#define printCV(v, d)       \
  printf("<");              \
  for (uint i=0; i<d; i++){ \
    printf("%s ", v[i]);    \
  }                         \
  printf("\n");

#define printV(v)                         \
  std::cout << "<";                       \
  for (auto a : v) std::cout << a << " "; \
  std::cout << ">" << endl;

#define printM(M, discr, size) \
for (int i=0; i<discr; i++){   \
  cout << "th" << i;           \
  for (int j=0; j<size; j++){  \
    cout << M[i][j] << "\t";   \
  }                            \
  cout << endl;                \
}

#define printVM(M, discr, size)       \
for (int i=0; i<discr; i++){          \
  std::cout << "th" << i;             \
  for (int j=0; j<size; j++){         \
    std::cout << std::setw(30);       \
    std::cout << M[i*size+j] << "\t"; \
  }                                   \
  std::cout << std::endl;             \
}

#define printCVM(M, discr, size)             \
for (int i=0; i<discr; i++){                 \
  printf("th%d", i);                         \
  for (int j=0; j<size; j++){                \
    printf("\t%-5f", (double)(M[i*size+j])); \
  }                                          \
  printf("\n");                              \
}

template<class T>
inline T ABS(T x, T y) {return (x>y ? (x-y) : (y-x)); }

template<class T>
inline bool eq(const T x, const T y, const T EPSI=std::numeric_limits<T>::epsilon()) {
  return ((ABS(x, y)>(EPSI)) ? false : true);
}

#endif //UTILS_HH
