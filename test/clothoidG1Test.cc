#include <clothoidG1.hh>
#include <utils.hh>
#include <fstream>

using namespace std;

//This is similar to what Contiki does. Since the beginning of the loop is similar for both BOOST and GoogleTest, then I can define it 
#define READ_FROM_FILE()                                                  \
  ifstream input("build_clothoid.txt");                                   \
    float x0, y0, th0, x1, y1, th1, k, dk, l;                             \
    int i=0;                                                              \
    while (input >> x0 >> y0 >> th0 >> x1 >> y1 >> th1 >> k >> dk >> l){  \
      i++;                                                                \
      Configuration2<float>ci(x0, y0, th0);                               \
      Configuration2<float>cf(x1, y1, th1);                               \
      ClothoidG1<float>c(ci, cf);

#define CLOSE_FILE() } input.close();

#if defined(BOOST)
#define BOOST_TEST_MODULE ClothoidG1
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(ClothoidG1Test){
//  READ_FROM_FILE()
//    if (!eq<double>(c.l(), l, 0.001)){ BOOST_ERROR("Length does not match"); }
//    if (!eq<double>(c.k(), k)){ BOOST_ERROR("Curvature does not match"); }
//    if (!eq<double>(c.dk(), dk, 0.000001)){ BOOST_ERROR("Curvature derivative does not match"); }
//  CLOSE_FILE()
}

#elif defined(GTEST)
#include <gtest/gtest.h>
TEST(ClothoidG1Test, G1Float){
//    READ_FROM_FILE()
//      if (eq<double>(c.l(), l, 0.1)){ FAIL() << "Length does not match " << c.l() << "!=" << l; }
//      if (!eq<double>(c.k(), k)){ BOOST_ERROR("Curvature does not match"); }
//      if (!eq<double>(c.dk(), dk, 1e-6)){ FAIL() << "Curvature derivative does not match " << c.dk() << "!=" << dk; }
//      EXPECT_FLOAT_EQ(c.l(), l);
//      EXPECT_FLOAT_EQ(c.k(), k);
//      EXPECT_FLOAT_EQ(c.dk(), dk);
//  CLOSE_FILE()
}

#endif
