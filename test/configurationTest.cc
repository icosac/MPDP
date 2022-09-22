#include <configuration.hh>
#include <utils.hh>
#include <cmath>

#if defined(BOOST)
#define BOOST_TEST_MODULE Configuration
#include <boost/test/unit_test.hpp>
BOOST_AUTO_TEST_SUITE(ConstructorsTest)
BOOST_AUTO_TEST_CASE(ConstructorsTestInt){
  Configuration2<int> c(1,2,3.14);
  if (!eq<int>(c.x(), 1)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<int>(c.y(), 2)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),3.14)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 0.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}

BOOST_AUTO_TEST_CASE(ConstructorsTestDouble){
  Configuration2<double> c((double)1.3,(double)2.4,3.14);
  if (!eq<double>(c.x(),1.3)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<double>(c.y(),2.4)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),3.14)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 0.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}

BOOST_AUTO_TEST_CASE(ConstructorsTestFloat){
  Configuration2<float> c((float)1.3,(float)2.4,3.14);
  if (!eq<float>(c.x(),1.3)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<float>(c.y(),2.4)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),3.14)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 0.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SettingTest)
BOOST_AUTO_TEST_CASE(SettingTestInt){
  Configuration2<int> c(1,2,3.14);
  c.x(2); c.y(c.y()+1); c.th(2*M_PI); c.k(1);
  if (!eq<int>(c.x(), 2)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<int>(c.y(), 3)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),2*M_PI)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 1.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}

BOOST_AUTO_TEST_CASE(SettingTestFloat){
  Configuration2<float> c(1.0,2.0,3.14);
  c.x(2.2); c.y(c.y()+1.3); c.th(2*M_PI); c.k(1);
  if (!eq<float>(c.x(), 2.2)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<float>(c.y(), 3.3)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),2*M_PI)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 1.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}

BOOST_AUTO_TEST_CASE(SettingTestDouble){
  Configuration2<double> c(1,2,3.14);
  c.x(2.0); c.y(c.y()+1.1); c.th(2*M_PI); c.k(1);
  if (!eq<double>(c.x(), 2.0)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<double>(c.y(), 3.1)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),2*M_PI)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 1.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CopyTest)
BOOST_AUTO_TEST_CASE(CopyDoubleToInt){
  Configuration2<int> c(1,2,3.14);
  Configuration2<double> c1(2.8, 3.4, M_PI, 1.0);
  c=c1;
  if (!eq<int>(c.x(), 2)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<int>(c.y(), 3)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),M_PI)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 1.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}

BOOST_AUTO_TEST_CASE(CopyFloatToInt){
  Configuration2<int> c(1,2,3.14);
  Configuration2<float> c1(2.8, 3.4, M_PI, 1.0);
  c=c1;
  if (!eq<int>(c.x(), 2)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<int>(c.y(), 3)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),M_PI)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 1.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}

BOOST_AUTO_TEST_CASE(CopyIntToInt){
  Configuration2<int> c(1,2,3.14);
  Configuration2<int> c1(2, 3, M_PI, 2.0);
  c=c1;
  if (!eq<int>(c.x(), 2)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<int>(c.y(), 3)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),M_PI)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 2.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}

BOOST_AUTO_TEST_CASE(CopyDoubleToDouble){
  Configuration2<double> c(1.1,2.2,3.14);
  Configuration2<double> c1(2.8, 3.4, M_PI, 1.0);
  c=c1;
  if (!eq<int>(c.x(), 2)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<int>(c.y(), 3)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),M_PI)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 1.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}

BOOST_AUTO_TEST_CASE(CopyFloatToDouble){
  Configuration2<double> c(1.1,2.2,3.14);
  Configuration2<float> c1(2.8, 3.4, M_PI, 1.0);
  c=c1;
  if (!eq<int>(c.x(), 2)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<int>(c.y(), 3)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),M_PI)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 1.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}

BOOST_AUTO_TEST_CASE(CopyIntToDouble){
  Configuration2<double> c(1.1,2.2,3.14);
  Configuration2<int> c1(2, 3, M_PI, 2.0);
  c=c1;
  if (!eq<int>(c.x(), 2)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<int>(c.y(), 3)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),M_PI)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 2.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}

BOOST_AUTO_TEST_CASE(CopyDoubleToFloat){
  Configuration2<float> c(1.1,2.2,3.14);
  Configuration2<double> c1(2.8, 3.4, M_PI, 1.0);
  c=c1;
  if (!eq<int>(c.x(), 2)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<int>(c.y(), 3)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),M_PI)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 1.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}

BOOST_AUTO_TEST_CASE(CopyFloatToFloat){
  Configuration2<float> c(1.1,2.2,3.14);
  Configuration2<float> c1(2.8, 3.4, M_PI, 1.0);
  c=c1;
  if (!eq<int>(c.x(), 2)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<int>(c.y(), 3)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),M_PI)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 1.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}

BOOST_AUTO_TEST_CASE(CopyIntToFloat){
  Configuration2<float> c(1.1,2.2,3.14);
  Configuration2<int> c1(2, 3, M_PI, 2.0);
  c=c1;
  if (!eq<int>(c.x(), 2)){ BOOST_ERROR("Configuration x does not match."); }
  if (!eq<int>(c.y(), 3)){ BOOST_ERROR("Configuration y does not match."); }
  if (!eq<Angle>(c.th(),M_PI)){ BOOST_ERROR("Configuration angle does not match."); }
  if (!eq<K_T>(c.k(), 2.0)){ BOOST_ERROR("Configuration curvature does not match."); }
}
BOOST_AUTO_TEST_SUITE_END()

#elif defined(GTEST)
#include <gtest/gtest.h>
TEST(Configuration2, constructorInt){
  Configuration2<int> c(1, 2, 3.14);
  EXPECT_EQ(1, c.x());
  EXPECT_EQ(2, c.y());
  EXPECT_FLOAT_EQ(3.14, c.th());
  EXPECT_FLOAT_EQ(0.0, c.k());
}

TEST(Configuration2, constructorDouble){
  Configuration2<double> c((double)1.0, (double)2.0, 3.14);
  EXPECT_FLOAT_EQ(1.0, c.x());
  EXPECT_FLOAT_EQ(2.0, c.y());
  EXPECT_FLOAT_EQ(3.14, c.th());
  EXPECT_FLOAT_EQ(0.0, c.k());
}

TEST(Configuration2, constructorFloat){
  Configuration2<float> c((float)1.0, (float)2.0, 3.14);
  EXPECT_FLOAT_EQ(1.0, c.x());
  EXPECT_FLOAT_EQ(2.0, c.y());
  EXPECT_FLOAT_EQ(3.14, c.th());
  EXPECT_FLOAT_EQ(0.0, c.k());
}

TEST(Configuration2, SettingTestInt){
  Configuration2<int> c(1,2,3.14);
  c.x(2); c.y(c.y()+1); c.th(2*M_PI); c.k(1);
  EXPECT_EQ(2, c.x());
  EXPECT_EQ(3, c.y());
  EXPECT_FLOAT_EQ(2*M_PI, c.th());
  EXPECT_FLOAT_EQ(1.0, c.k());
}

TEST(Configuration2, SettingTestDouble){
  Configuration2<double> c((double)1.3,(double)2.6,3.14);
  c.x(2.5); c.y(c.y()+1.3); c.th(2*M_PI); c.k(1.5);
  EXPECT_FLOAT_EQ(2.5, c.x()) << "Configuration x does not match.";
  EXPECT_FLOAT_EQ(3.9, c.y()) << "Configuration y does not match.";
  EXPECT_FLOAT_EQ(2*M_PI, c.th()) << "Configuration angle does not match.";
  EXPECT_FLOAT_EQ(1.5, c.k()) << "Configuration curvature does not match.";
}

TEST(Configuration2, SettingTestfloat){
  Configuration2<float> c((float)1.3,(float)2.6,3.14);
  c.x(2.5); c.y(c.y()+1.3); c.th(2*M_PI); c.k(1.5);
  EXPECT_FLOAT_EQ(2.5, c.x());
  EXPECT_FLOAT_EQ(3.9, c.y());
  EXPECT_FLOAT_EQ(2*M_PI, c.th());
  EXPECT_FLOAT_EQ(1.5, c.k());
}

TEST(Configuration2, CopyTestDoubleToInt){
  Configuration2<int> c(1,2,3.14);
  Configuration2<double> c1(2.2,3.1,2*M_PI,1.0);
  c=c1;
  std::cout << c << std::endl;
  EXPECT_EQ(2, c.x());
  EXPECT_EQ(3, c.y());
  EXPECT_FLOAT_EQ(2*M_PI, c.th());
  EXPECT_FLOAT_EQ(1.0, c.k());
}

TEST(Configuration2, CopyTestFloatToInt){
  Configuration2<int> c(1,2,3.14);
  Configuration2<float> c1((float)2.2,(float)3.3,2*M_PI,1.2);
  c=c1;
  EXPECT_EQ(2, c.x());
  EXPECT_EQ(3, c.y());
  EXPECT_FLOAT_EQ(2*M_PI, c.th());
  EXPECT_FLOAT_EQ(1.2, c.k());
}

TEST(Configuration2, CopyTestIntToInt){
  Configuration2<int> c(1,2,3.14);
  Configuration2<double> c1(2,3,2*M_PI,1.0);
  c=c1;
  std::cout << c << std::endl;
  EXPECT_EQ(2, c.x());
  EXPECT_EQ(3, c.y());
  EXPECT_FLOAT_EQ(2*M_PI, c.th());
  EXPECT_FLOAT_EQ(1.0, c.k());
}

TEST(Configuration2, CopyTestIntToDouble){
  Configuration2<int> c(1,2,3.14);
  Configuration2<double> c1(2.2,3.1,2*M_PI,1.0);
  c=c1;
  std::cout << c << std::endl;
  EXPECT_EQ(2, c.x());
  EXPECT_EQ(3, c.y());
  EXPECT_FLOAT_EQ(2*M_PI, c.th());
  EXPECT_FLOAT_EQ(1.0, c.k());
}

TEST(Configuration2, CopyTestFloatToDouble){
  Configuration2<double> c((double)2.2,(double)3.3,2*M_PI,1.2);
  Configuration2<float> c1((float)1.2,(float)3.6,2*M_PI,1.2);
  c=c1;
  EXPECT_FLOAT_EQ(1.2, c.x());
  EXPECT_FLOAT_EQ(3.6, c.y());
  EXPECT_FLOAT_EQ(2*M_PI, c.th());
  EXPECT_FLOAT_EQ(1.2, c.k());
}

TEST(Configuration2, CopyTestDoubleToDouble){
  Configuration2<double> c((double)2.2,(double)3.3,2*M_PI,1.2);
  Configuration2<double> c1((double)1.2,(double)3.6,2*M_PI,1.2);
  c=c1;
  EXPECT_FLOAT_EQ(1.2, c.x());
  EXPECT_FLOAT_EQ(3.6, c.y());
  EXPECT_FLOAT_EQ(2*M_PI, c.th());
  EXPECT_FLOAT_EQ(1.2, c.k());
}

TEST(Configuration2, CopyTestIntToFloat){
  Configuration2<float> c(1.4,2.3,3.14);
  Configuration2<int> c1(2,3,2*M_PI,1.0);
  c=c1;
  std::cout << c << std::endl;
  EXPECT_FLOAT_EQ(2.0, c.x());
  EXPECT_FLOAT_EQ(3.0, c.y());
  EXPECT_FLOAT_EQ(2*M_PI, c.th());
  EXPECT_FLOAT_EQ(1.0, c.k());
}

TEST(Configuration2, CopyTestFloatToFloat){
  Configuration2<float> c((float)2.2,(float)3.3,2*M_PI,1.2);
  Configuration2<float> c1((float)1.2,(float)3.6,2*M_PI,1.2);
  c=c1;
  EXPECT_FLOAT_EQ(1.2, c.x());
  EXPECT_FLOAT_EQ(3.6, c.y());
  EXPECT_FLOAT_EQ(2*M_PI, c.th());
  EXPECT_FLOAT_EQ(1.2, c.k());
}

TEST(Configuration2, CopyTestDoubleToFloat){
  Configuration2<float> c((float)2.2,(float)3.3,2*M_PI,1.2);
  Configuration2<double> c1((double)1.2,(double)3.6,2*M_PI,1.2);
  c=c1;
  EXPECT_FLOAT_EQ(1.2, c.x());
  EXPECT_FLOAT_EQ(3.6, c.y());
  EXPECT_FLOAT_EQ(2*M_PI, c.th());
  EXPECT_FLOAT_EQ(1.2, c.k());
}

#endif




