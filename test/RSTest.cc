/** 
 * @author Enrico Saccon <enrico.saccon@unitn.it>.
 * @file RS.cc
 * @brief 
 */

#include<iostream>
#include<sstream>
#include<fstream>

#include <rs.hh>
#include <utils.hh>

#define DISCR 720
#define EPSILON 1e-10

inline const char* toCString(std::stringstream msg){
  return msg.str().c_str();
}

const std::vector<LEN_T> KayaLengths={3.415578858075, 6.2780346, 11.916212654286, 7.467562181965};

#define READ_FROM_FILE_RS()                                                                                \
  std::ifstream input("test/RSDataset.txt");                                                               \
    real_type th0, th1, kmax, l;                                                                           \
    int nMan, nSeg;                                                                                        \
    std::string sMan;                                                                                      \
    int i=0;                                                                                               \
    while (input >> th0 >> th1 >> kmax >> l >> nMan >> sMan >> nSeg){                                      \
      i++;                                                                                                 \
      Configuration2 ci(-1, 0, th0);                                                                       \
      Configuration2 cf(1, 0, th1);

#define CLOSE_FILE_RS() } input.close();


#if defined(BOOST)
#define BOOST_TEST_MODULE RS
#include <boost/test/unit_test.hpp>
#include <boost/format.hpp>

BOOST_AUTO_TEST_SUITE(SingleRSTest)
BOOST_AUTO_TEST_CASE(RSP2P){
  READ_FROM_FILE_RS()
    Dubins<real_type> d(ci, cf, kmax);
    if (!eq(d.l(), l, 1e-03)){ BOOST_ERROR(boost::format("Length l %1% does not match %2%\n") % d.l() % l); }
    if (!eq(d.s1(), s1, 1e-03)){ BOOST_ERROR(boost::format("Length s1 %1% does not match %2%\n") % d.s1() % s1); }
    if (!eq(d.s2(), s2, 1e-03)){ BOOST_ERROR(boost::format("Length s2 %1% does not match %2%\n") % d.s2() % s2); }
    if (!eq(d.s3(), s3, 1e-03)){ BOOST_ERROR(boost::format("Length s3 %1% does not match %2%\n") % d.s3() % s3); }
    if (!eq(d.k1(), k1, EPSILON)){ BOOST_ERROR(boost::format("Curvature k1 %1% does not match %2%\n") % d.k1() % k1); }
    if (!eq(d.k2(), k2, EPSILON)){ BOOST_ERROR(boost::format("Curvature k2 %1% does not match %2%\n") % d.k2() % k2); }
    if (!eq(d.k3(), k3, EPSILON)){ BOOST_ERROR(boost::format("Curvature k3 %1% does not match %2%\n") % d.k3() % k3); }
  CLOSE_FILE_RS()
}
BOOST_AUTO_TEST_SUITE_END()

#elif defined(GTEST)
#include <gtest/gtest.h>

TEST(RSTest, RSP2P){//TODO Find bug for which ctest shows this test passed, but ./build/RSTest does not pass.
  READ_FROM_FILE_RS()
    RS rs(ci, cf, {kmax});
    rs.solve();
    EXPECT_NEAR(rs.l(), l, EPSILON) << "\nReeds-Shepp length not computed correctly: " << th0 << ", " << th1 << ", "<< kmax << ", " << nMan << ", " << l << "!=" << rs.l() << std::endl;;
//    EXPECT_EQ(rs.getNseg(), nSeg);

    if (nMan != -1 && nMan != rs.getNman()) {
      std::vector<double> debug = {};
      rs.reeds_shepp(-1, &debug);
      EXPECT_NEAR(debug[nMan], rs.l(), EPSILON) << "\nReeds-Shepp Nman not computed correctly: " << th0 << ", " << th1 << ", " << kmax << ", " << rs.l()
                                                << ", " << l << " " << nMan << "!= " << rs.getNman() << std::endl;
    }
  CLOSE_FILE_RS()
}

#endif


