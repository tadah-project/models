#include "catch2/catch.hpp"
#include "../descriptors/d2/d2_all.h"
#include "../../CORE/maths.h"
#include "../cutoffs/cutoffs.h"

template<typename D, typename C>
void test_d2(Config &c, double r) {

  D d(c);
  double rcut = c.get<double>("RCUT2B");
  C fc(rcut);
  aed_type2 aed(d.size());
  d.calc_aed(r,r*r,fc.calc(r),aed);
  fd_type fd(d.size());
  d.calc_dXijdri(r,r*r,fc.calc(r),fc.calc_prime(r),fd);
  fd_type fda(d.size());
  d.calc_fd_approx(&fc, r, fda, 1e-6);

  // Convert to std::vector, as catch2 matchers require these
  std::vector<double> fx(fd[0].data(),fd[0].data()+fd[0].size());
  std::vector<double> fy(fd[1].data(),fd[1].data()+fd[1].size());
  std::vector<double> fz(fd[2].data(),fd[2].data()+fd[2].size());

  std::vector<double> fxa(fda[0].data(),
                          fda[0].data()+fda[0].size());
  std::vector<double> fya(fda[1].data(),
                          fda[1].data()+fda[1].size());
  std::vector<double> fza(fda[2].data(),
                          fda[2].data()+fda[2].size());

  REQUIRE_THAT(fxa, Catch::Matchers::Approx(fx).margin(1e-5));
  // y and z should be zero
  REQUIRE_THAT(fya, Catch::Matchers::Approx(fy));
  REQUIRE_THAT(fza, Catch::Matchers::Approx(fz));

  aed_type2 aed2(d.size());
  fd_type fd2(d.size());
  d.calc_all(r,r*r,fc.calc(r),fc.calc_prime(r),aed2,fd2);
  std::vector<double> fx2(fd2[0].data(),fd2[0].data()+fd2[0].size());
  std::vector<double> fy2(fd2[1].data(),fd2[1].data()+fd2[1].size());
  std::vector<double> fz2(fd2[2].data(),fd2[2].data()+fd2[2].size());

  REQUIRE_THAT(fx2, Catch::Matchers::Approx(fx).epsilon(1e-12));
  REQUIRE_THAT(fy2, Catch::Matchers::Approx(fy).epsilon(1e-12));
  REQUIRE_THAT(fz2, Catch::Matchers::Approx(fz).epsilon(1e-12));
}

template<typename D>
void test_D2(Config &c) {
  //c.remove("VERBOSE");
  //c.add("VERBOSE", 1);
  // We cannot test r=rcut as finite difference will fail
  // but we can test r slightly smaller than rcut
  std::vector<double> rs = {2,2.0001,2.0123,5,11};
  std::vector<double> rcuts {10, 2.1, 2.5, 3.5,6.7};

  for (const auto &r: rs) {
    for (const auto &rcut: rcuts) {
      c.remove("RCUT2B");
      c.add("RCUT2B", rcut);
      test_d2<D,Cut_Dummy>(c,r);
      test_d2<D,Cut_Cos>(c,r);
      test_d2<D,Cut_Tanh>(c,r);
      test_d2<D,Cut_Poly2>(c,r);
    }
  }
  Config c2; // INIT2B=false
  D d(c2);
  REQUIRE(d.size()==0);
}

TEST_CASE( "Testing D2_Dummy", "[D2]" ) {
  Config c;
  double rcut=5.1;
  c.remove("INIT2B");
  c.add("INIT2B", "true");
  c.add("RCUT2B", rcut);
  test_d2<D2_Dummy,Cut_Dummy>(c,1);

  REQUIRE_NOTHROW(D2_Dummy());
  REQUIRE(D2_Dummy().label()=="D2_Dummy");
}
TEST_CASE( "Testing D2_LJ", "[D2]" ) {
  Config c;
  c.remove("INIT2B");
  c.add("INIT2B", "true");
  test_D2<D2_LJ>(c);

  D2_LJ d(c);
  REQUIRE(d.size()==2);
  REQUIRE(d.label()=="D2_LJ");
}
TEST_CASE( "Testing D2_BLIP", "[D2]" ) {
  using D = D2_Blip;
  Config c;
  c.remove("INIT2B");
  c.add("INIT2B", "true");
  c.add("SGRID2B","0.01 0.4 5.6");
  c.add("CGRID2B","0.0 1.2 3.4");
  test_D2<D>(c);
  D d(c);
  REQUIRE(d.size()==3);
  REQUIRE(d.label()=="D2_Blip");

  // test throws
  c.remove("SGRID2B");
  c.remove("CGRID2B");
  c.add("SGRID2B","0.4 5.6");
  c.add("CGRID2B","0.0 1.2 3.4");
  REQUIRE_THROWS(D(c));
  c.remove("SGRID2B");
  c.remove("CGRID2B");
  c.add("SGRID2B","0.4 5.6");
  c.add("CGRID2B","0.0");
  REQUIRE_THROWS(D(c));

  c.remove("SGRID2B");
  c.remove("CGRID2B");
  c.add("SGRID2B","0.0 0.4 5.6");
  c.add("CGRID2B","0.0 1.2 3.4");
  REQUIRE_THROWS(D(c));
}
TEST_CASE( "Testing D2_BP", "[D2]" ) {
  using D = D2_BP;
  Config c;
  c.remove("INIT2B");
  c.add("INIT2B", "true");
  c.add("SGRID2B","0.01 0.4");
  c.add("CGRID2B","0.0 1.2");
  test_D2<D>(c);
  D d(c);
  REQUIRE(d.size()==2);
  REQUIRE(d.label()=="D2_BP");
  // test throws
  c.remove("SGRID2B");
  c.remove("CGRID2B");
  c.add("SGRID2B","0.4 5.6");
  c.add("CGRID2B","0.0 1.2 3.4");
  REQUIRE_THROWS(D(c));
  c.remove("SGRID2B");
  c.remove("CGRID2B");
  c.add("SGRID2B","0.4 5.6");
  c.add("CGRID2B","0.0");
  REQUIRE_THROWS(D(c));

  c.remove("SGRID2B");
  c.remove("CGRID2B");
  c.add("SGRID2B","0.0 0.4");
  c.add("CGRID2B","0.0 1.2");
  REQUIRE_THROWS(D(c));
}
TEST_CASE( "Testing D2_EAM", "[D2]" ) {
}
TEST_CASE( "Testing D2_MIE", "[D2]" ) {
  Config c;
  c.remove("INIT2B");
  c.add("INIT2B", "true");
  c.add("SGRID2B","7 13");
  test_D2<D2_MIE>(c);
  D2_MIE d(c);
  REQUIRE(d.size()==2);
  REQUIRE(d.label()=="D2_MIE");

  // test throws
  c.remove("SGRID2B");
  c.add("SGRID2B","7");
  REQUIRE_THROWS(D2_MIE(c));

  c.remove("SGRID2B");
  c.add("SGRID2B","7 8 9");
  REQUIRE_THROWS(D2_MIE(c));

  c.remove("SGRID2B");
  c.add("SGRID2B","-8 9");
  REQUIRE_THROWS(D2_MIE(c));

  c.remove("SGRID2B");
  c.add("SGRID2B","8 -9");
  REQUIRE_THROWS(D2_MIE(c));

}
