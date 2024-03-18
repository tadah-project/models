#include "catch2/catch.hpp"
#include "../descriptors/dm/dm_all.h"
#include "../../CORE/maths.h"
#include "../cutoffs/cutoffs.h"
#include <vector>

template <typename DM>
void central_difference_mb(DM *dm, std::vector<Vec3d> &positions, 
                           size_t j, fd_type &fdij, double h=1e-10) {

  //size_t i=0;
  std::vector<size_t> indices = {0,1,2,3};
  indices.erase(indices.begin() + j); // remove jth index
  indices.erase(indices.begin()); // remove 0 index

  Vec3d &x1 = positions[0]; 

  Vec3d del13 = x1-positions[indices.back()];
  indices.pop_back();
  Vec3d del14 = x1-positions[indices.back()];
  indices.pop_back();

  for (int i=0; i<3; ++i) {
    Vec3d xjn = positions[j]; // copy
    Vec3d xjp = positions[j]; // copy
    xjn[i]-=h;
    xjp[i]+=h;

    Vec3d del12n = x1-xjn;
    Vec3d del12p = x1-xjp;
    aed_type2 rhoip;
    aed_type2 rhoin;

    dm->init_rhoi(rhoip);
    dm->init_rhoi(rhoin);

    dm->calc_rho_mb(del12p, rhoip);
    dm->calc_rho_mb(del13, rhoip);
    dm->calc_rho_mb(del14, rhoip);

    dm->calc_rho_mb(del12n, rhoin);
    dm->calc_rho_mb(del13, rhoin);
    dm->calc_rho_mb(del14, rhoin);

    aed_type2 aedip(fdij.rows());
    aed_type2 aedin(fdij.rows());
    aedip.set_zero();
    aedin.set_zero();
    dm->calc_aed(rhoip, aedip);
    dm->calc_aed(rhoin, aedin);

    fdij[i] = -(0.5/h)*(aedip - aedin);
  }
}

template<typename DM>
void test_case(Config &config_mb, 
               std::string label,
               size_t size) {
  // Define a simple system of 4 atoms
  // Will calculate forces acting on the first atom
  Vec3d x1(0.0,0.0,0.0);
  Vec3d x2(1.1,0.3,0.4);
  Vec3d x3(2.1,1.3,1.4);
  Vec3d x4(-1.11,0.33,-0.44);
  std::vector<Vec3d> positions = {x1,x2,x3,x4};

  Vec3d del12 = x1-x2;
  Vec3d del13 = x1-x3;
  Vec3d del14 = x1-x4;

  double rsq12 = del12*del12;
  double rsq13 = del13*del13;
  double rsq14 = del14*del14;
  double r12 = sqrt(rsq12);
  double r13 = sqrt(rsq13);
  double r14 = sqrt(rsq14);

  DM dm(config_mb);
  aed_type2 aed1(dm.size());
  aed1.set_zero();
  rho_type rho1;
  rho1.set_zero();
  fd_type fd12(dm.size());
  fd_type fd13(dm.size());
  fd_type fd14(dm.size());
  fd12.set_zero();
  fd13.set_zero();
  fd14.set_zero();
  dm.init_rhoi(rho1);
  dm.calc_rho_mb(del12, rho1);
  dm.calc_rho_mb(del13, rho1);
  dm.calc_rho_mb(del14, rho1);
  dm.calc_aed(rho1,aed1);

  int mode = 
  dm.calc_dXijdri(r12,rsq12,del12,1.0,0.0,rho1,fd12);
  dm.calc_dXijdri(r13,rsq13,del13,1.0,0.0,rho1,fd13);
  dm.calc_dXijdri(r14,rsq14,del14,1.0,0.0,rho1,fd14);
  if (mode==0) {
    fd12[0] *= del12[0]/r12;
    fd13[0] *= del13[0]/r13;
    fd14[0] *= del14[0]/r14;
  }

  // Calculate approximation
  fd_type fd12a(dm.size());
  fd_type fd13a(dm.size());
  fd_type fd14a(dm.size());
  fd12a.set_zero();
  fd13a.set_zero();
  fd14a.set_zero();
  central_difference_mb(&dm,positions, 1, fd12a, 1e-8);
  central_difference_mb(&dm,positions, 2, fd13a, 1e-8);
  central_difference_mb(&dm,positions, 3, fd14a, 1e-8);

  std::vector<double> fx12(fd12[0].data(),fd12[0].data()+fd12[0].size());
  std::vector<double> fy12(fd12[1].data(),fd12[1].data()+fd12[1].size());
  std::vector<double> fz12(fd12[2].data(),fd12[2].data()+fd12[2].size());
  std::vector<double> fx12a(fd12a[0].data(),fd12a[0].data()+fd12a[0].size());
  std::vector<double> fy12a(fd12a[1].data(),fd12a[1].data()+fd12a[1].size());
  std::vector<double> fz12a(fd12a[2].data(),fd12a[2].data()+fd12a[2].size());
  REQUIRE_THAT(fx12a, Catch::Matchers::Approx(fx12).margin(1e-5));
  REQUIRE_THAT(fy12a, Catch::Matchers::Approx(fy12).margin(1e-5));
  REQUIRE_THAT(fz12a, Catch::Matchers::Approx(fz12).margin(1e-5));

  std::vector<double> fx13(fd13[0].data(),fd13[0].data()+fd13[0].size());
  std::vector<double> fy13(fd13[1].data(),fd13[1].data()+fd13[1].size());
  std::vector<double> fz13(fd13[2].data(),fd13[2].data()+fd13[2].size());
  std::vector<double> fx13a(fd13a[0].data(),fd13a[0].data()+fd13a[0].size());
  std::vector<double> fy13a(fd13a[1].data(),fd13a[1].data()+fd13a[1].size());
  std::vector<double> fz13a(fd13a[2].data(),fd13a[2].data()+fd13a[2].size());
  REQUIRE_THAT(fx13a, Catch::Matchers::Approx(fx13).margin(1e-5));
  REQUIRE_THAT(fy13a, Catch::Matchers::Approx(fy13).margin(1e-5));
  REQUIRE_THAT(fz13a, Catch::Matchers::Approx(fz13).margin(1e-5));

  std::vector<double> fx14(fd14[0].data(),fd14[0].data()+fd14[0].size());
  std::vector<double> fy14(fd14[1].data(),fd14[1].data()+fd14[1].size());
  std::vector<double> fz14(fd14[2].data(),fd14[2].data()+fd14[2].size());
  std::vector<double> fx14a(fd14a[0].data(),fd14a[0].data()+fd14a[0].size());
  std::vector<double> fy14a(fd14a[1].data(),fd14a[1].data()+fd14a[1].size());
  std::vector<double> fz14a(fd14a[2].data(),fd14a[2].data()+fd14a[2].size());
  REQUIRE_THAT(fx14a, Catch::Matchers::Approx(fx14).margin(1e-5));
  REQUIRE_THAT(fy14a, Catch::Matchers::Approx(fy14).margin(1e-5));
  REQUIRE_THAT(fz14a, Catch::Matchers::Approx(fz14).margin(1e-5));

  REQUIRE(dm.label() == label);
  REQUIRE(dm.size() == size);
}
TEST_CASE( "Testing DM_Dummy", "[DM]" ) {
}
TEST_CASE( "Testing DM_EAD 1", "[DM]" ) {
  Config config_mb;
  config_mb.add("RCUTMB",10);
  config_mb.remove("INITMB");
  config_mb.add("INITMB","true");
  config_mb.add("AGRIDMB","2");
  config_mb.add("CGRIDMB","1 2");
  config_mb.add("SGRIDMB","0.1 0.2");
  test_case<DM_EAD>(config_mb, "DM_EAD", 3*2);
}
TEST_CASE( "Testing DM_EAD 2", "[DM]" ) {
  Config config_mb;
  config_mb.add("RCUTMB",6);
  config_mb.remove("INITMB");
  config_mb.add("INITMB","true");
  config_mb.add("AGRIDMB","3");
  config_mb.add("CGRIDMB","1 2 3");
  config_mb.add("SGRIDMB","0.1 0.2 0.9");
  test_case<DM_EAD>(config_mb, "DM_EAD", 4*3);
}
TEST_CASE( "Testing DM_Blip 1", "[DM]" ) {
  Config config_mb;
  config_mb.add("RCUTMB",10);
  config_mb.remove("INITMB");
  config_mb.add("INITMB","true");
  config_mb.add("AGRIDMB","2");
  config_mb.add("CGRIDMB","1 2");
  config_mb.add("SGRIDMB","0.1 0.2");
  test_case<DM_Blip>(config_mb, "DM_Blip", 3*2);
}
TEST_CASE( "Testing DM_Blip 2", "[DM]" ) {
  Config config_mb;
  config_mb.add("RCUTMB",6);
  config_mb.remove("INITMB");
  config_mb.add("INITMB","true");
  config_mb.add("AGRIDMB","3");
  config_mb.add("CGRIDMB","1 2 3");
  config_mb.add("SGRIDMB","0.1 0.2 0.9");
  test_case<DM_Blip>(config_mb, "DM_Blip", 4*3);
}
