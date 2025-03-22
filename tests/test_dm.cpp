#include "catch2/catch.hpp"
#include <tadah/core/core_types.h>
#include <tadah/models/descriptors/dm/dm_all.h>
#include <tadah/core/maths.h>
#include <tadah/models/cutoffs.h>
#include <cmath>
#include <vector>

using namespace tadah::core;
using namespace tadah::models;
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
    aed_type rhoip;
    aed_type rhoin;

    dm->init_rhoi(rhoip);
    dm->init_rhoi(rhoin);

    dm->calc_rho_mb(1,1,del12p, rhoip);
    dm->calc_rho_mb(1,1,del13,  rhoip);
    dm->calc_rho_mb(1,1,del14,  rhoip);

    dm->calc_rho_mb(1,1,del12n, rhoin);
    dm->calc_rho_mb(1,1,del13, rhoin ); 
    dm->calc_rho_mb(1,1,del14, rhoin );

    aed_type aedip(fdij.rows());
    aed_type aedin(fdij.rows());
    aedip.set_zero();
    aedin.set_zero();
    dm->calc_aed(rhoip, aedip);
    dm->calc_aed(rhoin, aedin);

    fdij[i] = -(0.5/h)*(aedip - aedin);
  }
}

template<typename DM, typename CUT>
void test_case(Context &context, 
               std::string label,
               size_t size) {
  context.add("ATOM", "H");
  context.add("WATOM", 1);
  context.add("TYPEMB", "H");
  context.add("TYPEMB", "H");
  // Define a simple system of 4 atoms
  // Will calculate forces acting on the first atom
  Vec3d x1(0.0,0.0,0.0);
  Vec3d x2(2.1,3.3,4.4);
  Vec3d x3(5.1,4.3,3.4);
  Vec3d x4(-4.11,5.33,-3.44);
  std::vector<Vec3d> positions = {x1,x2,x3,x4};

  Vec3d del12 = x1-x2; Vec3d del21 = x2-x1;
  Vec3d del13 = x1-x3; Vec3d del23 = x2-x3;
  Vec3d del14 = x1-x4; Vec3d del24 = x2-x4;

  Vec3d del31 = x3-x1; Vec3d del41 = x4-x1;
  Vec3d del32 = x3-x2; Vec3d del42 = x4-x2;
  Vec3d del34 = x3-x4; Vec3d del43 = x4-x3;

  double rsq12 = del12*del12;  double rsq21 = del21*del21;
  double rsq13 = del13*del13;  double rsq23 = del23*del23;
  double rsq14 = del14*del14;  double rsq24 = del24*del24;
  double r12 = sqrt(rsq12);    double r21 = sqrt(rsq21);
  double r13 = sqrt(rsq13);    double r23 = sqrt(rsq23);
  double r14 = sqrt(rsq14);    double r24 = sqrt(rsq24);

  double rsq31 = del31*del31;  double rsq41 = del41*del41;
  double rsq32 = del32*del32;  double rsq42 = del42*del42;
  double rsq34 = del34*del34;  double rsq43 = del43*del43;
  double r31 = sqrt(rsq31);    double r41 = sqrt(rsq41);
  double r32 = sqrt(rsq32);    double r42 = sqrt(rsq42);
  double r34 = sqrt(rsq34);    double r43 = sqrt(rsq43);

  CUT fc(context.get<double>("RCUTMB"));
  DM dm(context);
  dm.set_fcut(&fc,false);
  aed_type aed1(dm.size());
  aed_type aed2(dm.size());
  aed_type aed3(dm.size());
  aed_type aed4(dm.size());
  aed1.set_zero();
  aed2.set_zero();
  aed3.set_zero();
  aed4.set_zero();
  rho_type rho1;
  rho_type rho2;
  rho_type rho3;
  rho_type rho4;
  rho1.set_zero();
  rho2.set_zero();
  rho3.set_zero();
  rho4.set_zero();
  fd_type fd12(dm.size());
  fd_type fd13(dm.size());
  fd_type fd14(dm.size());

  fd_type fd21(dm.size());
  fd_type fd23(dm.size());
  fd_type fd24(dm.size());

  fd_type fd31(dm.size());
  fd_type fd32(dm.size());
  fd_type fd34(dm.size());

  fd_type fd41(dm.size());
  fd_type fd42(dm.size());
  fd_type fd43(dm.size());

  fd12.set_zero();
  fd13.set_zero();
  fd14.set_zero();

  fd21.set_zero();
  fd23.set_zero();
  fd24.set_zero();

  fd31.set_zero();
  fd32.set_zero();
  fd34.set_zero();

  fd41.set_zero();
  fd42.set_zero();
  fd43.set_zero();

  dm.init_rhoi(rho1);
  dm.init_rhoi(rho2);
  dm.init_rhoi(rho3);
  dm.init_rhoi(rho4);

  dm.calc_rho_mb(1,1,del12, rho1);
  dm.calc_rho_mb(1,1,del13, rho1);
  dm.calc_rho_mb(1,1,del14, rho1);

  dm.calc_rho_mb(1,1,del21, rho2);
  dm.calc_rho_mb(1,1,del23, rho2);
  dm.calc_rho_mb(1,1,del24, rho2);

  dm.calc_rho_mb(1,1,del31, rho3);
  dm.calc_rho_mb(1,1,del32, rho3);
  dm.calc_rho_mb(1,1,del34, rho3);

  dm.calc_rho_mb(1,1,del41, rho4);
  dm.calc_rho_mb(1,1,del42, rho4);
  dm.calc_rho_mb(1,1,del43, rho4);

  dm.calc_aed(rho1,aed1);
  dm.calc_aed(rho2,aed2);
  dm.calc_aed(rho3,aed3);
  dm.calc_aed(rho4,aed4);

  dm.calc_dXijdri(1,1,r12,rsq12,del12,rho1,fd12);
  dm.calc_dXijdri(1,1,r13,rsq13,del13,rho1,fd13);
  dm.calc_dXijdri(1,1,r14,rsq14,del14,rho1,fd14);

  dm.calc_dXijdri(1,1,r21,rsq21,del21,rho2,fd21);
  dm.calc_dXijdri(1,1,r23,rsq23,del23,rho2,fd23);
  dm.calc_dXijdri(1,1,r24,rsq24,del24,rho2,fd24);

  dm.calc_dXijdri(1,1,r31,rsq31,del31,rho3,fd31);
  dm.calc_dXijdri(1,1,r32,rsq32,del32,rho3,fd32);
  dm.calc_dXijdri(1,1,r34,rsq34,del34,rho3,fd34);

  dm.calc_dXijdri(1,1,r41,rsq41,del41,rho4,fd41);
  dm.calc_dXijdri(1,1,r42,rsq42,del42,rho4,fd42);
  dm.calc_dXijdri(1,1,r43,rsq43,del43,rho4,fd43);

  fd_type FD12(dm.size());
  fd_type FD13(dm.size());
  fd_type FD14(dm.size());

  fd_type FD21(dm.size());
  fd_type FD23(dm.size());
  fd_type FD24(dm.size());

  fd_type FD31(dm.size());
  fd_type FD32(dm.size());
  fd_type FD34(dm.size());

  fd_type FD41(dm.size());
  fd_type FD42(dm.size());
  fd_type FD43(dm.size());

  FD12.set_zero();
  FD13.set_zero();
  FD14.set_zero();

  FD21.set_zero();
  FD23.set_zero();
  FD24.set_zero();

  FD31.set_zero();
  FD32.set_zero();
  FD34.set_zero();

  FD41.set_zero();
  FD42.set_zero();
  FD43.set_zero();
  dm.calc_dXijdri_dXjidri(1,1,r12,rsq12,del12,rho1,rho2,FD12,1);
  dm.calc_dXijdri_dXjidri(1,1,r13,rsq13,del13,rho1,rho3,FD13,1);
  dm.calc_dXijdri_dXjidri(1,1,r14,rsq14,del14,rho1,rho4,FD14,1);

  dm.calc_dXijdri_dXjidri(1,1,r21,rsq21,del21,rho2,rho1,FD21,1);
  dm.calc_dXijdri_dXjidri(1,1,r23,rsq23,del23,rho2,rho3,FD23,1);
  dm.calc_dXijdri_dXjidri(1,1,r24,rsq24,del24,rho2,rho4,FD24,1);

  dm.calc_dXijdri_dXjidri(1,1,r31,rsq31,del31,rho3,rho1,FD31,1);
  dm.calc_dXijdri_dXjidri(1,1,r32,rsq32,del32,rho3,rho2,FD32,1);
  dm.calc_dXijdri_dXjidri(1,1,r34,rsq34,del34,rho3,rho4,FD34,1);

  dm.calc_dXijdri_dXjidri(1,1,r41,rsq41,del41,rho4,rho1,FD41,1);
  dm.calc_dXijdri_dXjidri(1,1,r42,rsq42,del42,rho4,rho2,FD42,1);
  dm.calc_dXijdri_dXjidri(1,1,r43,rsq43,del43,rho4,rho3,FD43,1);

  REQUIRE(FD12.isApprox(fd12-fd21,1e-6));
  REQUIRE(FD13.isApprox(fd13-fd31,1e-6));
  REQUIRE(FD14.isApprox(fd14-fd41,1e-6));

  REQUIRE(FD21.isApprox(fd21-fd12,1e-6));
  REQUIRE(FD23.isApprox(fd23-fd32,1e-6));
  REQUIRE(FD24.isApprox(fd24-fd42,1e-6));

  REQUIRE(FD31.isApprox(fd31-fd13,1e-6));
  REQUIRE(FD32.isApprox(fd32-fd23,1e-6));
  REQUIRE(FD34.isApprox(fd34-fd43,1e-6));

  REQUIRE(FD41.isApprox(fd41-fd14,1e-6));
  REQUIRE(FD42.isApprox(fd42-fd24,1e-6));
  REQUIRE(FD43.isApprox(fd43-fd34,1e-6));

  // Calculate approximation
  fd_type fd12a(dm.size());
  fd_type fd13a(dm.size());
  fd_type fd14a(dm.size());
  fd12a.set_zero();
  fd13a.set_zero();
  fd14a.set_zero();
  double margin=1e-5;
  central_difference_mb(&dm,positions, 1, fd12a, 1e-8);
  central_difference_mb(&dm,positions, 2, fd13a, 1e-8);
  central_difference_mb(&dm,positions, 3, fd14a, 1e-8);

  std::vector<double> fx12(fd12[0].data(),fd12[0].data()+fd12[0].size());
  std::vector<double> fy12(fd12[1].data(),fd12[1].data()+fd12[1].size());
  std::vector<double> fz12(fd12[2].data(),fd12[2].data()+fd12[2].size());
  std::vector<double> fx12a(fd12a[0].data(),fd12a[0].data()+fd12a[0].size());
  std::vector<double> fy12a(fd12a[1].data(),fd12a[1].data()+fd12a[1].size());
  std::vector<double> fz12a(fd12a[2].data(),fd12a[2].data()+fd12a[2].size());
  REQUIRE_THAT(fx12a, Catch::Matchers::Approx(fx12).margin(margin));
  REQUIRE_THAT(fy12a, Catch::Matchers::Approx(fy12).margin(margin));
  REQUIRE_THAT(fz12a, Catch::Matchers::Approx(fz12).margin(margin));

  std::vector<double> fx13(fd13[0].data(),fd13[0].data()+fd13[0].size());
  std::vector<double> fy13(fd13[1].data(),fd13[1].data()+fd13[1].size());
  std::vector<double> fz13(fd13[2].data(),fd13[2].data()+fd13[2].size());
  std::vector<double> fx13a(fd13a[0].data(),fd13a[0].data()+fd13a[0].size());
  std::vector<double> fy13a(fd13a[1].data(),fd13a[1].data()+fd13a[1].size());
  std::vector<double> fz13a(fd13a[2].data(),fd13a[2].data()+fd13a[2].size());
  REQUIRE_THAT(fx13a, Catch::Matchers::Approx(fx13).margin(margin));
  REQUIRE_THAT(fy13a, Catch::Matchers::Approx(fy13).margin(margin));
  REQUIRE_THAT(fz13a, Catch::Matchers::Approx(fz13).margin(margin));

  std::vector<double> fx14(fd14[0].data(),fd14[0].data()+fd14[0].size());
  std::vector<double> fy14(fd14[1].data(),fd14[1].data()+fd14[1].size());
  std::vector<double> fz14(fd14[2].data(),fd14[2].data()+fd14[2].size());
  std::vector<double> fx14a(fd14a[0].data(),fd14a[0].data()+fd14a[0].size());
  std::vector<double> fy14a(fd14a[1].data(),fd14a[1].data()+fd14a[1].size());
  std::vector<double> fz14a(fd14a[2].data(),fd14a[2].data()+fd14a[2].size());
  REQUIRE_THAT(fx14a, Catch::Matchers::Approx(fx14).margin(margin));
  REQUIRE_THAT(fy14a, Catch::Matchers::Approx(fy14).margin(margin));
  REQUIRE_THAT(fz14a, Catch::Matchers::Approx(fz14).margin(margin));

  REQUIRE(dm.label() == label);
  REQUIRE(dm.size() == size);
}

// Fixture to initialize the PeriodicTable
struct TestFixture {
    TestFixture() {
        PeriodicTable::initialize();
    }
};
TEST_CASE_METHOD(TestFixture, "Testing DM_EAD 1", "[DM]") {
    Context context;
    context.add("RCUTMB", 10);
    context.remove("INITMB");
    context.add("INITMB", "true");
    context.add("TYPEMB", "DM_EAD");
    context.add("TYPEMB", 2);
    context.add("TYPEMB", 2);
    context.add("TYPEMB", 2);
    context.add("CGRIDMB", "1 2");
    context.add("SGRIDMB", "0.1 0.2");
    test_case<DM_EAD, Cut_Dummy>(context, "DM_EAD", 3 * 2);
    test_case<DM_EAD, Cut_Cos>(context, "DM_EAD", 3 * 2);
    test_case<DM_EAD, Cut_Tanh>(context, "DM_EAD", 3 * 2);
    test_case<DM_EAD, Cut_PolyS>(context, "DM_EAD", 3 * 2);
}

TEST_CASE_METHOD(TestFixture, "Testing DM_EAD 2", "[DM]") {
    Context context;
    context.add("RCUTMB", 6);
    context.remove("INITMB");
    context.add("INITMB", "true");
    context.add("TYPEMB", "DM_EAD");
    context.add("TYPEMB", 3);
    context.add("TYPEMB", 3);
    context.add("TYPEMB", 3);
    context.add("CGRIDMB", "1 2 3");
    context.add("SGRIDMB", "0.1 0.2 0.9");
    test_case<DM_EAD, Cut_Dummy>(context, "DM_EAD", 4 * 3);
    test_case<DM_EAD, Cut_Cos>(context, "DM_EAD", 4 * 3);
    test_case<DM_EAD, Cut_Tanh>(context, "DM_EAD", 4 * 3);
    test_case<DM_EAD, Cut_PolyS>(context, "DM_EAD", 4 * 3);
}

TEST_CASE_METHOD(TestFixture, "Testing DM_Blip 1", "[DM]") {
    Context context;
    context.add("RCUTMB", 10);
    context.remove("INITMB");
    context.add("INITMB", "true");
    context.add("TYPEMB", "DM_Blip");
    context.add("TYPEMB", 2);
    context.add("TYPEMB", 2);
    context.add("TYPEMB", 2);
    context.add("CGRIDMB", "1 2");
    context.add("SGRIDMB", "0.1 0.2");
    test_case<DM_Blip, Cut_Dummy>(context, "DM_Blip", 3 * 2);
    test_case<DM_Blip, Cut_Cos>(context, "DM_Blip", 3 * 2);
    test_case<DM_Blip, Cut_Tanh>(context, "DM_Blip", 3 * 2);
    test_case<DM_Blip, Cut_PolyS>(context, "DM_Blip", 3 * 2);
}

TEST_CASE_METHOD(TestFixture, "Testing DM_Blip 2", "[DM]") {
    Context context;
    context.add("RCUTMB", 6);
    context.remove("INITMB");
    context.add("INITMB", "true");
    context.add("TYPEMB", "DM_Blip");
    context.add("TYPEMB", 3);
    context.add("TYPEMB", 3);
    context.add("TYPEMB", 3);
    context.add("CGRIDMB", "1 2 3");
    context.add("SGRIDMB", "0.1 0.2 0.9");
    test_case<DM_Blip, Cut_Dummy>(context, "DM_Blip", 4 * 3);
    test_case<DM_Blip, Cut_Cos>(context, "DM_Blip", 4 * 3);
    test_case<DM_Blip, Cut_Tanh>(context, "DM_Blip", 4 * 3);
    test_case<DM_Blip, Cut_PolyS>(context, "DM_Blip", 4 * 3);
}

TEST_CASE_METHOD(TestFixture, "Testing DM_EAM 1", "[DM]") {
    Context context;
    context.add("RCUTMB", 5.58190742024);
    context.remove("INITMB");
    context.add("INITMB", "true");
    context.add("TYPEMB", "DM_EAM");
    context.add("TYPEMB", "tests_data/eam.alloy");
    test_case<DM_EAM, Cut_Dummy>(context, "DM_EAM", 1);

    // Throws
    Context c_err;
    c_err.remove("INITMB");
    c_err.add("INITMB", "true");
    c_err.add("TYPEMB", "DM_EAM");
    c_err.add("TYPEMB", "path/to/nowhere");
    REQUIRE_THROWS(DM_EAM(c_err));
}

TEST_CASE_METHOD(TestFixture, "Testing DM_EAM 2", "[DM]") {
    Context context;
    context.add("RCUTMB", 5.58190742024);
    context.remove("INITMB");
    context.add("INITMB", "true");
    context.add("TYPEMB", "DM_EAM");
    context.add("TYPEMB", "tests_data/eam.alloy");
    test_case<DM_EAM, Cut_Dummy>(context, "DM_EAM", 1);
}

TEST_CASE_METHOD(TestFixture, "Testing DM_EAM 3", "[DM]") {
    Context context;
    context.add("RCUTMB", 5.58190742024);
    context.remove("INITMB");
    context.add("INITMB", "true");
    context.add("TYPEMB", "DM_EAM");
    context.add("TYPEMB", "tests_data/eam.alloy");
    test_case<DM_EAM, Cut_Dummy>(context, "DM_EAM", 1);

    // Throws
    Context c_err;
    c_err.remove("INITMB");
    c_err.add("INITMB", "true");
    c_err.add("SETFL", "path/to/nowhere");
    REQUIRE_THROWS(DM_EAM(c_err));
}

TEST_CASE_METHOD(TestFixture, "Testing DM_Dummy", "[DM]") {
    Context c;
    c.add("RCUTMB", 5.58190742024);
    c.add("RCTYPEMB", "Cut_Cos");
    c.remove("INITMB");
    c.add("INITMB", "true");
    c.add("TYPEMB", "DM_Dummy");

    REQUIRE_NOTHROW(DM_Dummy());
    REQUIRE_NOTHROW(DM_Dummy(c));

    DM_Dummy d(c);
    REQUIRE(d.size() == 0);
    REQUIRE(d.label() == "DM_Dummy");
    aed_type aed;
    rho_type rho1;
    fd_type fd;
    Vec3d del;
    REQUIRE_NOTHROW(d.calc_aed(rho1, aed));
    REQUIRE_NOTHROW(d.init_rhoi(rho1));
    REQUIRE_NOTHROW(d.calc_rho(0,0, 0, 0, del, rho1));
    REQUIRE_NOTHROW(d.calc_dXijdri(1,1, 0, 0, del, rho1, fd));
    REQUIRE_NOTHROW(d.calc_dXijdri_dXjidri(1, 1, 0, 0, del, rho1, rho1, fd, 1));
    REQUIRE(d.rhoi_size() == 0);
    REQUIRE(d.rhoip_size() == 0);
}
