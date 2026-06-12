// test_d2_knot5.cpp
//
// One-sided quintic knot family (Knot5): primitive math audits (value /
// derivative FD parity, C^4 knot continuity, eta gauge identity), the
// D2_Knot5 descriptor (forces vs FD, closed-form dtheta parity vs
// FD-perturbed contexts, factory dispatch) and DM_REAM bcode 2.

#include "catch2/catch.hpp"

#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/core/periodic_table.h>
#include <tadah/core/registry.h>
#include <tadah/models/cutoffs.h>
#include <tadah/models/descriptors/d2/d2_knot5.h>
#include <tadah/models/descriptors/d2/d2_base.h>
#include <tadah/models/descriptors/dm/dm_ream.h>
#include <tadah/models/descriptors/d_basis_functions.h>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace {

using tadah::core::aed_type;
using tadah::core::fd_type;
using tadah::core::Context;
using tadah::core::PeriodicTable;
using tadah::models::Cut_Cos;
using tadah::models::D2_Knot5;
using tadah::models::DM_REAM;
using namespace tadah::models;  // K5L/K5R primitives

struct Knot5Fixture {
  Knot5Fixture() { PeriodicTable::initialize(); }
};

Context make_knot5_context(const std::vector<double>& cgrid, double rcut) {
  Context c;
  c.add("ATOM", "H");
  c.add("WATOM", 1);
  c.add("INIT2B", "true");
  c.add("VERBOSE", 0);
  {
    std::ostringstream ss;
    for (std::size_t i = 0; i < cgrid.size(); ++i) ss << (i ? " " : "") << "1.0";
    c.add("SGRID2B", ss.str());
  }
  {
    std::ostringstream ss;
    for (std::size_t i = 0; i < cgrid.size(); ++i)
      ss << (i ? " " : "") << cgrid[i];
    c.add("CGRID2B", ss.str());
  }
  c.add("RCUT2B", rcut);
  c.add("TYPE2B", "D2_Knot5");
  c.add("TYPE2B", 2);
  c.add("TYPE2B", 2);
  c.add("TYPE2B", "H");
  c.add("TYPE2B", "H");
  return c;
}

aed_type eval_aed(const std::vector<double>& cgrid, double rcut, double rij,
                  double sgrid_k = 1.0, int pert_idx = -1,
                  bool pert_cgrid = true, double h = 0.0) {
  // optionally perturb cgrid[pert_idx] (or sgrid) by h for FD references
  Context c;
  c.add("ATOM", "H");
  c.add("WATOM", 1);
  c.add("INIT2B", "true");
  c.add("VERBOSE", 0);
  {
    std::ostringstream ss;
    ss.precision(17);  // keep the 1e-6 FD perturbation (default 6 digits drops it)
    for (std::size_t i = 0; i < cgrid.size(); ++i) {
      double v = sgrid_k;
      if (!pert_cgrid && static_cast<int>(i) == pert_idx) v += h;
      ss << (i ? " " : "") << v;
    }
    c.add("SGRID2B", ss.str());
  }
  {
    std::ostringstream ss;
    ss.precision(17);
    for (std::size_t i = 0; i < cgrid.size(); ++i) {
      double v = cgrid[i];
      if (pert_cgrid && static_cast<int>(i) == pert_idx) v += h;
      ss << (i ? " " : "") << v;
    }
    c.add("CGRID2B", ss.str());
  }
  c.add("RCUT2B", rcut);
  c.add("TYPE2B", "D2_Knot5");
  c.add("TYPE2B", 2);
  c.add("TYPE2B", 2);
  c.add("TYPE2B", "H");
  c.add("TYPE2B", "H");
  D2_Knot5 d(c);
  Cut_Cos fc(rcut);
  d.set_fcut(&fc, false);
  aed_type aed(d.size());
  aed.set_zero();
  d.calc_aed(1, 1, rij, rij * rij, aed);
  return aed;
}

}  // namespace

TEST_CASE_METHOD(Knot5Fixture, "K5 primitives: values, FD parity, gauge",
                 "[knot5]") {
  // value identities
  REQUIRE(K5L(-1.0) == Approx(1.0));            // (t-r)=1
  REQUIRE(K5L(0.0) == 0.0);                     // at the knot
  REQUIRE(K5L(0.5) == 0.0);                     // beyond (right of) the knot
  REQUIRE(K5R(1.0) == Approx(1.0));
  REQUIRE(K5R(-0.5) == 0.0);

  // derivative convention: dK5(x_c, eta) == eta * d/dx_c K5(x_c)
  const double h = 1e-6;
  for (double x : {-1.7, -0.9, -0.3, 0.2, 0.8, 1.4}) {
    for (double eta : {1.0, 2.5}) {
      const double fdL = eta * (K5L(x + h) - K5L(x - h)) / (2 * h);
      REQUIRE(dK5L(x, eta) == Approx(fdL).margin(1e-7));
      const double fdR = eta * (K5R(x + h) - K5R(x - h)) / (2 * h);
      REQUIRE(dK5R(x, eta) == Approx(fdR).margin(1e-7));
      const double fd2L = (dK5L(x + h, 1.0) - dK5L(x - h, 1.0)) / (2 * h);
      REQUIRE(d2K5L(x) == Approx(fd2L).margin(1e-5));
      const double fd2R = (dK5R(x + h, 1.0) - dK5R(x - h, 1.0)) / (2 * h);
      REQUIRE(d2K5R(x) == Approx(fd2R).margin(1e-5));
    }
  }

  // C^4 at the knot: value..4th derivative vanish from the active side
  const double e = 1e-3;
  REQUIRE(K5L(-e) == Approx(0.0).margin(1e-14));      // ~e^5
  REQUIRE(dK5L(-e, 1.0) == Approx(0.0).margin(1e-11)); // ~5e^4
  REQUIRE(d2K5L(-e) == Approx(0.0).margin(3e-8));      // 20e^3 = 2e-8

  // eta gauge: c * K5L(eta*(x-mu)) == (c*eta^5) * K5L(x-mu)
  const double xm = -0.37, eta = 1.9;
  REQUIRE(K5L(eta * xm) ==
          Approx(std::pow(eta, 5) * K5L(xm)).epsilon(1e-12));
  REQUIRE(K5R(eta * (-xm)) ==
          Approx(std::pow(eta, 5) * K5R(-xm)).epsilon(1e-12));
}

TEST_CASE_METHOD(Knot5Fixture, "D2_Knot5: factory dispatch + value + force FD",
                 "[knot5]") {
  const double rcut = 6.0;
  const std::vector<double> cgrid = {2.5, 3.5, 4.5, 5.5};
  Context c = make_knot5_context(cgrid, rcut);

  // string-factory dispatch (the path DC_Selector / LAMMPS uses)
  using tadah::core::factory;
  auto* d_ptr = factory<tadah::models::D2_Base, const Context&>("D2_Knot5", c);
  REQUIRE(d_ptr != nullptr);
  REQUIRE(d_ptr->label() == "D2_Knot5");
  delete d_ptr;

  D2_Knot5 d(c);
  Cut_Cos fc(rcut);
  d.set_fcut(&fc, false);

  const double rij = 3.0;
  aed_type aed(d.size());
  aed.set_zero();
  d.calc_aed(1, 1, rij, rij * rij, aed);
  // hand-computed: (t - r)_+^5 * fc(r) per knot
  for (std::size_t k = 0; k < cgrid.size(); ++k) {
    const double expect = K5L(1.0 * (rij - cgrid[k])) * fc.calc(rij);
    REQUIRE(aed(k) == Approx(expect).margin(1e-12));
  }

  // radial force column vs FD of calc_aed in rij
  fd_type fd(d.size());
  fd.set_zero();
  d.calc_dXijdri(1, 1, rij, rij * rij, fd);
  const double h = 1e-6;
  aed_type ap(d.size()), am(d.size());
  ap.set_zero(); am.set_zero();
  d.calc_aed(1, 1, rij + h, (rij + h) * (rij + h), ap);
  d.calc_aed(1, 1, rij - h, (rij - h) * (rij - h), am);
  for (std::size_t k = 0; k < cgrid.size(); ++k)
    REQUIRE(fd(k, 0) == Approx((ap(k) - am(k)) / (2 * h)).margin(1e-6));

  // calc_all consistency with calc_aed + calc_dXijdri
  aed_type aed2(d.size()); aed2.set_zero();
  fd_type fd2(d.size()); fd2.set_zero();
  d.calc_all(1, 1, rij, rij * rij, aed2, fd2);
  for (std::size_t k = 0; k < cgrid.size(); ++k) {
    REQUIRE(aed2(k) == Approx(aed(k)));
    REQUIRE(fd2(k, 0) == Approx(fd(k, 0)));
  }
}

TEST_CASE_METHOD(Knot5Fixture, "D2_Knot5: dtheta parity vs FD contexts",
                 "[knot5]") {
  const double rcut = 6.0;
  const std::vector<double> cgrid = {2.5, 3.5, 4.5};
  const double rij = 3.1;
  Context c = make_knot5_context(cgrid, rcut);
  D2_Knot5 d(c);
  Cut_Cos fc(rcut);
  d.set_fcut(&fc, false);

  const double h = 1e-6;
  for (int k = 1; k <= static_cast<int>(cgrid.size()); ++k) {
    // CGRID2B
    aed_type dD(d.size());
    dD.set_zero();
    REQUIRE(d.calc_dDdtheta("CGRID2B", k, 1, 1, rij, rij * rij, dD));
    aed_type p = eval_aed(cgrid, rcut, rij, 1.0, k - 1, true, h);
    aed_type m = eval_aed(cgrid, rcut, rij, 1.0, k - 1, true, -h);
    REQUIRE(dD(k - 1) ==
            Approx((p(k - 1) - m(k - 1)) / (2 * h)).margin(1e-6));
    // SGRID2B
    dD.set_zero();
    REQUIRE(d.calc_dDdtheta("SGRID2B", k, 1, 1, rij, rij * rij, dD));
    p = eval_aed(cgrid, rcut, rij, 1.0, k - 1, false, h);
    m = eval_aed(cgrid, rcut, rij, 1.0, k - 1, false, -h);
    REQUIRE(dD(k - 1) ==
            Approx((p(k - 1) - m(k - 1)) / (2 * h)).margin(1e-6));
  }
}

TEST_CASE_METHOD(Knot5Fixture, "DM_REAM Knot5 (bcode 2): basis orientation",
                 "[knot5]") {
  Context c;
  c.add("VERBOSE", 0);
  c.add("INITMB", "true");
  c.add("RCUTMB", 6.0);
  c.add("RCTYPEMB", "Cut_Cos");
  c.add("TYPEMB", "DM_REAM");
  c.add("TYPEMB", 2);          // Knot5
  c.add("TYPEMB", "H");
  c.add("TYPEMB", "H");
  c.add("CGRIDRHOMB", "3.0 4.5 6.0");
  c.add("SGRIDRHOMB", "1.0 1.0 1.0");
  c.add("AGRIDMB",    "0.40 0.25 0.12");
  c.add("CGRIDFMB", "0.0 0.5 1.0");
  c.add("SGRIDFMB", "1.0 1.0 1.0");
  DM_REAM dm(c);
  Cut_Cos fc(6.0);
  dm.set_fcut(&fc, false);

  // density: LEFT-active radial psi_n = (mu_n - r)_+^5 * fc
  const double r = 3.5;
  tadah::core::rho_type rhoi;
  dm.init_rhoi(rhoi);
  tadah::core::Vec3d del(r, 0.0, 0.0);
  dm.calc_rho_mb(1, 1, del, rhoi);
  const double expect_rho = 0.40 * K5L(r - 3.0, fc.calc(r)) +
                            0.25 * K5L(r - 4.5, fc.calc(r)) +
                            0.12 * K5L(r - 6.0, fc.calc(r));
  REQUIRE(rhoi(0) == Approx(expect_rho).margin(1e-12));
  REQUIRE(K5L(r - 3.0) == 0.0);  // knot below r: contributes nothing (left-active)
  REQUIRE(K5L(r - 4.5) > 0.0);   // knot above r: active

  // embedding: RIGHT-active g_m = (rho - mu_m)_+^5
  aed_type aed(dm.size());
  aed.set_zero();
  dm.calc_aed(rhoi, aed);
  const double rho_i = rhoi(0);
  REQUIRE(aed(0) == Approx(K5R(rho_i - 0.0)).margin(1e-12));
  REQUIRE(aed(1) == Approx(K5R(rho_i - 0.5)).margin(1e-12));
  REQUIRE(aed(2) == Approx(K5R(rho_i - 1.0)).margin(1e-12));
}


TEST_CASE_METHOD(Knot5Fixture, "D2_Knot5 mixed degrees (DGRID2B)", "[knot5][knotmix]") {
  const double rcut = 6.0;
  Context c;
  c.add("ATOM", "H");
  c.add("WATOM", 1);
  c.add("INIT2B", "true");
  c.add("VERBOSE", 0);
  c.add("CGRID2B", "2.5 3.5 4.5");
  c.add("SGRID2B", "1.0 1.0 1.0");
  c.add("DGRID2B", "5 3 5");          // middle atom CUBIC
  c.add("RCUT2B", rcut);
  c.add("TYPE2B", "D2_Knot5");
  c.add("TYPE2B", 2);
  c.add("TYPE2B", 2);
  c.add("TYPE2B", "H");
  c.add("TYPE2B", "H");
  D2_Knot5 d(c);
  Cut_Cos fc(rcut);
  d.set_fcut(&fc, false);

  const double rij = 3.0;
  aed_type aed(d.size());
  aed.set_zero();
  d.calc_aed(1, 1, rij, rij * rij, aed);
  REQUIRE(aed(0) == Approx(K5L(rij - 2.5, fc.calc(rij))).margin(1e-12));
  REQUIRE(aed(1) == Approx(K3L(rij - 3.5, fc.calc(rij))).margin(1e-12));  // cubic
  REQUIRE(aed(2) == Approx(K5L(rij - 4.5, fc.calc(rij))).margin(1e-12));

  // force FD parity with mixed degrees
  fd_type fd(d.size());
  fd.set_zero();
  d.calc_dXijdri(1, 1, rij, rij * rij, fd);
  const double h = 1e-6;
  aed_type ap(d.size()), am(d.size());
  ap.set_zero(); am.set_zero();
  d.calc_aed(1, 1, rij + h, (rij + h) * (rij + h), ap);
  d.calc_aed(1, 1, rij - h, (rij - h) * (rij - h), am);
  for (std::size_t k = 0; k < 3; ++k)
    REQUIRE(fd(k, 0) == Approx((ap(k) - am(k)) / (2 * h)).margin(1e-6));

  // invalid degree rejected
  Context bad;
  bad.add("ATOM", "H");
  bad.add("WATOM", 1);
  bad.add("INIT2B", "true");
  bad.add("VERBOSE", 0);
  bad.add("CGRID2B", "2.5 3.5");
  bad.add("SGRID2B", "1.0 1.0");
  bad.add("DGRID2B", "5 4");          // 4 is not a knot degree
  bad.add("RCUT2B", rcut);
  bad.add("TYPE2B", "D2_Knot5");
  bad.add("TYPE2B", 2);
  bad.add("TYPE2B", 2);
  bad.add("TYPE2B", "H");
  bad.add("TYPE2B", "H");
  REQUIRE_THROWS_WITH(D2_Knot5(bad),
                      Catch::Matchers::Contains("must be 3 or 5"));
}

TEST_CASE_METHOD(Knot5Fixture, "DM_REAM mixed degrees (DGRID*)", "[knot5][knotmix]") {
  Context c;
  c.add("VERBOSE", 0);
  c.add("INITMB", "true");
  c.add("RCUTMB", 6.0);
  c.add("RCTYPEMB", "Cut_Cos");
  c.add("TYPEMB", "DM_REAM");
  c.add("TYPEMB", 2);
  c.add("TYPEMB", "H");
  c.add("TYPEMB", "H");
  c.add("CGRIDRHOMB", "3.0 4.5 6.0");
  c.add("SGRIDRHOMB", "1.0 1.0 1.0");
  c.add("AGRIDMB",    "0.40 0.25 0.12");
  c.add("DGRIDRHOMB", "5 3 5");
  c.add("CGRIDFMB", "0.0 0.5 1.0");
  c.add("SGRIDFMB", "1.0 1.0 1.0");
  c.add("DGRIDFMB", "3 5 5");
  DM_REAM dm(c);
  Cut_Cos fc(6.0);
  dm.set_fcut(&fc, false);

  const double r = 3.5;
  tadah::core::rho_type rhoi;
  dm.init_rhoi(rhoi);
  tadah::core::Vec3d del(r, 0.0, 0.0);
  dm.calc_rho_mb(1, 1, del, rhoi);
  const double expect_rho = 0.40 * K5L(r - 3.0, fc.calc(r)) +
                            0.25 * K3L(r - 4.5, fc.calc(r)) +   // cubic
                            0.12 * K5L(r - 6.0, fc.calc(r));
  REQUIRE(rhoi(0) == Approx(expect_rho).margin(1e-12));

  aed_type aed(dm.size());
  aed.set_zero();
  dm.calc_aed(rhoi, aed);
  const double rho_i = rhoi(0);
  REQUIRE(aed(0) == Approx(K3R(rho_i - 0.0)).margin(1e-12));    // cubic
  REQUIRE(aed(1) == Approx(K5R(rho_i - 0.5)).margin(1e-12));
  REQUIRE(aed(2) == Approx(K5R(rho_i - 1.0)).margin(1e-12));
}
