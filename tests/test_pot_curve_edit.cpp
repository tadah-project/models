#include "catch2/catch.hpp"
#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/core/logger.h>
#include <tadah/core/number_format.h>
#include <tadah/core/periodic_table.h>
#include <tadah/models/eam/pot_curve_edit.h>
#include <tadah/models/descriptors/basis_dispatch.h>
#include <tadah/models/descriptors/d2/d2_knot5.h>
#include <tadah/models/descriptors/dm/dm_ream.h>
#include <tadah/models/descriptors/dm/embf0.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace tadah::models;
using namespace tadah::models::potedit;
using tadah::core::Context;

namespace {

struct EditFixture {
  EditFixture() { tadah::core::PeriodicTable::initialize(); }
};

void add_common(Context &c, const std::string &el) {
  c.add("MODEL", "M_BLR");
  c.add("MODEL", "BF_Linear");
  c.add("ATOM", el);
  c.add("WATOM", 1.0);
  c.add("BIAS", "false");
  c.add("NORM", "false");
  c.add("FORCE", "false");
  c.add("STRESS", "false");
  c.add("LAMBDA", 0);
  c.add("OALGO", 2);
}

// Knot-family pot mirroring the `tadah refit` output shape (deliberately
// UNSORTED grids — layouts are positional).
Context make_knot_pot() {
  Context c;
  add_common(c, "Ta");
  c.add("INIT2B", "true");
  c.add("TYPE2B", "D2_Knot5");
  c.add("TYPE2B", 3);
  c.add("TYPE2B", 0);
  c.add("TYPE2B", "Ta");
  c.add("TYPE2B", "Ta");
  for (double v : {5.0, 3.0, 4.2}) c.add("CGRID2B", v);
  c.add("RCUT2B", 5.0);  // == max knot (self-consistent)
  c.add("RCTYPE2B", "Cut_Dummy");
  c.add("INITMB", "true");
  c.add("TYPEMB", "DM_REAM");
  c.add("TYPEMB", 2);  // bcode Knot5
  c.add("TYPEMB", 2);  // nC
  c.add("TYPEMB", 0);  // nS
  c.add("TYPEMB", 2);  // nAmp
  c.add("TYPEMB", 3);  // nEc
  c.add("TYPEMB", 0);  // nEs
  c.add("TYPEMB", 0);  // nKd
  c.add("TYPEMB", 0);  // nKe
  c.add("TYPEMB", "Ta");
  c.add("TYPEMB", "Ta");
  for (double v : {4.8, 2.5}) c.add("CGRIDMB", v);
  for (double v : {0.01, -1.2}) c.add("AMPGRIDMB", v);
  for (double v : {0.9, -2.9, 0.1}) c.add("CEMBFUNC", v);
  c.add("RCUTMB", 4.8);
  c.add("RCTYPEMB", "Cut_Dummy");
  for (double v : {1.5, -900.0, 12.0, 220.0, -160.0, 5300.0})
    c.add("WEIGHTS", v);
  c.add("FIXINDEX", "1-6");
  for (double v : {1.5, -900.0, 12.0, 220.0, -160.0, 5300.0})
    c.add("FIXWEIGHT", v);
  return c;
}

// Blip-family pot (bcode 0, width keys everywhere).
Context make_blip_pot() {
  Context c;
  add_common(c, "Cu");
  c.add("INIT2B", "true");
  c.add("TYPE2B", "D2_Blip");
  c.add("TYPE2B", 3);
  c.add("TYPE2B", 3);
  c.add("TYPE2B", "Cu");
  c.add("TYPE2B", "Cu");
  for (double v : {5.0, 3.0, 4.2}) c.add("CGRID2B", v);
  for (double v : {1.2, 1.0, 1.4}) c.add("SGRID2B", v);
  c.add("RCUT2B", 5.0);
  c.add("RCTYPE2B", "Cut_Dummy");
  c.add("INITMB", "true");
  c.add("TYPEMB", "DM_REAM");
  c.add("TYPEMB", 0);  // bcode Blip
  c.add("TYPEMB", 2);
  c.add("TYPEMB", 2);
  c.add("TYPEMB", 2);
  c.add("TYPEMB", 3);
  c.add("TYPEMB", 3);
  c.add("TYPEMB", 0);
  c.add("TYPEMB", 0);
  c.add("TYPEMB", "Cu");
  c.add("TYPEMB", "Cu");
  for (double v : {4.8, 2.5}) c.add("CGRIDMB", v);
  for (double v : {1.0, 1.1}) c.add("SGRIDMB", v);
  for (double v : {0.01, -1.2}) c.add("AMPGRIDMB", v);
  for (double v : {0.9, 2.9, 0.1}) c.add("CEMBFUNC", v);
  for (double v : {1.5, 2.0, 1.8}) c.add("SEMBFUNC", v);
  c.add("RCUTMB", 4.8);
  c.add("RCTYPEMB", "Cut_Dummy");
  for (double v : {1.5, -900.0, 12.0, 220.0, -160.0, 5300.0})
    c.add("WEIGHTS", v);
  c.add("FIXINDEX", "1-6");
  for (double v : {1.5, -900.0, 12.0, 220.0, -160.0, 5300.0})
    c.add("FIXWEIGHT", v);
  return c;
}

// Pure two-body knot pot (the HPO_two_body_fit seed shape: TYPE2B only, no
// TYPEMB/INITMB anywhere, WEIGHTS = one per pair basis function).
Context make_pure2b_pot() {
  Context c;
  add_common(c, "Kr");
  c.add("INIT2B", "true");
  c.add("TYPE2B", "D2_Knot5");
  c.add("TYPE2B", 3);
  c.add("TYPE2B", 0);
  c.add("TYPE2B", "Kr");
  c.add("TYPE2B", "Kr");
  for (double v : {5.0, 3.0, 4.2}) c.add("CGRID2B", v);
  c.add("RCUT2B", 5.0);
  c.add("RCTYPE2B", "Cut_Dummy");
  for (double v : {1.5, -900.0, 12.0}) c.add("WEIGHTS", v);
  c.add("FIXINDEX", "1-3");
  for (double v : {1.5, -900.0, 12.0}) c.add("FIXWEIGHT", v);
  return c;
}

// Knot pair channel + a FOREIGN (non-DM_REAM) many-body descriptor: the
// pair block is the WEIGHTS prefix; the many-body tail must ride along
// opaquely and every MB key must survive untouched.
Context make_foreign_mb_pot() {
  Context c;
  add_common(c, "Ta");
  c.add("INIT2B", "true");
  c.add("TYPE2B", "D2_Knot5");
  c.add("TYPE2B", 3);
  c.add("TYPE2B", 0);
  c.add("TYPE2B", "Ta");
  c.add("TYPE2B", "Ta");
  for (double v : {5.0, 3.0, 4.2}) c.add("CGRID2B", v);
  c.add("RCUT2B", 5.0);
  c.add("RCTYPE2B", "Cut_Dummy");
  c.add("INITMB", "true");
  c.add("TYPEMB", "DM_Blip");
  c.add("TYPEMB", 2);
  c.add("TYPEMB", 2);
  c.add("TYPEMB", "Ta");
  c.add("TYPEMB", "Ta");
  c.add("CGRIDMB", "2.5 4.0");
  c.add("SGRIDMB", "1.0 1.2");
  c.add("RCUTMB", 4.6);
  c.add("RCTYPEMB", "Cut_Cos");
  // WEIGHTS = [3 pair | 2 many-body]
  for (double v : {1.5, -900.0, 12.0, 7.0, -3.0}) c.add("WEIGHTS", v);
  c.add("FIXINDEX", "1-5");
  for (double v : {1.5, -900.0, 12.0, 7.0, -3.0}) c.add("FIXWEIGHT", v);
  return c;
}

// Knot pot with per-function embedding degrees (the -2.9 knot is cubic).
Context make_knot3_pot() {
  Context c = make_knot_pot();
  c.remove("TYPEMB");
  c.add("TYPEMB", "DM_REAM");
  for (int v : {2, 2, 0, 2, 3, 0, 0, 3}) c.add("TYPEMB", v);  // nKe = 3
  c.add("TYPEMB", "Ta");
  c.add("TYPEMB", "Ta");
  for (double v : {5.0, 3.0, 5.0}) c.add("KDEGEMB", v);
  return c;
}

Context with_embf0(Context c, const std::string &mode) {
  c.add("EMBF0", mode);
  return c;
}

std::string stored_embf0_text(const Context &c) {
  return c.stored_value("EMBF0").value_or("");
}

// The embedding energy the runtime evaluates: w . calc_aed of a DM_REAM
// built from a (rebuilt) potential Context, at each density in rhos, with w
// the F coefficients (summed in order).
std::vector<double> ream_energies(Context pot, const std::vector<double> &w,
                                  const std::vector<double> &rhos) {
  pot.add("VERBOSE", 0);
  DM_REAM dm(pot);
  REQUIRE(dm.size() == w.size());
  std::vector<double> E;
  for (double rho : rhos) {
    tadah::core::rho_type rhoi;
    dm.init_rhoi(rhoi);
    rhoi(0) = rho;
    tadah::core::aed_type aed(dm.size());
    aed.set_zero();
    dm.calc_aed(rhoi, aed);
    double e = 0.0;
    for (std::size_t m = 0; m < dm.size(); ++m) e += w[m] * aed(m);
    E.push_back(e);
  }
  return E;
}

// Max |difference| of value/deriv/deriv2 between two curve views over a
// dense sample of [lo, hi].
double curve_maxdiff(const CurveView &a, const CurveView &b, double lo,
                     double hi, int N = 2000) {
  double mx = 0.0;
  for (int i = 0; i < N; ++i) {
    const double x = lo + (hi - lo) * i / (N - 1);
    mx = std::max(mx, std::fabs(a.value(x) - b.value(x)));
    mx = std::max(mx, std::fabs(a.deriv(x) - b.deriv(x)));
    mx = std::max(mx, std::fabs(a.deriv2(x) - b.deriv2(x)));
  }
  return mx;
}

}  // namespace

TEST_CASE_METHOD(EditFixture, "potedit: parse->rebuild round-trip (no edit)",
                 "[potedit]") {
  for (Context pot : {make_knot_pot(), make_blip_pot()}) {
    PotEAMView v = parse_pot_view(pot);
    Context out = rebuild_pot_context(v);
    for (const char *key :
         {"CGRID2B", "SGRID2B", "CGRIDMB", "SGRIDMB", "AMPGRIDMB", "CEMBFUNC",
          "SEMBFUNC", "WEIGHTS", "FIXWEIGHT", "TYPE2B", "TYPEMB"}) {
      INFO("key=" << key);
      REQUIRE(out.contains(key) == pot.contains(key));
      if (!pot.contains(key)) continue;
      REQUIRE(out.size(key) == pot.size(key));
      for (std::size_t i = 0; i < pot.size(key); ++i)
        REQUIRE(out.get<std::string>(key, i) == pot.get<std::string>(key, i));
    }
    REQUIRE(out.get<double>("RCUT2B") == pot.get<double>("RCUT2B"));
    REQUIRE(out.get<double>("RCUTMB") == pot.get<double>("RCUTMB"));
    REQUIRE(out.get<std::string>("FIXINDEX") ==
            pot.get<std::string>("FIXINDEX"));
  }
}

TEST_CASE_METHOD(EditFixture, "potedit: parse rejects non-refit shapes",
                 "[potedit]") {
  Context c = make_knot_pot();
  c.remove("TYPE2B");
  c.add("TYPE2B", "D2_BP");
  c.add("TYPE2B", 3);
  c.add("TYPE2B", 3);
  c.add("TYPE2B", "Ta");
  c.add("TYPE2B", "Ta");
  REQUIRE_THROWS_WITH(parse_pot_view(c),
                      Catch::Contains("unsupported two-body"));

  Context m = make_knot_pot();
  Context bad;  // no TYPE2B/TYPEMB at all
  REQUIRE_THROWS_WITH(parse_pot_view(bad), Catch::Contains("TYPE2B"));

  // WEIGHTS size mismatch
  Context w = make_knot_pot();
  w.add("WEIGHTS", 1.0);
  REQUIRE_THROWS_WITH(parse_pot_view(w), Catch::Contains("WEIGHTS"));
}

TEST_CASE_METHOD(EditFixture,
                 "potedit: ADD is bit-identical on every channel and family",
                 "[potedit]") {
  struct Case {
    const char *tag;
    bool knot;
    Channel ch;
    double centre;
  };
  // Blip centres sit inward enough that the default eta keeps the new
  // support inside RCUT (support crossing RCUT is an error since the
  // bump-support gate).
  for (const Case &cs : {Case{"knot-phi", true, Channel::Phi, 3.6},
                         Case{"knot-rho", true, Channel::Rho, 3.9},
                         Case{"knot-F", true, Channel::F, 0.4},
                         Case{"blip-phi", false, Channel::Phi, 3.3},
                         Case{"blip-rho", false, Channel::Rho, 2.6},
                         Case{"blip-F", false, Channel::F, 1.4}}) {
    INFO("case=" << cs.tag);
    Context pot = cs.knot ? make_knot_pot() : make_blip_pot();
    PotEAMView v = parse_pot_view(pot);
    CurveView &cv = cs.ch == Channel::Phi ? v.phi
                  : cs.ch == Channel::Rho ? v.rho
                                          : v.F;
    const CurveView before = cv;
    const bool radial = cs.ch != Channel::F;
    AddReport rep = add_basis(cv, AddSpec{cs.centre, {}, {}}, radial);
    REQUIRE(rep.index == before.n());
    REQUIRE(cv.coeff.back() == 0.0);
    // bit-identical curve, derivatives INCLUDING beyond rcut / max knot
    REQUIRE(curve_maxdiff(before, cv, 0.2, 8.0) == 0.0);

    // survives the full rebuild -> parse round-trip
    Context out = rebuild_pot_context(v);
    PotEAMView v2 = parse_pot_view(out);
    const CurveView &cv2 = cs.ch == Channel::Phi ? v2.phi
                         : cs.ch == Channel::Rho ? v2.rho
                                                 : v2.F;
    REQUIRE(curve_maxdiff(before, cv2, 0.2, 8.0) == 0.0);
  }
}

TEST_CASE_METHOD(EditFixture, "potedit: ADD bookkeeping and validation",
                 "[potedit]") {
  SECTION("knot deg 3 materialises the degree array") {
    PotEAMView v = parse_pot_view(make_knot_pot());
    AddReport rep = add_basis(v.phi, AddSpec{3.6, {}, 3}, true);
    REQUIRE(rep.deg_used == 3);
    REQUIRE(v.phi.deg == std::vector<int>{5, 5, 5, 3});
    Context out = rebuild_pot_context(v);
    REQUIRE(out.size("KDEG2B") == 4);
    REQUIRE(out.get<int>("TYPE2B", 1) == 4);  // nC
    REQUIRE(out.get<int>("TYPE2B", 2) == 4);  // nKdeg
    REQUIRE(out.size("WEIGHTS") == 7);
    REQUIRE(out.get<std::string>("FIXINDEX") == "1-7");
    REQUIRE(out.get<double>("WEIGHTS", 3) == 0.0);   // appended phi coeff
    REQUIRE(out.get<double>("WEIGHTS", 4) == 220.0); // F block follows
  }
  SECTION("knot deg 5 keeps the degree array absent") {
    PotEAMView v = parse_pot_view(make_knot_pot());
    add_basis(v.phi, AddSpec{3.6, {}, {}}, true);
    REQUIRE(v.phi.deg.empty());
    Context out = rebuild_pot_context(v);
    REQUIRE(!out.contains("KDEG2B"));
  }
  SECTION("blip eta defaults to the mean of the two nearest widths") {
    PotEAMView v = parse_pot_view(make_blip_pot());
    AddReport rep = add_basis(v.phi, AddSpec{3.2, {}, {}}, true);
    REQUIRE(rep.eta_used == Approx(0.5 * (1.0 + 1.4)));  // mu 3.0 and 4.2
    REQUIRE(v.phi.eta.back() == Approx(1.2));
  }
  SECTION("blip eta override") {
    PotEAMView v = parse_pot_view(make_blip_pot());
    AddReport rep = add_basis(v.phi, AddSpec{3.55, 2.5, {}}, true);
    REQUIRE(rep.eta_used == 2.5);
  }
  SECTION("knot beyond the largest knot grows the derived cutoff") {
    PotEAMView v = parse_pot_view(make_knot_pot());
    AddReport rep = add_basis(v.phi, AddSpec{5.6, {}, {}}, true);
    REQUIRE_FALSE(rep.warnings.empty());
    Context out = rebuild_pot_context(v);
    REQUIRE(out.get<double>("RCUT2B") == Approx(5.6));
  }
  SECTION("blip support entirely beyond RCUT is rejected") {
    PotEAMView v = parse_pot_view(make_blip_pot());
    REQUIRE_THROWS_WITH(add_basis(v.phi, AddSpec{6.5, 2.0, {}}, true),
                        Catch::Contains("entirely beyond"));
  }
  SECTION("blip support crossing RCUT errors under the hard cutoff") {
    PotEAMView v = parse_pot_view(make_blip_pot());
    REQUIRE_THROWS_WITH(add_basis(v.phi, AddSpec{4.8, 2.0, {}}, true),
                        Catch::Contains("crosses RCUT"));
  }
  SECTION("blip support crossing RCUT warns under a smooth cutoff") {
    PotEAMView v = parse_pot_view(make_blip_pot());
    AddReport rep = add_basis(v.phi, AddSpec{4.8, 2.0, {}}, true,
                              /*hard_cutoff=*/false);
    REQUIRE_FALSE(rep.warnings.empty());
    Context out = rebuild_pot_context(v);
    REQUIRE(out.get<double>("RCUT2B") == 5.0);  // RCUT never auto-changed
  }
  SECTION("blip support ending exactly at RCUT passes silently") {
    PotEAMView v = parse_pot_view(make_blip_pot());
    // centre + w/eta = 4.0 + 2/2.0 = 5.0 == RCUT2B: C^2-smooth at the edge.
    AddReport rep = add_basis(v.phi, AddSpec{4.0, 2.0, {}}, true);
    REQUIRE(rep.warnings.empty());
  }
  SECTION("invalid specs throw") {
    PotEAMView vk = parse_pot_view(make_knot_pot());
    REQUIRE_THROWS_WITH(add_basis(vk.phi, AddSpec{3.6, 1.5, {}}, true),
                        Catch::Contains("pure gauge"));
    REQUIRE_THROWS_WITH(add_basis(vk.phi, AddSpec{-1.0, {}, {}}, true),
                        Catch::Contains("positive"));
    REQUIRE_THROWS_WITH(add_basis(vk.phi, AddSpec{3.6, {}, 4}, true),
                        Catch::Contains("3 or 5"));
    PotEAMView vb = parse_pot_view(make_blip_pot());
    REQUIRE_THROWS_WITH(add_basis(vb.phi, AddSpec{3.6, {}, 3}, true),
                        Catch::Contains("knot family"));
  }
}

TEST_CASE_METHOD(EditFixture,
                 "potedit: bump_channel_halfwidth and the cutoff-kind probe",
                 "[potedit][rcut]") {
  SECTION("half-width table") {
    Context blip = make_blip_pot();  // D2_Blip + DM_REAM bcode 0
    REQUIRE(bump_channel_halfwidth(blip, "2B") == 2.0);
    REQUIRE(bump_channel_halfwidth(blip, "MB") == 2.0);

    Context blip5 = make_blip_pot();
    blip5.remove("TYPE2B");
    blip5.add("TYPE2B", "D2_Blip5");
    blip5.add("TYPE2B", 3);
    blip5.add("TYPE2B", 3);
    blip5.add("TYPE2B", "Cu");
    blip5.add("TYPE2B", "Cu");
    REQUIRE(bump_channel_halfwidth(blip5, "2B") == 3.0);

    Context mb5 = make_blip_pot();
    mb5.remove("TYPEMB");
    mb5.add("TYPEMB", "DM_REAM");
    mb5.add("TYPEMB", 3);  // bcode Blip5
    for (int v : {2, 2, 2, 3, 3, 0, 0}) mb5.add("TYPEMB", v);
    mb5.add("TYPEMB", "Cu");
    mb5.add("TYPEMB", "Cu");
    REQUIRE(bump_channel_halfwidth(mb5, "MB") == 3.0);

    Context knot = make_knot_pot();  // D2_Knot5 + DM_REAM bcode 2
    REQUIRE_FALSE(bump_channel_halfwidth(knot, "2B").has_value());
    REQUIRE_FALSE(bump_channel_halfwidth(knot, "MB").has_value());
    Context empty;
    REQUIRE_FALSE(bump_channel_halfwidth(empty, "2B").has_value());
  }
  SECTION("cutoff-kind probe") {
    Context blip = make_blip_pot();
    REQUIRE(channel_cutoff_is_discontinuous(blip, "2B"));  // Cut_Dummy
    Context foreign = make_foreign_mb_pot();
    REQUIRE_FALSE(channel_cutoff_is_discontinuous(foreign, "MB"));  // Cut_Cos
    Context absent = make_blip_pot();
    absent.remove("RCTYPE2B");  // descriptor ctor default is Cut_Dummy
    REQUIRE(channel_cutoff_is_discontinuous(absent, "2B"));
    Context smooth = make_blip_pot();
    smooth.remove("RCTYPE2B");
    smooth.add("RCTYPE2B", "Cut_Poly");
    REQUIRE_FALSE(channel_cutoff_is_discontinuous(smooth, "2B"));
  }
}

TEST_CASE_METHOD(EditFixture, "potedit: REMOVE selection and compensation",
                 "[potedit]") {
  SECTION("nearest-centre selection, compensation beats plain removal") {
    PotEAMView v = parse_pot_view(make_knot_pot());
    const CurveView before = v.phi;
    RemoveReport rep = remove_basis(v.phi, 4.1, true);  // nearest: mu = 4.2
    REQUIRE(rep.removed_mu == 4.2);
    REQUIRE(rep.removed_coeff == 12.0);
    REQUIRE(rep.compensated);
    REQUIRE(rep.adjusted.size() == 2);
    REQUIRE(rep.resid_comp <= rep.resid_plain);
    REQUIRE(rep.vresid_comp < rep.vresid_plain);
    REQUIRE(v.phi.n() == 2);

    // the compensated curve is closer to the original than plain removal
    CurveView plain = before;
    plain.mu.erase(plain.mu.begin() + 2);
    plain.coeff.erase(plain.coeff.begin() + 2);
    const double lo = 0.15 * before.rcut, hi = 4.2;
    double se_comp = 0.0, se_plain = 0.0;
    for (int i = 0; i < 512; ++i) {
      const double x = lo + (hi - lo) * i / 511.0;
      se_comp += std::pow(v.phi.value(x) - before.value(x), 2);
      se_plain += std::pow(plain.value(x) - before.value(x), 2);
    }
    REQUIRE(se_comp < se_plain);
  }
  SECTION("zero-coefficient removal is exact and skips compensation") {
    PotEAMView v = parse_pot_view(make_knot_pot());
    add_basis(v.phi, AddSpec{3.6, {}, {}}, true);  // coeff 0
    const CurveView before = v.phi;
    RemoveReport rep = remove_basis(v.phi, 3.6, true);
    REQUIRE(rep.removed_coeff == 0.0);
    REQUIRE_FALSE(rep.compensated);
    REQUIRE(curve_maxdiff(before, v.phi, 0.2, 8.0) == 0.0);
  }
  SECTION("single compensator when only one BF remains") {
    CurveView c;
    c.family = BasisKind::Knot5;
    c.mu = {3.0, 4.5};
    c.coeff = {2.0, -1.0};
    c.rcut = 4.5;
    RemoveReport rep = remove_basis(c, 4.5, true);
    REQUIRE(rep.compensated);
    REQUIRE(rep.adjusted.size() == 1);
    REQUIRE(c.n() == 1);
  }
  SECTION("outermost-knot removal warns and shrinks the derived cutoff") {
    PotEAMView v = parse_pot_view(make_knot_pot());
    RemoveReport rep = remove_basis(v.phi, 5.0, true);
    bool mentions_outermost = false;
    for (const auto &w : rep.warnings)
      if (w.find("OUTERMOST") != std::string::npos) mentions_outermost = true;
    REQUIRE(mentions_outermost);
    Context out = rebuild_pot_context(v);
    REQUIRE(out.get<double>("RCUT2B") == Approx(4.2));
  }
  SECTION("embedding (right-active) removal compensates too") {
    // The fixture has no EMBF0 line -> ZERO: the -2.9 knot has g(0) = 2.9^5,
    // so the region is the widened F domain (see the [embf0] cases).
    PotEAMView v = parse_pot_view(make_knot_pot());
    RemoveReport rep = remove_basis(v.F, -2.9, false);
    REQUIRE(rep.removed_mu == -2.9);
    REQUIRE(rep.compensated);
    REQUIRE(rep.resid_comp <= rep.resid_plain);
    REQUIRE(rep.region_widened);
  }
  SECTION("embedding (right-active) removal, EMBF0 FREE (raw basis)") {
    PotEAMView v = parse_pot_view(with_embf0(make_knot_pot(), "FREE"));
    REQUIRE_FALSE(v.F.zero_ref);
    RemoveReport rep = remove_basis(v.F, -2.9, false);
    REQUIRE(rep.removed_mu == -2.9);
    REQUIRE(rep.compensated);
    REQUIRE(rep.resid_comp <= rep.resid_plain);
    REQUIRE_FALSE(rep.region_widened);
    REQUIRE(rep.notes.empty());
  }
  SECTION("degree array entry removed and all-quintic collapses to absent") {
    PotEAMView v = parse_pot_view(make_knot_pot());
    add_basis(v.phi, AddSpec{3.6, {}, 3}, true);
    REQUIRE(v.phi.deg.size() == 4);
    remove_basis(v.phi, 3.6, true);
    REQUIRE(v.phi.deg.empty());  // remaining degrees are all 5
    Context out = rebuild_pot_context(v);
    REQUIRE(!out.contains("KDEG2B"));
  }
  SECTION("removing the last basis function throws") {
    CurveView c;
    c.family = BasisKind::Knot5;
    c.mu = {3.0};
    c.coeff = {1.0};
    c.rcut = 3.0;
    REQUIRE_THROWS_WITH(remove_basis(c, 3.0, true),
                        Catch::Contains("last basis function"));
  }
}

TEST_CASE_METHOD(EditFixture, "potedit: derived/normalized knot RCUT",
                 "[potedit]") {
  SECTION("knot channels derive, blip channels return nullopt") {
    Context k = make_knot_pot();
    REQUIRE(derived_knot_rcut(k, "2B").value() == Approx(5.0));
    REQUIRE(derived_knot_rcut(k, "MB").value() == Approx(4.8));
    Context b = make_blip_pot();
    REQUIRE(!derived_knot_rcut(b, "2B").has_value());
    REQUIRE(!derived_knot_rcut(b, "MB").has_value());
  }
  SECTION("normalize rewrites a stale stated RCUT and is idempotent") {
    Context k = make_knot_pot();
    k.remove("RCUT2B");
    k.add("RCUT2B", 5.3);  // stale: table cutoff instead of max knot
    REQUIRE(normalize_knot_rcut(k, "test"));
    REQUIRE(k.get<double>("RCUT2B") == Approx(5.0));
    REQUIRE(k.get<double>("RCUTMB") == Approx(4.8));  // was consistent
    REQUIRE(k.get<double>("RCUT2BMAX") == Approx(5.0));  // internal key
    REQUIRE_FALSE(normalize_knot_rcut(k, "test"));  // second call: no-op
  }
  SECTION("blip pots are never modified") {
    Context b = make_blip_pot();
    REQUIRE_FALSE(normalize_knot_rcut(b, "test"));
    REQUIRE(b.get<double>("RCUT2B") == 5.0);
  }
  SECTION("a consistent stated RCUT is left alone") {
    Context k = make_knot_pot();
    REQUIRE_FALSE(normalize_knot_rcut(k, "test"));
    REQUIRE(k.get<double>("RCUT2B") == 5.0);
  }
}

TEST_CASE_METHOD(EditFixture,
                 "potedit: normalize_knot_rcut rewrites ANY difference "
                 "exactly; WARNING only above 1e-9 relative",
                 "[potedit][lossless]") {
  // Every consumer (train, predict, refit --edit, hpo) must evaluate the
  // same cutoff as hpo's per-evaluation derivation, so no dead band: a
  // stated value 3e-13 off (the <= 1.5.0 12-decimal CGRID gap) is fixed
  // too, quietly (INFO); a real mismatch still warns.
  auto &log = tadah::core::Logger::instance();
  const auto old_level = log.getLogLevel();
  log.setLogLevel(tadah::core::LogLevel::WARNING);
  auto normalize_capturing = [](Context &k, std::string &out) {
    std::ostringstream captured;
    std::streambuf *old_buf = std::cout.rdbuf(captured.rdbuf());
    const bool changed = normalize_knot_rcut(k, "test");
    std::cout.rdbuf(old_buf);
    out = captured.str();
    return changed;
  };
  std::string text;
  for (double off : {3e-13, -3e-13, 1e-11, -1e-11}) {
    INFO("offset " << off);
    Context k = make_knot_pot();
    k.remove("RCUT2B");
    k.add("RCUT2B", 5.0 + off);
    REQUIRE(normalize_capturing(k, text));
    CHECK(k.get<double>("RCUT2B") == 5.0);    // exact, not Approx
    CHECK(k.get<double>("RCUT2BMAX") == 5.0);
    CHECK(text.find("WARNING") == std::string::npos);
    REQUIRE_FALSE(normalize_capturing(k, text));  // idempotent
  }
  for (double off : {1e-8, -1e-8, 0.3}) {
    INFO("offset " << off);
    Context k = make_knot_pot();
    k.remove("RCUT2B");
    k.add("RCUT2B", 5.0 + off);
    REQUIRE(normalize_capturing(k, text));
    CHECK(k.get<double>("RCUT2B") == 5.0);
    CHECK(text.find("RCUT2B") != std::string::npos);
    // the stated value is printed exactly (it may differ in digit 13)
    CHECK(text.find(tadah::core::format_roundtrip(5.0 + off)) !=
          std::string::npos);
  }
  log.setLogLevel(old_level);
}

TEST_CASE_METHOD(EditFixture,
                 "potedit: normalize_knot_rcut under a smooth cutoff warns "
                 "that pair_style tadah reads the stated RCUT",
                 "[potedit][lossless]") {
  // Under a smooth RCTYPE the cutoff enters f_c(r) of every pair, and the
  // LAMMPS plugin (never modified) reads the stated RCUT: a legacy knot pot
  // with a 3e-13 gap then predicts differently in LAMMPS than in
  // train/predict/hpo. That must be visible at the default level; under
  // Cut_Dummy (prediction-exact) it stays INFO without the note.
  auto &log = tadah::core::Logger::instance();
  const auto old_level = log.getLogLevel();
  log.setLogLevel(tadah::core::LogLevel::WARNING);
  auto normalize_capturing = [](Context &k, std::string &out) {
    std::ostringstream captured;
    std::streambuf *old_buf = std::cout.rdbuf(captured.rdbuf());
    const bool changed = normalize_knot_rcut(k, "test");
    std::cout.rdbuf(old_buf);
    out = captured.str();
    return changed;
  };
  std::string text;
  for (const char *rctype : {"Cut_Cos", "Cut_Poly"}) {
    for (double off : {3e-13, -3e-13, 0.3}) {
      INFO(rctype << " offset " << off);
      Context k = make_knot_pot();
      k.remove("RCTYPE2B");
      k.add("RCTYPE2B", rctype);
      k.remove("RCUT2B");
      k.add("RCUT2B", 5.0 + off);
      REQUIRE(normalize_capturing(k, text));
      CHECK(k.get<double>("RCUT2B") == 5.0);
      CHECK(text.find("WARNING") != std::string::npos);
      CHECK(text.find("pair_style tadah") != std::string::npos);
      CHECK(text.find(tadah::core::format_roundtrip(5.0 + off)) !=
            std::string::npos);
      REQUIRE_FALSE(normalize_capturing(k, text));  // idempotent, silent
      CHECK(text.empty());
    }
  }
  // the many-body knot channel follows its own RCTYPEMB
  {
    Context k = make_knot_pot();
    k.remove("RCTYPEMB");
    k.add("RCTYPEMB", "Cut_Cos");
    k.remove("RCUTMB");
    k.add("RCUTMB", 4.8 + 3e-13);
    REQUIRE(normalize_capturing(k, text));
    CHECK(text.find("WARNING") != std::string::npos);
    CHECK(text.find("RCUTMB") != std::string::npos);
    CHECK(text.find("pair_style tadah") != std::string::npos);
  }
  // hard cutoff: sub-1e-9 gap stays below WARNING, and no LAMMPS note
  {
    Context k = make_knot_pot();  // RCTYPE2B Cut_Dummy
    k.remove("RCUT2B");
    k.add("RCUT2B", 5.0 + 3e-13);
    log.setLogLevel(tadah::core::LogLevel::INFO);
    REQUIRE(normalize_capturing(k, text));
    CHECK(text.find("INFO") != std::string::npos);
    CHECK(text.find("WARNING") == std::string::npos);
    CHECK(text.find("pair_style tadah") == std::string::npos);
  }
  log.setLogLevel(old_level);
}

TEST_CASE_METHOD(EditFixture, "potedit: pure two-body potentials",
                 "[potedit]") {
  SECTION("parse classifies and splits the weights") {
    PotEAMView v = parse_pot_view(make_pure2b_pot());
    REQUIRE(v.mb == MBKind::None);
    REQUIRE(v.phi.n() == 3);
    REQUIRE(v.phi.coeff == std::vector<double>{1.5, -900.0, 12.0});
    REQUIRE(v.tail.empty());
  }
  SECTION("round-trip without edit is text-identical") {
    Context pot = make_pure2b_pot();
    Context out = rebuild_pot_context(parse_pot_view(pot));
    for (const char *key : {"CGRID2B", "WEIGHTS", "FIXWEIGHT", "TYPE2B"}) {
      REQUIRE(out.size(key) == pot.size(key));
      for (std::size_t i = 0; i < pot.size(key); ++i)
        REQUIRE(out.get<std::string>(key, i) == pot.get<std::string>(key, i));
    }
    // no many-body TERM appears out of nowhere. (Context::contains
    // materialises blueprint defaults, so probing INITMB during the parse
    // leaves an explicit `INITMB false` behind — semantically identical to
    // its absence, and it keeps the pot pair-only for every consumer.)
    for (const char *key : {"TYPEMB", "CGRIDMB", "AMPGRIDMB",
                            "CEMBFUNC", "RCUTMB"}) {
      INFO("key=" << key);
      REQUIRE(!out.contains(key));
    }
    REQUIRE_FALSE(out.get<bool>("INITMB"));
  }
  SECTION("ADD beyond the largest knot: bit-identical + cutoff grows") {
    PotEAMView v = parse_pot_view(make_pure2b_pot());
    const CurveView before = v.phi;
    add_basis(v.phi, AddSpec{6.8, {}, {}}, true);
    REQUIRE(curve_maxdiff(before, v.phi, 0.2, 9.0) == 0.0);
    Context out = rebuild_pot_context(v);
    REQUIRE(out.get<double>("RCUT2B") == Approx(6.8));
    REQUIRE(out.size("WEIGHTS") == 4);
    REQUIRE(out.get<double>("WEIGHTS", 3) == 0.0);
    REQUIRE(out.get<std::string>("FIXINDEX") == "1-4");
    REQUIRE(out.get<int>("TYPE2B", 1) == 4);
  }
  SECTION("REMOVE compensates as usual") {
    PotEAMView v = parse_pot_view(make_pure2b_pot());
    RemoveReport rep = remove_basis(v.phi, 4.1, true);
    REQUIRE(rep.removed_mu == 4.2);
    REQUIRE(rep.compensated);
    Context out = rebuild_pot_context(v);
    REQUIRE(out.size("WEIGHTS") == 2);
    REQUIRE(out.get<std::string>("FIXINDEX") == "1-2");
  }
  SECTION("derived/normalized RCUT works without a many-body term") {
    Context k = make_pure2b_pot();
    k.remove("RCUT2B");
    k.add("RCUT2B", 6.3);  // stale headroom-style value
    REQUIRE(derived_knot_rcut(k, "2B").value() == Approx(5.0));
    REQUIRE(!derived_knot_rcut(k, "MB").has_value());
    REQUIRE(normalize_knot_rcut(k, "test"));
    REQUIRE(k.get<double>("RCUT2B") == Approx(5.0));
  }
}

TEST_CASE_METHOD(EditFixture,
                 "potedit: foreign many-body descriptors (pair-only editing)",
                 "[potedit]") {
  PotEAMView v = parse_pot_view(make_foreign_mb_pot());
  REQUIRE(v.mb == MBKind::Foreign);
  REQUIRE(v.mb_name == "DM_Blip");
  REQUIRE(v.phi.coeff == std::vector<double>{1.5, -900.0, 12.0});
  REQUIRE(v.tail == std::vector<double>{7.0, -3.0});

  const CurveView before = v.phi;
  add_basis(v.phi, AddSpec{3.6, {}, {}}, true);
  REQUIRE(curve_maxdiff(before, v.phi, 0.2, 8.0) == 0.0);
  Context out = rebuild_pot_context(v);
  // pair block grew; the opaque many-body tail rides along unchanged
  REQUIRE(out.size("WEIGHTS") == 6);
  REQUIRE(out.get<double>("WEIGHTS", 3) == 0.0);
  REQUIRE(out.get<double>("WEIGHTS", 4) == 7.0);
  REQUIRE(out.get<double>("WEIGHTS", 5) == -3.0);
  REQUIRE(out.get<std::string>("FIXINDEX") == "1-6");
  // every many-body key survives verbatim
  const Context ref = make_foreign_mb_pot();
  for (const char *key : {"TYPEMB", "CGRIDMB", "SGRIDMB", "RCUTMB",
                          "RCTYPEMB", "INITMB"}) {
    REQUIRE(out.size(key) == ref.size(key));
    for (std::size_t i = 0; i < ref.size(key); ++i)
      REQUIRE(out.get<std::string>(key, i) == ref.get<std::string>(key, i));
  }
}

TEST_CASE_METHOD(EditFixture, "potedit: BIAS/NORM potentials are rejected",
                 "[potedit]") {
  Context b = make_pure2b_pot();
  b.remove("BIAS");
  b.add("BIAS", "true");
  REQUIRE_THROWS_WITH(parse_pot_view(b), Catch::Contains("BIAS"));
  Context n = make_pure2b_pot();
  n.remove("NORM");
  n.add("NORM", "true");
  REQUIRE_THROWS_WITH(parse_pot_view(n), Catch::Contains("NORM"));
}

TEST_CASE_METHOD(EditFixture,
                 "potedit: rebuilt contexts construct real descriptors",
                 "[potedit]") {
  PotEAMView v = parse_pot_view(make_knot_pot());
  add_basis(v.phi, AddSpec{3.6, {}, 3}, true);
  add_basis(v.rho, AddSpec{4.9, {}, {}}, true);
  remove_basis(v.F, 0.9, false);
  Context out = rebuild_pot_context(v);  // internal self-check already ran
  out.add("VERBOSE", 0);
  REQUIRE_NOTHROW(D2_Knot5{out});
  REQUIRE_NOTHROW(DM_REAM{out});
  // TYPEMB count tokens follow the edits
  REQUIRE(out.get<int>("TYPEMB", 2) == 3);  // nC (rho grew)
  REQUIRE(out.get<int>("TYPEMB", 4) == 3);  // nAmp
  REQUIRE(out.get<int>("TYPEMB", 5) == 2);  // nEc (F shrank)
  REQUIRE(out.get<double>("RCUTMB") == Approx(4.9));
  REQUIRE(stored_embf0_text(out) == "ZERO");
}

TEST_CASE_METHOD(EditFixture,
                 "potedit: EMBF0 — DM_REAM potentials always declare it, "
                 "others never do", "[potedit][embf0]") {
  SECTION("legacy DM_REAM potential (no key): ZERO, written explicitly") {
    for (Context pot : {make_knot_pot(), make_blip_pot()}) {
      REQUIRE_FALSE(pot.has_stored("EMBF0"));
      PotEAMView v = parse_pot_view(pot);
      REQUIRE(v.embf0 == EmbF0Mode::Zero);
      REQUIRE(v.F.zero_ref);
      REQUIRE_FALSE(v.phi.zero_ref);
      REQUIRE_FALSE(v.rho.zero_ref);
      const Context out = rebuild_pot_context(v);
      REQUIRE(stored_embf0_text(out) == "ZERO");
      REQUIRE(out.size("EMBF0") == 1);
      REQUIRE(parse_pot_view(out).embf0 == EmbF0Mode::Zero);
    }
  }
  SECTION("EMBF0 FREE is preserved (canonical spelling)") {
    for (Context pot : {with_embf0(make_knot_pot(), "FREE"),
                        with_embf0(make_blip_pot(), "free")}) {
      PotEAMView v = parse_pot_view(pot);
      REQUIRE(v.embf0 == EmbF0Mode::Free);
      REQUIRE_FALSE(v.F.zero_ref);
      const Context out = rebuild_pot_context(v);
      REQUIRE(stored_embf0_text(out) == "FREE");
      REQUIRE(out.size("EMBF0") == 1);
    }
  }
  SECTION("explicit EMBF0 ZERO is preserved") {
    PotEAMView v = parse_pot_view(with_embf0(make_blip_pot(), "ZERO"));
    REQUIRE(v.F.zero_ref);
    REQUIRE(stored_embf0_text(rebuild_pot_context(v)) == "ZERO");
  }
  SECTION("pure two-body and foreign many-body potentials never gain it") {
    for (Context pot : {make_pure2b_pot(), make_foreign_mb_pot()}) {
      const Context out = rebuild_pot_context(parse_pot_view(pot));
      REQUIRE_FALSE(out.has_stored("EMBF0"));
    }
    // a stray key on a potential without DM_REAM is dropped
    for (Context pot : {with_embf0(make_pure2b_pot(), "FREE"),
                        with_embf0(make_foreign_mb_pot(), "ZERO")}) {
      const Context out = rebuild_pot_context(parse_pot_view(pot));
      REQUIRE_FALSE(out.has_stored("EMBF0"));
    }
  }
  SECTION("an invalid EMBF0 is rejected") {
    REQUIRE_THROWS_WITH(parse_pot_view(with_embf0(make_knot_pot(), "ON")),
                        Catch::Contains("EMBF0 must be ZERO or FREE"));
  }
}

TEST_CASE_METHOD(EditFixture,
                 "potedit: the F view evaluates what DM_REAM evaluates "
                 "(EMBF0 ZERO and FREE)", "[potedit][embf0]") {
  const std::vector<double> xs = {0.0, 1e-9, 0.05, 0.1, 0.3, 0.9,
                                  1.4, 2.2,  3.1,  4.5, 6.0};
  // knot fixtures carry c = -2.9 (g(0) = 2.9^5; 2.9^3 when cubic); the blip
  // fixture has two bumps covering 0.
  for (const char *family : {"knot", "knot3", "blip"}) {
    const std::string f(family);
    const Context base = f == "knot"    ? make_knot_pot()
                         : f == "knot3" ? make_knot3_pot()
                                        : make_blip_pot();
    const PotEAMView vz = parse_pot_view(base);
    const PotEAMView vf = parse_pot_view(with_embf0(base, "FREE"));
    REQUIRE(vz.F.zero_ref);
    REQUIRE_FALSE(vf.F.zero_ref);
    const std::vector<double> Ez =
        ream_energies(rebuild_pot_context(vz), vz.F.coeff, xs);
    const std::vector<double> Ef =
        ream_energies(rebuild_pot_context(vf), vf.F.coeff, xs);

    // exactly 0 at rho == 0 under ZERO (view and runtime); raw F(0) != 0
    REQUIRE(vz.F.value(0.0) == 0.0);
    REQUIRE(Ez[0] == 0.0);
    const double F0 = vf.F.value(0.0);
    REQUIRE(F0 != 0.0);

    for (std::size_t i = 0; i < xs.size(); ++i) {
      const double x = xs[i];
      INFO(family << " x=" << x);
      // independent sums of the raw / zero-referenced embedding basis
      double raw = 0.0, zref = 0.0;
      const CurveView &F = vz.F;
      for (std::size_t k = 0; k < F.n(); ++k) {
        const double g = embedding_value(F.family, x, F.width(k), F.mu[k],
                                         F.degree(k));
        const double g0 = embedding_value(F.family, 0.0, F.width(k), F.mu[k],
                                          F.degree(k));
        raw += F.coeff[k] * g;
        zref += F.coeff[k] * (g - g0);
      }
      const double tz = 1e-12 * std::max(1.0, std::fabs(zref));
      const double tf = 1e-12 * std::max(1.0, std::fabs(raw));
      CHECK(std::fabs(vz.F.value(x) - zref) <= tz);
      CHECK(std::fabs(vz.F.value(x) - Ez[i]) <= tz);
      CHECK(std::fabs(vf.F.value(x) - raw) <= tf);
      CHECK(std::fabs(vf.F.value(x) - Ef[i]) <= tf);
      // the two references differ by the constant raw F(0) ...
      CHECK(std::fabs((vf.F.value(x) - vz.F.value(x)) - F0) <=
            1e-12 * std::max({1.0, std::fabs(raw), std::fabs(F0)}));
      // ... and share every derivative exactly
      CHECK(vz.F.deriv(x) == vf.F.deriv(x));
      CHECK(vz.F.deriv2(x) == vf.F.deriv2(x));
    }
  }
}

TEST_CASE_METHOD(EditFixture, "potedit: REMOVE under EMBF0 ZERO",
                 "[potedit][embf0]") {
  SECTION("a bump covering 0 widens the region; the compensated curve "
          "stays close to the original zero-referenced curve") {
    CurveView zero;
    zero.family = BasisKind::Blip;
    zero.zero_ref = true;
    zero.mu = {0.0, 0.1, 0.2, 1.0, 2.0};
    zero.eta = {1.0, 1.0, 1.0, 1.0, 1.0};
    zero.coeff = {3.0, -2.0, 1.5, 0.7, -0.4};
    const CurveView before = zero;

    RemoveReport rep = remove_basis(zero, 0.1, false);
    REQUIRE(rep.removed_mu == 0.1);
    REQUIRE(rep.region_widened);
    REQUIRE(rep.notes.size() == 1);
    // largest upper support: 2.0 + 2/1.0
    CHECK(rep.notes[0].find("widened to the F domain [0, 4]") !=
          std::string::npos);
    REQUIRE(rep.compensated);
    std::vector<std::size_t> adj = rep.adjusted;
    std::sort(adj.begin(), adj.end());
    REQUIRE(adj == std::vector<std::size_t>{0, 2});  // the two neighbours
    REQUIRE(rep.resid_comp <= rep.resid_plain);
    REQUIRE(rep.vresid_comp < rep.vresid_plain);
    REQUIRE(zero.value(0.0) == 0.0);

    // plain removal leaves the constant c_r g_r(0) beyond the bump; the
    // compensated curve stays close over the whole F domain (tail included)
    CurveView plain = before;
    plain.mu.erase(plain.mu.begin() + 1);
    plain.eta.erase(plain.eta.begin() + 1);
    plain.coeff.erase(plain.coeff.begin() + 1);
    double dev_comp = 0.0, dev_plain = 0.0;
    for (int i = 0; i < 2001; ++i) {
      const double x = 6.0 * i / 2000.0;
      dev_comp = std::max(dev_comp, std::fabs(zero.value(x) - before.value(x)));
      dev_plain =
          std::max(dev_plain, std::fabs(plain.value(x) - before.value(x)));
    }
    CHECK(dev_plain > 1.0);
    CHECK(dev_comp < 0.05 * dev_plain);

    // FREE twin: the bump's own support (HEAD), no note
    CurveView free = before;
    free.zero_ref = false;
    RemoveReport rf = remove_basis(free, 0.1, false);
    REQUIRE(rf.compensated);
    REQUIRE_FALSE(rf.region_widened);
    REQUIRE(rf.notes.empty());
  }
  SECTION("a COMPENSATOR covering 0 widens the region too; a support-edge "
          "g(0) at rounding level does not") {
    // dx = 0.1, eta = 10 (support +-0.2): removing mu = 0.2 (g(0) = B(-2)
    // = 0) adjusts the neighbours 0.1 (g(0) = B(-1) = 0.25) and 0.3, whose
    // g(0) != 0 shifts F by a constant beyond 0.3 -- the whole F domain.
    CurveView zero;
    zero.family = BasisKind::Blip;
    zero.zero_ref = true;
    for (int k = 0; k <= 20; ++k) {
      zero.mu.push_back(0.1 * k);
      zero.eta.push_back(10.0);
      zero.coeff.push_back(std::sin(0.7 * k) + 0.3);
    }
    const CurveView before = zero;
    RemoveReport rep = remove_basis(zero, 0.2, false);
    REQUIRE(rep.removed_mu == 0.2);
    REQUIRE(rep.compensated);
    REQUIRE(rep.region_widened);
    REQUIRE(rep.notes.size() == 1);
    CHECK(rep.notes[0].find("a compensating function covers rho = 0") !=
          std::string::npos);
    REQUIRE(zero.value(0.0) == 0.0);
    // the compensation objective (value/1st/2nd rows) over the widened
    // region, which includes the constant tail the compensators' g(0) would
    // otherwise leave unaccounted
    REQUIRE(rep.resid_comp <= rep.resid_plain);

    // FREE twin: local region as in 1.5.0
    CurveView free = before;
    free.zero_ref = false;
    RemoveReport rf = remove_basis(free, 0.2, false);
    REQUIRE_FALSE(rf.region_widened);

    // mu = 0.5: removed and compensators (0.4, 0.6) all have g(0) = 0
    CurveView far = before;
    RemoveReport rr = remove_basis(far, 0.5, false);
    REQUIRE_FALSE(rr.region_widened);
    REQUIRE(rr.notes.empty());
  }
  SECTION("fixtures: g(0) != 0 widens under ZERO, never under FREE") {
    struct Case {
      const char *tag;
      bool knot;
      double centre;
      bool widened_zero;
    };
    for (const Case &cs : {Case{"knot -2.9 (g(0) = 2.9^5)", true, -2.9, true},
                           // g(0) = 0, but its compensator -2.9 covers 0
                           Case{"knot 0.9 (compensator -2.9 covers 0)", true,
                                0.9, true},
                           Case{"blip 0.1 (covers 0)", false, 0.1, true},
                           // clear of 0, but its compensators 0.9 and 0.1
                           // cover 0 (the non-widening case is the dense
                           // grid section above)
                           Case{"blip 2.9 (compensators cover 0)", false, 2.9,
                                true}}) {
      INFO(cs.tag);
      const Context base = cs.knot ? make_knot_pot() : make_blip_pot();
      PotEAMView vz = parse_pot_view(base);
      RemoveReport rz = remove_basis(vz.F, cs.centre, false);
      REQUIRE(rz.removed_mu == cs.centre);
      REQUIRE(rz.region_widened == cs.widened_zero);
      REQUIRE(rz.notes.size() == (cs.widened_zero ? 1u : 0u));
      REQUIRE(rz.resid_comp <= rz.resid_plain);
      REQUIRE(vz.F.value(0.0) == 0.0);
      // the rebuilt potential still evaluates F(0) = 0 at runtime
      const std::vector<double> E =
          ream_energies(rebuild_pot_context(vz), vz.F.coeff, {0.0, 0.7});
      REQUIRE(E[0] == 0.0);
      CHECK(std::fabs(E[1] - vz.F.value(0.7)) <=
            1e-12 * std::max(1.0, std::fabs(E[1])));

      PotEAMView vf = parse_pot_view(with_embf0(base, "FREE"));
      RemoveReport rf = remove_basis(vf.F, cs.centre, false);
      REQUIRE(rf.removed_mu == cs.centre);
      REQUIRE_FALSE(rf.region_widened);
      REQUIRE(rf.notes.empty());
      REQUIRE(rf.resid_comp <= rf.resid_plain);
    }
  }
}

TEST_CASE_METHOD(EditFixture,
                 "potedit: ADD keeps DM_REAM energies bitwise (EMBF0 ZERO "
                 "and FREE)", "[potedit][embf0]") {
  const std::vector<double> rhos = {0.0, 1e-9, 0.05, 0.3, 0.9,
                                    1.4, 2.2,  3.1,  4.5};
  struct Case {
    const char *tag;
    bool knot;
    double centre;
  };
  for (const char *mode : {"ZERO", "FREE"}) {
    for (const Case &cs : {Case{"knot, centre > 0", true, 0.4},
                           Case{"knot, centre < 0 (g(0) != 0)", true, -1.5},
                           Case{"blip, bump clear of 0", false, 2.4},
                           Case{"blip, bump covering 0", false, 0.3}}) {
      INFO(mode << ": " << cs.tag);
      PotEAMView v = parse_pot_view(
          with_embf0(cs.knot ? make_knot_pot() : make_blip_pot(), mode));
      REQUIRE(v.F.zero_ref == (std::string(mode) == "ZERO"));
      const CurveView before = v.F;
      const std::vector<double> E0 =
          ream_energies(rebuild_pot_context(v), v.F.coeff, rhos);
      add_basis(v.F, AddSpec{cs.centre, {}, {}}, false);
      REQUIRE(v.F.coeff.back() == 0.0);
      const std::vector<double> E1 =
          ream_energies(rebuild_pot_context(v), v.F.coeff, rhos);
      for (std::size_t i = 0; i < rhos.size(); ++i) {
        INFO("rho=" << rhos[i]);
        REQUIRE(E1[i] == E0[i]);
      }
      REQUIRE(curve_maxdiff(before, v.F, 0.0, 6.0) == 0.0);
    }
  }
}
