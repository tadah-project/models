#include "catch2/catch.hpp"
#include <tadah/core/periodic_table.h>
#include <tadah/models/eam/eam_setfl_io.h>
#include <tadah/models/eam/eam_refit.h>
#include <tadah/models/eam/eam_lattice.h>
#include <tadah/models/descriptors/d_basis_functions.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

using namespace tadah::models;

namespace {
struct RefitFixture { RefitFixture() { tadah::core::PeriodicTable::initialize(); } };

// Minimal single-element setfl (header + tiny tabulated block) with a chosen
// lattice token and lattice constant. Enough for read_setfl(); do NOT call
// build_splines() on it (nr=nrho=2 is below the spline stencil).
void write_min_setfl(const std::string &path, const std::string &lat,
                     double a0) {
  std::ofstream os(path);
  os << "# minimal setfl for tests\n# comment 2\n# comment 3\n"
     << "1 Cu\n"
     << "2 0.5 2 0.5 1.0\n"
     << "29 63.55 " << a0 << " " << lat << "\n"
     << "0.0 0.1\n0.2 0.3\n0.4 0.5\n";
}

// Relative RMS error of f_fit vs f_tab sampled over [a,b], normalised by the
// peak-to-peak amplitude of f_tab over the same range.
template <class FFit, class FTab>
double rel_rmse(FFit fit, FTab tab, double a, double b, int N) {
  double se = 0.0, lo = 1e300, hi = -1e300;
  for (int i = 0; i < N; ++i) {
    const double x = a + (b - a) * i / (N - 1);
    const double t = tab(x);
    const double d = fit(x) - t;
    se += d * d;
    lo = std::min(lo, t); hi = std::max(hi, t);
  }
  const double amp = (hi - lo) > 0 ? (hi - lo) : 1.0;
  return std::sqrt(se / N) / amp;
}
} // namespace

TEST_CASE_METHOD(RefitFixture, "EAMRefit reproduces the Ta2 EAM functions", "[refit]") {
  EAMSetflIO eam;
  eam.read_setfl("tests_data/eam.alloy");
  eam.build_splines();

  EAMRefitConfig cfg;          // defaults
  EAMRefit refit(eam, cfg);
  refit.fit();

  const double rcut = eam.rcut();
  const double r_lo = 0.20 * rcut;
  const double rhomax = eam.rhomax();

  const double phi_rel = rel_rmse([&](double r){ return refit.phi_fit(r); },
                                  [&](double r){ return eam.phi_tab(r); },
                                  r_lo, rcut, 400);
  const double rho_rel = rel_rmse([&](double r){ return refit.rho_fit(r); },
                                  [&](double r){ return eam.rho_tab(r); },
                                  r_lo, rcut, 400);
  // F fidelity contract (default windowed placement): TIGHT inside the
  // physically reachable window [0, rho_reach], crude-but-bounded on the
  // unreachable tail (it keeps only a small uniform basis by design).
  const double F_rel_win = rel_rmse([&](double x){ return refit.F_fit(x); },
                                    [&](double x){ return eam.F_tab(x); },
                                    0.0, refit.rho_reach(), 400);
  const double F_rel_full = rel_rmse([&](double x){ return refit.F_fit(x); },
                                     [&](double x){ return eam.F_tab(x); },
                                     0.0, rhomax, 400);

  std::cout << "[refit] phi rel_rmse=" << phi_rel
            << " rho rel_rmse=" << rho_rel
            << " F rel_rmse win/full=" << F_rel_win << "/" << F_rel_full
            << std::endl;
  std::cout << "[refit] phi maxabs=" << refit.phi_stats().maxabs
            << " rho maxabs=" << refit.rho_stats().maxabs
            << " F maxabs=" << refit.F_stats().maxabs << std::endl;

  REQUIRE(phi_rel < 0.02);
  REQUIRE(rho_rel < 0.015);
  REQUIRE(F_rel_win < 0.002);
  REQUIRE(F_rel_full < 0.10);
}

TEST_CASE_METHOD(RefitFixture,
                 "EAMRefit windowed placement: reachability + F(0) pin",
                 "[refit]") {
  EAMSetflIO eam;
  eam.read_setfl("tests_data/eam.alloy");
  eam.build_splines();

  EAMRefitConfig cfg;  // defaults: placement = "window", constrain_F0 = true
  EAMRefit refit(eam, cfg);
  refit.fit();

  // reachability diagnostics populated and sane
  REQUIRE(refit.lattice_known());
  REQUIRE(refit.r_reach() > 0.0);
  REQUIRE(refit.r_reach() < eam.rcut());
  REQUIRE(refit.rho_reach() > 0.0);
  REQUIRE(refit.rho_reach() <= eam.rhomax());

  // F window fidelity is tight; the unreachable tail is crude BY DESIGN
  const double F_win = rel_rmse([&](double x){ return refit.F_fit(x); },
                                [&](double x){ return eam.F_tab(x); },
                                0.0, refit.rho_reach(), 400);
  REQUIRE(F_win < 2e-3);

  // isolated-atom reference: F(0) pinned to the table's value
  REQUIRE(std::fabs(refit.F_fit(0.0) - eam.F_tab(0.0)) < 1e-4);

  // legacy uniform placement still works
  EAMRefitConfig ucfg;
  ucfg.placement = "uniform";
  EAMRefit urefit(eam, ucfg);
  urefit.fit();
  const double F_uni = rel_rmse([&](double x){ return urefit.F_fit(x); },
                                [&](double x){ return eam.F_tab(x); },
                                0.0, eam.rhomax(), 400);
  REQUIRE(F_uni < 0.015);
}

TEST_CASE_METHOD(RefitFixture,
                 "EAMRefit compact-mode injection round-trips",
                 "[refit]") {
  EAMSetflIO eam;
  eam.read_setfl("tests_data/eam.alloy");
  eam.build_splines();
  EAMRefitConfig cfg;
  cfg.n2b = 5; cfg.nrho = 4; cfg.nF = 6;  // injected sizes
  EAMRefit refit(eam, cfg);
  refit.analyze_reachability();
  REQUIRE(refit.r_reach() > 0.0);

  std::vector<double> mu2{2.0, 2.8, 3.6, 4.4, 5.2}, et2(5, 1.2),
      c2{1.0, -0.5, 0.2, -0.1, 0.05};
  std::vector<double> mur{2.2, 3.1, 4.0, 4.9}, etr(4, 1.0),
      cr{0.4, 0.3, 0.2, 0.1};
  std::vector<double> muF{0.0, 0.6, 1.2, 1.8, 2.4, 3.0}, etF(6, 1.5),
      cF{0.0, -1.0, -1.5, -1.2, -0.6, 0.2};
  refit.set_blip_fit(EAMRefit::Fn::PHI, mu2, et2, c2);
  refit.set_blip_fit(EAMRefit::Fn::RHO, mur, etr, cr);
  refit.set_blip_fit(EAMRefit::Fn::F, muF, etF, cF);
  refit.finalize_injected(2.0, 3.0);
  refit.set_rho_join(3.0);

  // pot context carries the injected grids and WEIGHTS = [c2b | cmb]
  tadah::core::Context c = refit.build_pot_context();
  REQUIRE(c.size("CGRID2B") == 5);
  REQUIRE(c.size("AMPGRIDMB") == 4);
  REQUIRE(c.size("CEMBFUNC") == 6);
  REQUIRE(c.size("WEIGHTS") == 11);
  REQUIRE(c.get<double>("WEIGHTS", 0) == Approx(1.0));
  REQUIRE(c.get<double>("WEIGHTS", 5) == Approx(0.0));

  // F_setfl: identical to F_fit inside the join, original+CONSTANT offset
  // beyond it, C0 at the join
  REQUIRE(refit.F_setfl(2.0) == Approx(refit.F_fit(2.0)));
  REQUIRE(refit.F_setfl(3.0) == Approx(refit.F_fit(3.0)).margin(1e-12));
  const double off1 = refit.F_setfl(4.0) - eam.F_tab(4.0);
  const double off2 = refit.F_setfl(6.0) - eam.F_tab(6.0);
  REQUIRE(off1 == Approx(off2).margin(1e-10));

  // stats were computed on the requested windows
  REQUIRE(refit.phi_stats().n == 1000);
  REQUIRE(refit.F_stats().n == 1000);
}

TEST_CASE_METHOD(RefitFixture,
                 "EAMRefit setfl continuation below r_lo follows the original",
                 "[refit]") {
  EAMSetflIO eam;
  eam.read_setfl("tests_data/eam.alloy");
  eam.build_splines();
  EAMRefitConfig cfg;
  EAMRefit refit(eam, cfg);
  refit.fit();

  const double rl = refit.r_lo();
  // C0 at the join
  REQUIRE(refit.phi_setfl(rl) == Approx(refit.phi_fit(rl)).margin(1e-10));
  REQUIRE(refit.rho_setfl(rl) == Approx(refit.rho_fit(rl)).margin(1e-10));
  // below the join: original shape with a CONSTANT phi offset
  const double off1 = refit.phi_setfl(0.5 * rl) - eam.phi_tab(0.5 * rl);
  const double off2 = refit.phi_setfl(0.8 * rl) - eam.phi_tab(0.8 * rl);
  REQUIRE(off1 == Approx(off2).margin(1e-8));
  // z2r finite at the origin (the writer tabulates z2_setfl from r=0)
  REQUIRE(std::isfinite(refit.z2_setfl(0.0)));
  // identity above the join
  REQUIRE(refit.phi_setfl(1.5 * rl) == Approx(refit.phi_fit(1.5 * rl)));
}

TEST_CASE_METHOD(RefitFixture, "EAMRefit assembles a structurally valid pot context", "[refit]") {
  EAMSetflIO eam;
  eam.read_setfl("tests_data/eam.alloy");
  eam.build_splines();

  EAMRefitConfig cfg;
  EAMRefit refit(eam, cfg);
  refit.fit();

  tadah::core::Context c = refit.build_pot_context();
  // WATOM must be pinned to 1.0: when absent, the HPO pipeline fills it from
  // the blueprint default (atomic number), scaling DM_REAM's density by Z.
  REQUIRE(c.contains("WATOM"));
  REQUIRE(c.get<double>("WATOM", 0) == Approx(1.0));
  REQUIRE(c.get<std::string>("MODEL", 0) == "M_BLR");
  REQUIRE(c.get<std::string>("MODEL", 1) == "BF_Linear");
  REQUIRE(c.get<std::string>("TYPE2B", 0) == "D2_Blip");
  REQUIRE(c.get<std::string>("TYPEMB", 0) == "DM_REAM");
  REQUIRE(c.size("CGRID2B") == static_cast<std::size_t>(cfg.n2b));
  REQUIRE(c.size("AMPGRIDMB") == static_cast<std::size_t>(cfg.nrho));
  REQUIRE(c.size("CEMBFUNC") == static_cast<std::size_t>(cfg.nF));
  // WEIGHTS = [c2b | cmb] = n2b + nF
  REQUIRE(c.size("WEIGHTS") ==
          static_cast<std::size_t>(cfg.n2b + cfg.nF));
}


TEST_CASE_METHOD(RefitFixture,
                 "EAMRefit honours injected knot degrees and emits DGRID keys",
                 "[refit][knotmix]") {
  EAMSetflIO eam;
  eam.read_setfl("tests_data/eam.alloy");
  eam.build_splines();
  EAMRefitConfig cfg;
  cfg.basis = "KnotMix";
  EAMRefit refit(eam, cfg);
  refit.analyze_reachability();

  const std::vector<double> mu = {3.0, 4.0, 5.0};
  const std::vector<double> eta = {1.0, 1.0, 1.0};
  const std::vector<double> c = {0.2, -0.1, 0.05};
  const std::vector<int> deg = {5, 3, 5};
  refit.set_blip_fit(EAMRefit::Fn::PHI, mu, eta, c, deg);
  refit.set_blip_fit(EAMRefit::Fn::RHO, mu, eta, c);       // all 5
  refit.set_blip_fit(EAMRefit::Fn::F, {0.1, 0.5}, {1.0, 1.0},
                     {0.3, -0.2}, {3, 5});
  refit.finalize_injected(2.0, 1.0);

  // (a) degree honoured in evaluation: middle phi atom must be CUBIC
  const double r = 3.5;
  const double expect = c[0] * K5L(r - mu[0]) +
                        c[1] * K3L(r - mu[1]) +
                        c[2] * K5L(r - mu[2]);
  REQUIRE(refit.phi_fit(r) == Approx(expect).margin(1e-12));

  // (b) F right-active with cubic first atom
  const double rho = 0.7;
  const double expectF = 0.3 * K3R(rho - 0.1) +
                         -0.2 * K5R(rho - 0.5);
  REQUIRE(refit.F_fit(rho) == Approx(expectF).margin(1e-12));

  // (c) writer: TYPE2B family name, bcode, DGRID emission (and absence)
  tadah::core::Context pot = refit.build_pot_context();
  REQUIRE(pot.get<std::string>("TYPE2B") == "D2_Knot5");
  REQUIRE(pot.get<int>("TYPEMB", 1) == 2);
  REQUIRE(pot.size("KDEG2B") == 3);
  tadah::core::v_type d2b(pot.size("KDEG2B"));
  pot.get("KDEG2B", d2b);
  REQUIRE(d2b[0] == 5.0);
  REQUIRE(d2b[1] == 3.0);
  REQUIRE(d2b[2] == 5.0);
  REQUIRE(!pot.contains("KDEGMB"));    // all-quintic: no key (back-compat)
  REQUIRE(pot.size("KDEGEMB") == 2);
  // knot family: NO width keys anywhere (eta is a pure gauge)
  REQUIRE(!pot.contains("SGRID2B"));
  REQUIRE(!pot.contains("SGRIDMB"));
  REQUIRE(!pot.contains("SEMBFUNC"));
  // TYPE2B contract: D2_Knot5 <nC> <nKdeg> El El
  REQUIRE(pot.get<int>("TYPE2B", 1) == 3);
  REQUIRE(pot.get<int>("TYPE2B", 2) == 3);
}

TEST_CASE_METHOD(RefitFixture, "EAMRefit emits per-family TYPE2B names",
                 "[refit]") {
  EAMSetflIO eam;
  eam.read_setfl("tests_data/eam.alloy");
  eam.build_splines();
  for (const auto &[basis, t2b, bcode] :
       {std::tuple<std::string, std::string, int>{"Blip", "D2_Blip", 0},
        {"Blip5", "D2_Blip5", 3},
        {"Knot5", "D2_Knot5", 2},
        {"Gaussian", "D2_BP", 1}}) {
    EAMRefitConfig cfg;
    cfg.basis = basis;
    EAMRefit refit(eam, cfg);
    // inject rather than fit(): the dense LSQ path is Blip/Gaussian-only —
    // this test exercises the WRITER's family dispatch
    refit.set_blip_fit(EAMRefit::Fn::PHI, {3.0, 4.0}, {1.0, 1.0},
                       {0.1, 0.2});
    refit.set_blip_fit(EAMRefit::Fn::RHO, {3.0, 4.0}, {1.0, 1.0},
                       {0.1, 0.2});
    refit.set_blip_fit(EAMRefit::Fn::F, {0.2, 0.6}, {1.0, 1.0},
                       {0.1, 0.2});
    refit.finalize_injected(2.0, 1.0);
    tadah::core::Context pot = refit.build_pot_context();
    INFO("basis=" << basis);
    REQUIRE(pot.get<std::string>("TYPE2B") == t2b);
    REQUIRE(pot.get<int>("TYPEMB", 1) == bcode);
    if (basis == "Knot5") {
      REQUIRE(!pot.contains("SGRID2B"));
      REQUIRE(!pot.contains("SGRIDMB"));
    } else {
      REQUIRE(pot.size("SGRID2B") == 2);
      REQUIRE(pot.size("SGRIDMB") == 2);
    }
  }
}

TEST_CASE_METHOD(RefitFixture,
                 "EAMRefit derives knot-family RCUT from the largest knot",
                 "[refit]") {
  EAMSetflIO eam;
  eam.read_setfl("tests_data/eam.alloy");
  eam.build_splines();

  // Knot family: the basis vanishes at its largest knot, so the emitted
  // RCUT2B/RCUTMB must be the largest knot position, NOT the table cutoff.
  {
    EAMRefitConfig cfg;
    cfg.basis = "Knot5";
    EAMRefit refit(eam, cfg);
    refit.set_blip_fit(EAMRefit::Fn::PHI, {3.0, 4.7}, {1.0, 1.0}, {0.1, 0.2});
    refit.set_blip_fit(EAMRefit::Fn::RHO, {2.8, 4.1}, {1.0, 1.0}, {0.1, 0.2});
    refit.set_blip_fit(EAMRefit::Fn::F, {0.2, 0.6}, {1.0, 1.0}, {0.1, 0.2});
    refit.finalize_injected(2.0, 1.0);
    tadah::core::Context pot = refit.build_pot_context();
    REQUIRE(pot.get<double>("RCUT2B") == Approx(4.7));
    REQUIRE(pot.get<double>("RCUTMB") == Approx(4.1));
  }
  // Bump families keep the table cutoff: the hard truncation at rcut is
  // part of the fitted model.
  {
    EAMRefitConfig cfg;
    cfg.basis = "Blip";
    EAMRefit refit(eam, cfg);
    refit.set_blip_fit(EAMRefit::Fn::PHI, {3.0, 4.7}, {1.0, 1.0}, {0.1, 0.2});
    refit.set_blip_fit(EAMRefit::Fn::RHO, {2.8, 4.1}, {1.0, 1.0}, {0.1, 0.2});
    refit.set_blip_fit(EAMRefit::Fn::F, {0.2, 0.6}, {1.0, 1.0}, {0.1, 0.2});
    refit.finalize_injected(2.0, 1.0);
    tadah::core::Context pot = refit.build_pot_context();
    REQUIRE(pot.get<double>("RCUT2B") == Approx(eam.rcut()));
    REQUIRE(pot.get<double>("RCUTMB") == Approx(eam.rcut()));
  }
}

TEST_CASE("EAMSetflIO normalises the lattice label to lower case",
          "[refit][setfl]") {
  const std::string p = "tmp_min_upper.eam.alloy";
  write_min_setfl(p, "FCC", 3.615);
  EAMSetflIO eam;
  eam.read_setfl(p);
  std::remove(p.c_str());
  REQUIRE(eam.lattice() == "fcc");        // canonical, used by the crystal code
  REQUIRE(eam.lattice_raw() == "FCC");    // original case preserved
  REQUIRE(eam_lattice::supported(eam.lattice()));
  REQUIRE(eam.lattice_param() == Approx(3.615));
  REQUIRE_NOTHROW(
      eam_lattice::require_crystal(eam.lattice_raw(), eam.lattice_param()));
}

TEST_CASE("EAMSetflIO: an unspecified lattice (NONE / a0=0) fails the crystal "
          "gate", "[refit][setfl]") {
  const std::string p = "tmp_min_none.eam.fs";
  write_min_setfl(p, "NONE", 0.0);
  EAMSetflIO eam;
  eam.read_setfl(p);                       // parsing the file itself is fine
  std::remove(p.c_str());
  REQUIRE(eam.lattice() == "none");
  REQUIRE(eam.lattice_raw() == "NONE");
  REQUIRE_FALSE(eam_lattice::supported(eam.lattice()));
  REQUIRE(eam.lattice_param() == Approx(0.0));
  // the gate the refit engine applies must reject it
  REQUIRE_THROWS(
      eam_lattice::require_crystal(eam.lattice_raw(), eam.lattice_param()));
}

TEST_CASE_METHOD(RefitFixture,
                 "crystal_epa: fitted EAM reproduces the cohesive energy",
                 "[refit]") {
  // The self-test's cohesive point: crystal_epa with the original tabulated
  // functions vs the fitted native functions, on the header lattice at a0.
  EAMSetflIO eam;
  eam.read_setfl("tests_data/eam.alloy");  // Ta, bcc, a0=3.304
  eam.build_splines();
  EAMRefitConfig cfg;                      // dense Blip defaults
  EAMRefit refit(eam, cfg);
  refit.fit();

  const std::string lat = eam.lattice();
  const double a0 = eam.lattice_param(), rc = eam.rcut();
  REQUIRE(eam_lattice::supported(lat));

  const double e_orig = eam_lattice::crystal_epa(
      lat, a0, rc, [&](double r) { return eam.phi_tab(r); },
      [&](double r) { return eam.rho_tab(r); },
      [&](double x) { return eam.F_tab(x); });
  const double e_fit = eam_lattice::crystal_epa(
      lat, a0, rc, [&](double r) { return refit.phi_fit(r); },
      [&](double r) { return refit.rho_fit(r); },
      [&](double x) { return refit.F_fit(x); });

  REQUIRE(std::isfinite(e_orig));
  REQUIRE(std::isfinite(e_fit));
  const double delta_mev = 1000.0 * std::fabs(e_fit - e_orig);
  std::cout << "[refit] cohesive original=" << e_orig << " fit=" << e_fit
            << " delta=" << delta_mev << " meV/atom" << std::endl;
  REQUIRE(delta_mev < 1.0);    // measured ~0.01 meV/atom (dense Blip on Ta);
                               // a broken self-test/fit lands orders larger
}
