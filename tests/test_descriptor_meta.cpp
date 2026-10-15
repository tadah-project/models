// Catch2 v2 tests for the per-descriptor metadata schema of the HPO
// bound-validation pipeline. Verifies:
//   * D_Base::meta() default returns nullopt
//   * D2_Blip::meta() returns a DescriptorMeta with the expected
//     parameter set, suggester output, and column counter
//   * DM_Blip::meta() likewise + cross-param size check
//   * DM_mEAD<F_mEnv<Cut_SinInv>>::meta() (i.e. DM_mBlipEnv) flags a
//     collapsed ρ: a radial GTO centre (CGRIDMB) outside the data's
//     neighbour distances makes every atom's density vanish

#include "catch2/catch.hpp"

#include <tadah/models/descriptors/d_base.h>
#include <tadah/models/descriptors/d2/d2_blip.h>
#include <tadah/models/descriptors/dm/dm_blip.h>
#include <tadah/models/descriptors/dm/dm_mEAD.h>
#include <tadah/models/descriptors/dm/dm_func.h>
#include <tadah/models/descriptors/dm/dm_all.h>
#include <tadah/models/cutoffs.h>

#include <tadah/core/context.h>
#include <tadah/core/stats/data_stats.h>

#include <algorithm>
#include <cmath>
#include <memory>

using namespace tadah::models;
using tadah::core::Context;
using tadah::core::stats::DataStats;
using tadah::core::stats::PerPairStats;

// ----------------------------------------------------------------------
// Helper: build a synthetic DataStats with a single Ta-Ta pair whose
// r_ij percentiles bracket typical bcc Ta neighbour distances.
// ----------------------------------------------------------------------
static DataStats make_ta_like_stats() {
  DataStats d;
  d.n_atoms = 1000;
  d.n_structures = 50;
  d.elements = {"Ta"};
  PerPairStats pp;
  pp.elements = {"Ta", "Ta"};
  pp.n_pairs  = 50000;
  pp.r_min    = 1.61;
  pp.r_max    = 6.5;
  // p1 .. p99 (mean, std)
  pp.r_percentiles = {2.33, 2.55, 2.85, 3.0, 5.03, 5.5, 6.30, 6.4, 6.48,
                      4.5, 1.5};
  // One KDE peak near 2.86 with FWHM ~ 0.3
  tadah::core::stats::KdePeak peak;
  peak.position = 2.86;
  peak.density  = 1.0;
  peak.fwhm     = 0.30;
  pp.r_peaks    = {peak};
  d.per_pair    = {pp};
  return d;
}

static DataStats make_ta_like_stats_with_collapsed_rho() {
  auto d = make_ta_like_stats();
  // ρ-collapse: every percentile near zero
  d.rho_percentiles_per_descriptor["DM_mBlipEnv"] = {
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-10, 0.0, 0.0};
  return d;
}

static DataStats make_ta_like_stats_with_healthy_rho() {
  auto d = make_ta_like_stats();
  d.rho_percentiles_per_descriptor["DM_mBlipEnv"] = {
      // p1, p5, p10, p25, p50, p75, p90, p95, p99, mean, std
      1.0, 2.0, 4.6, 5.5, 5.81, 6.2, 6.68, 6.9, 7.2,
      5.5, 1.0};
  return d;
}

static Context make_context_with_optim_keys() {
  Context c;
  c.add<bool>("INIT2B", true);
  c.add<bool>("INITMB", true);
  c.add<double>("RCUT2B", 5.6);
  c.add<double>("RCUTMB", 5.6);
  c.add<double>("RCUTENV", 5.0);
  c.add<std::string>("RCTYPE2B", "Cut_Poly");
  c.add<std::string>("RCTYPEMB", "Cut_Poly");
  // Vector grids — 7 entries each.
  for (int i = 0; i < 7; ++i) {
    c.add<double>("CGRID2B", 0.5 * (i + 1));
    c.add<double>("SGRID2B", 1.0);
    c.add<double>("CGRIDMB", 0.5 * (i + 1));
    c.add<double>("SGRIDMB", 1.0);
    c.add<double>("CEMBFUNC", 1.1 + 0.1 * i);
    c.add<double>("SEMBFUNC", 1.0);
  }
  c.add<std::string>("TYPE2B", "D2_Blip"); c.add<int>("TYPE2B", 7);
  c.add<int>("TYPE2B", 7); c.add<std::string>("TYPE2B", "Ta");
  c.add<std::string>("TYPE2B", "Ta");
  c.add<std::string>("TYPEMB", "DM_Blip"); c.add<int>("TYPEMB", 0);
  return c;
}

// ----------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------

TEST_CASE("D_Base::meta default is nullopt for descriptors that do not override") {
  // DM_Dummy is registered but its metadata has not been ported.
  // factory<>() requires a definition of registry — DM_Dummy lives under
  // Registry<D_Base> which is already defined in d2_all.cpp. Construct
  // and assert default behaviour.
  auto* p = tadah::core::factory<D_Base>("DM_Dummy");
  REQUIRE(p != nullptr);
  REQUIRE_FALSE(p->meta().has_value());
  delete p;
}

TEST_CASE("D2_Blip::meta returns a populated DescriptorMeta") {
  D2_Blip d;
  auto m = d.meta();
  REQUIRE(m.has_value());
  REQUIRE(m->label == "D2_Blip");
  REQUIRE(m->params.size() == 3);
  // RCUT2B, CGRID2B, SGRID2B in this order
  REQUIRE(m->params[0].key == "RCUT2B");
  REQUIRE(m->params[1].key == "CGRID2B");
  REQUIRE(m->params[2].key == "SGRID2B");
  // Hard limits: SGRID2B must be > 0
  REQUIRE(m->params[2].hard.lo.has_value());
  REQUIRE(m->params[2].hard.lo.value() > 0.0);
  // RCUT2B is critical
  REQUIRE(m->params[0].crit == Criticality::Critical);
  // CGRID2B is per-pair
  REQUIRE(m->params[1].scope == PairScope::PerElementPair);
  // expected_columns lambda is set
  REQUIRE(static_cast<bool>(m->expected_columns));
}

TEST_CASE("D2_Blip suggester anchors RCUT2B to the data's r_ij percentiles") {
  D2_Blip d;
  auto m = d.meta();
  REQUIRE(m.has_value());
  Context ctx = make_context_with_optim_keys();
  SeedView sv(ctx);
  auto data = make_ta_like_stats();
  auto bracket = m->params[0].suggest(data, sv, -1);
  // Lower bound should sit between 2NN peak (none in our synthetic
  // KdePeaks; only peak at 2.86) and p99.
  REQUIRE(bracket.first  > 0.0);
  REQUIRE(bracket.second  > bracket.first);
  REQUIRE(bracket.second <= data.per_pair[0].r_percentiles[8]);  // <= p99
}

TEST_CASE("D2_Blip cross-param check fails on size mismatch") {
  D2_Blip d;
  auto m = d.meta();
  REQUIRE(m.has_value());
  REQUIRE(m->cross_param_checks.size() == 1);

  Context ctx = make_context_with_optim_keys();
  // Make CGRID2B 8 entries while SGRID2B is 7 -> mismatch
  ctx.add<double>("CGRID2B", 99.0);
  SeedView sv(ctx);
  auto v = m->cross_param_checks[0](DataStats{}, sv);
  REQUIRE(is_fail(v));
  REQUIRE(verdict_message(v).find("CGRID2B") != std::string::npos);
}

TEST_CASE("DM_Blip::meta returns expected_columns scaling with TYPEMB Lmax") {
  DM_Blip d;
  auto m = d.meta();
  REQUIRE(m.has_value());
  REQUIRE(m->label == "DM_Blip");
  Context ctx = make_context_with_optim_keys();
  SeedView sv(ctx);
  // 7 CGRIDMB entries x (Lmax+1) = 7 * 1 = 7 (Lmax=0 in the ctx above)
  REQUIRE(m->expected_columns(sv) == 7);
}

TEST_CASE("DM_mBlipEnv::meta: CGRIDMB is critical and the cross-param check "
          "fails a collapsed ρ") {
  // The descriptor type is DM_mEAD<F_mEnv<Cut_SinInv>>. We exercise it
  // via meta() through a DM_mEAD<F_mEnv<Cut_SinInv>> instance directly.
  DM_mEAD<F_mEnv<Cut_SinInv>> d;
  auto m = d.meta();
  REQUIRE(m.has_value());
  REQUIRE(m->label == "DM_mBlipEnv");
  REQUIRE(m->params.size() == 6);     // RCUTMB, CGRIDMB, SGRIDMB, RCUTENV, CEMBFUNC, SEMBFUNC

  // CGRIDMB is Critical: the pre-flight audit warns when it is left out
  // of the OPTIM search.
  auto cgridmb = std::find_if(m->params.begin(), m->params.end(),
      [](const ParamMeta& p){ return p.key == "CGRIDMB"; });
  REQUIRE(cgridmb != m->params.end());
  REQUIRE(cgridmb->crit == Criticality::Critical);

  // Cross-param check fires Fail when ρ percentiles all collapsed.
  REQUIRE(!m->cross_param_checks.empty());
  Context ctx = make_context_with_optim_keys();
  SeedView sv(ctx);
  auto v_collapse = m->cross_param_checks[0](make_ta_like_stats_with_collapsed_rho(), sv);
  REQUIRE(is_fail(v_collapse));
  REQUIRE(verdict_message(v_collapse).find("collapse") != std::string::npos);

  // And returns Ok when ρ is healthy.
  auto v_ok = m->cross_param_checks[0](make_ta_like_stats_with_healthy_rho(), sv);
  REQUIRE(is_ok(v_ok));
}

TEST_CASE("DM_mBlipEnv per-param check rejects out-of-range CGRIDMB") {
  DM_mEAD<F_mEnv<Cut_SinInv>> d;
  auto m = d.meta();
  REQUIRE(m.has_value());

  auto cgridmb = std::find_if(m->params.begin(), m->params.end(),
      [](const ParamMeta& p){ return p.key == "CGRIDMB"; });
  REQUIRE(cgridmb != m->params.end());
  REQUIRE(static_cast<bool>(cgridmb->check));

  // Empty context with CGRIDMB = 0.0, a centre outside the data's r_ij
  // support (r_min=1.61)
  Context ctx;
  ctx.add<double>("CGRIDMB", 0.0);
  SeedView sv(ctx);
  auto v = cgridmb->check(make_ta_like_stats(), sv);
  REQUIRE(is_fail(v));
  REQUIRE(verdict_message(v).find("CGRIDMB") != std::string::npos);
  REQUIRE(verdict_message(v).find("ρ") != std::string::npos);
}

TEST_CASE("DM_mBlipEnv suggester places CGRIDMB on the first KDE peak") {
  DM_mEAD<F_mEnv<Cut_SinInv>> d;
  auto m = d.meta();
  REQUIRE(m.has_value());

  auto cgridmb = std::find_if(m->params.begin(), m->params.end(),
      [](const ParamMeta& p){ return p.key == "CGRIDMB"; });
  REQUIRE(cgridmb != m->params.end());
  REQUIRE(static_cast<bool>(cgridmb->suggest));

  Context ctx = make_context_with_optim_keys();
  SeedView sv(ctx);
  auto bracket = cgridmb->suggest(make_ta_like_stats(), sv, -1);
  // Peak at 2.86, fwhm 0.30 -> bracket [2.56, 3.16]
  REQUIRE(bracket.first  == Approx(2.56).margin(0.01));
  REQUIRE(bracket.second == Approx(3.16).margin(0.01));
}

TEST_CASE("DM_mBlipEnv expected_columns equals CEMBFUNC size") {
  DM_mEAD<F_mEnv<Cut_SinInv>> d;
  auto m = d.meta();
  REQUIRE(m.has_value());
  Context ctx = make_context_with_optim_keys();    // 7 CEMBFUNC entries
  SeedView sv(ctx);
  REQUIRE(m->expected_columns(sv) == 7);
}

// ----------------------------------------------------------------
// Coverage audit
//
// Assert that every registered D_Base subclass either has a populated
// meta() or is on the exempt list (Dummy / mJoin meta-descriptors).
// This is the regression that catches future descriptor additions
// being landed without a `meta()` declaration.
// ----------------------------------------------------------------

#include <set>
#include <tadah/core/registry.h>

TEST_CASE("descriptor coverage: every registered descriptor has meta()") {
  // Descriptors that intentionally return `nullopt` (placeholders /
  // meta-descriptors that delegate to nested children).
  const std::set<std::string> exempt = {
    "D2_Dummy", "D2_mJoin",
    "DM_Dummy", "DM_mJoin",
  };

  auto names = tadah::core::Registry<D_Base>::names();
  REQUIRE(!names.empty());                     // sanity

  std::vector<std::string> missing;
  for (const auto& name : names) {
    if (exempt.count(name)) continue;
    std::unique_ptr<D_Base> d(tadah::core::factory<D_Base>(name));
    if (!d->meta().has_value()) missing.push_back(name);
  }

  CAPTURE(missing);
  REQUIRE(missing.empty());
}

// ----------------------------------------------------------------
// Width bounds follow the units of the width coordinate. The compact
// bump families (D2_Blip, D2_Blip5, DM_Blip) evaluate b(eta (r - c)) with
// support |eta (r - c)| < 2: eta is an inverse distance (Å^-1), so the
// suggested bounds scale as 1/length when every distance in the data is
// scaled. The Gaussian families (D2_BP, DM_EAD) evaluate exp(-eta (r -
// c)^2): eta is an inverse squared distance (Å^-2), scaling as 1/length^2.
// ----------------------------------------------------------------

namespace {
DataStats scaled_stats(const DataStats &d, double s) {
  DataStats out = d;
  for (auto &pp : out.per_pair) {
    pp.r_min *= s;
    pp.r_max *= s;
    for (auto &v : pp.r_percentiles) v *= s;
    for (auto &peak : pp.r_peaks) {
      peak.position *= s;
      peak.fwhm *= s;
    }
  }
  return out;
}
}  // namespace

/* Pins: for each family the suggested SGRID2B/SGRIDMB bounds of data whose
 * distances are doubled are the bounds of the original data divided by 2
 * (bump families, Å^-1) or by 4 (Gaussian families, Å^-2), matching the
 * units the metadata declares. Before the fix the bump families declared
 * Å^-1 but returned 0.25/s^2 .. 25/s^2 (D2) and 0.25/w^2 .. 9/w^2 (DM), so
 * doubling the distances divided their bounds by 4. */
TEST_CASE("Suggested width bounds scale with the units of the width coordinate",
          "[unified-dm][metadata]") {
  const DataStats base = make_ta_like_stats();
  const DataStats twice = scaled_stats(base, 2.0);
  Context ctx = make_context_with_optim_keys();
  SeedView sv(ctx);
  struct Family { const char *label; const char *key; const char *units; int power; };
  for (const Family &f : {Family{"D2_Blip", "SGRID2B", "Å^-1", 1},
                          Family{"D2_Blip5", "SGRID2B", "Å^-1", 1},
                          Family{"DM_Blip", "SGRIDMB", "Å^-1", 1},
                          Family{"D2_BP", "SGRID2B", "Å^-2", 2},
                          Family{"DM_EAD", "SGRIDMB", "Å^-2", 2}}) {
    CAPTURE(f.label);
    std::unique_ptr<D_Base> d(tadah::core::factory<D_Base>(f.label));
    auto meta = d->meta();
    REQUIRE(meta.has_value());
    const ParamMeta *width = nullptr;
    for (const auto &p : meta->params)
      if (p.key == f.key) width = &p;
    REQUIRE(width != nullptr);
    REQUIRE(width->units == f.units);
    const auto b1 = width->suggest(base, sv, 0);
    const auto b2 = width->suggest(twice, sv, 0);
    REQUIRE(b1.first > 0.0);
    REQUIRE(b1.second > b1.first);
    const double factor = std::pow(2.0, f.power);
    REQUIRE(b2.first == Approx(b1.first / factor).epsilon(1e-12));
    REQUIRE(b2.second == Approx(b1.second / factor).epsilon(1e-12));
  }
}
