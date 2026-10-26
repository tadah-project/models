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
#include <sstream>

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

/* Pins: a narrow radial GTO at the origin fails the CGRIDMB check: with
 * CGRIDMB 0 and SGRIDMB 6, exp(-6 r_min^2) is 1.8e-7 at the smallest
 * neighbour distance (r_min 1.61), so ρ is ~0 for every atom; the message
 * names CGRIDMB and ρ.
 * Before the fix: passes as well (a pin: the check compared the centre
 * alone with [r_min, r_max], so a centre at 0 failed whatever its width;
 * this case had no SGRIDMB). */
TEST_CASE("DM_mBlipEnv per-param check rejects out-of-range CGRIDMB") {
  DM_mEAD<F_mEnv<Cut_SinInv>> d;
  auto m = d.meta();
  REQUIRE(m.has_value());

  auto cgridmb = std::find_if(m->params.begin(), m->params.end(),
      [](const ParamMeta& p){ return p.key == "CGRIDMB"; });
  REQUIRE(cgridmb != m->params.end());
  REQUIRE(static_cast<bool>(cgridmb->check));

  // CGRIDMB = 0.0 with a narrow width: below 1e-6 at r_min = 1.61
  Context ctx;
  ctx.add<double>("CGRIDMB", 0.0);
  ctx.add<double>("SGRIDMB", 6.0);
  SeedView sv(ctx);
  auto v = cgridmb->check(make_ta_like_stats(), sv);
  REQUIRE(is_fail(v));
  REQUIRE(verdict_message(v).find("CGRIDMB") != std::string::npos);
  REQUIRE(verdict_message(v).find("ρ") != std::string::npos);
}

/* Pins: without a seed centre (no CGRIDMB) the CGRIDMB suggestion of
 * DM_mBlipEnv is the first KDE peak +- its FWHM, 2.56 to 3.16.
 * Before the fix: passes as well (a pin; the case read the seed's CGRIDMB
 * 0.5, which got the first-peak range too: the seed's centre was not
 * read, and a centre of 0.5 now gets a range about it, as the case below
 * pins). */
TEST_CASE("DM_mBlipEnv suggester places CGRIDMB on the first KDE peak") {
  DM_mEAD<F_mEnv<Cut_SinInv>> d;
  auto m = d.meta();
  REQUIRE(m.has_value());

  auto cgridmb = std::find_if(m->params.begin(), m->params.end(),
      [](const ParamMeta& p){ return p.key == "CGRIDMB"; });
  REQUIRE(cgridmb != m->params.end());
  REQUIRE(static_cast<bool>(cgridmb->suggest));

  Context ctx = make_context_with_optim_keys();
  ctx.remove("CGRIDMB");
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
 * doubling the distances divided their bounds by 4. The DM_Poly families
 * follow their radial basis: DM_PS_SQ (blip) Å^-1, DM_GPS_SQ (Gaussian)
 * Å^-2; before that fix DM_PS_SQ declared the Gaussian's Å^-2 and returned
 * its 0.25/w^2 .. 9/w^2 for a blip width. */
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
                          Family{"DM_EAD", "SGRIDMB", "Å^-2", 2},
                          Family{"DM_PS_SQ", "SGRIDMB", "Å^-1", 1},
                          Family{"DM_GPS_SQ", "SGRIDMB", "Å^-2", 2}}) {
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

// ----------------------------------------------------------------
// The reach of a radial density channel. A density sums its radial
// function over the neighbours within its cutoff, so a channel is judged
// by what it gives at the neighbour distances of the data, not by where
// its centre lies: a Gaussian exp(-SGRIDMB (r - CGRIDMB)^2) misses them
// when it is below 1e-6 (the floor of the ρ-collapse check) at every
// distance from r_min to the smaller of r_max and RCUTMB, and a blip
// b(SGRIDMB (r - CGRIDMB)) when its support, |r - CGRIDMB| < 2/SGRIDMB,
// holds none of them.
// ----------------------------------------------------------------

namespace {

/* The CGRIDMB check of the descriptor @p label on make_ta_like_stats()
 * (r_min 1.61, r_max 6.5) for one channel: centre @p c, width @p eta,
 * RCUTMB 5.6. */
HealthVerdict reach_verdict(const std::string &label, double c, double eta) {
  std::unique_ptr<D_Base> d(tadah::core::factory<D_Base>(label));
  auto meta = d->meta();
  REQUIRE(meta.has_value());
  const ParamMeta *centre = nullptr;
  for (const auto &p : meta->params)
    if (p.key == "CGRIDMB") centre = &p;
  REQUIRE(centre != nullptr);
  REQUIRE(static_cast<bool>(centre->check));
  Context ctx;
  ctx.add<double>("CGRIDMB", c);
  ctx.add<double>("SGRIDMB", eta);
  ctx.add<double>("RCUTMB", 5.6);
  SeedView sv(ctx);
  return centre->check(make_ta_like_stats(), sv);
}

}  // namespace

/* Pins: for the Gaussian densities (DM_mBlipEnv, DM_EAD, DM_mRLR, DM_GPS_SQ,
 * DM_GOW_SQ) the CGRIDMB check passes a centre at the origin whose width
 * reaches the neighbours (CGRIDMB 0, SGRIDMB 0.2: exp(-0.2 * 1.61^2) = 0.6,
 * the Ta seed of the HPO tests) and a centre beyond the cutoff whose tail
 * reaches them (CGRIDMB 8, SGRIDMB 0.2: 0.32 at RCUTMB 5.6), and fails a
 * narrow Gaussian at the origin (SGRIDMB 6: 1.8e-7 at r_min) and a centre
 * far outside (CGRIDMB 30), each message naming CGRIDMB and SGRIDMB.
 * Before the fix: CGRIDMB 0 and 8 failed whatever the width ("CGRIDMB=0.000000
 * is outside the dataset's r_ij range [...] -- ρ will be ~0 for every atom";
 * the centre alone was compared with [r_min, r_max]), and no message named
 * SGRIDMB. */
TEST_CASE("Density reach: a Gaussian channel fails only when it is below "
          "1e-6 at every neighbour distance",
          "[metadata][density_reach]") {
  for (const char *label :
       {"DM_mBlipEnv", "DM_EAD", "DM_mRLR", "DM_GPS_SQ", "DM_GOW_SQ"}) {
    CAPTURE(label);
    CHECK(is_ok(reach_verdict(label, 0.0, 0.2)));
    CHECK(is_ok(reach_verdict(label, 8.0, 0.2)));
    for (const auto &[c, eta] : {std::pair{0.0, 6.0}, std::pair{30.0, 0.2}}) {
      CAPTURE(c, eta);
      const HealthVerdict v = reach_verdict(label, c, eta);
      CHECK(is_fail(v));
      CHECK(verdict_message(v).find("CGRIDMB") != std::string::npos);
      CHECK(verdict_message(v).find("SGRIDMB") != std::string::npos);
    }
  }
}

/* Pins: for the blip densities (DM_Blip, DM_PS_SQ, DM_OW_SQ) the CGRIDMB
 * check passes a centre below r_min whose support reaches it (CGRIDMB 1,
 * SGRIDMB 1: |r - 1| < 2 holds 1.61 to 3) and one inside the data (CGRIDMB
 * 3), and fails a blip at the origin whose support ends before r_min
 * (SGRIDMB 2: |r| < 1) and one centred far outside (CGRIDMB 30), each
 * message naming CGRIDMB and SGRIDMB.
 * Before the fix: CGRIDMB 1 failed (below r_min 1.61, its support
 * ignored), and no message named SGRIDMB. */
TEST_CASE("Density reach: a blip channel fails only when its support holds "
          "no neighbour distance",
          "[metadata][density_reach]") {
  for (const char *label : {"DM_Blip", "DM_PS_SQ", "DM_OW_SQ"}) {
    CAPTURE(label);
    CHECK(is_ok(reach_verdict(label, 1.0, 1.0)));
    CHECK(is_ok(reach_verdict(label, 3.0, 1.0)));
    for (const auto &[c, eta] : {std::pair{0.0, 2.0}, std::pair{30.0, 1.0}}) {
      CAPTURE(c, eta);
      const HealthVerdict v = reach_verdict(label, c, eta);
      CHECK(is_fail(v));
      CHECK(verdict_message(v).find("CGRIDMB") != std::string::npos);
      CHECK(verdict_message(v).find("SGRIDMB") != std::string::npos);
    }
  }
}

// ----------------------------------------------------------------
// The suggested ranges of a Gaussian channel whose seed centre is not at
// the origin: the centre is searched about the seed's own centre, with the
// width of the neighbour shell nearest to it.
// ----------------------------------------------------------------

namespace {

/* Neighbour distances with three shells, at 2.86 A (FWHM 0.3), 3.30 A
 * (FWHM 0.24) and 4.67 A (FWHM 0.36); r_min 2.25 A, r_max 6.5 A. */
DataStats three_shell_stats() {
  DataStats d = make_ta_like_stats();
  auto &pp = d.per_pair[0];
  pp.r_min = 2.25;
  tadah::core::stats::KdePeak a, b, c;
  a.position = 2.86;
  a.density = 1.0;
  a.fwhm = 0.30;
  b.position = 3.30;
  b.density = 0.7;
  b.fwhm = 0.24;
  c.position = 4.67;
  c.density = 0.5;
  c.fwhm = 0.36;
  pp.r_peaks = {a, b, c};
  return d;
}

const ParamMeta *param_of(const DescriptorMeta &m, const char *key) {
  for (const auto &p : m.params)
    if (p.key == key) return &p;
  return nullptr;
}

}  // namespace

/* Pins: for a Gaussian channel of DM_EAD whose seed centre is not 0, the
 * CGRIDMB suggestion is the seed's centre +- the FWHM of the neighbour
 * shell nearest to it, its low end not below r_min and its high end not
 * above RCUTMB unless the seed's centre is (then the range ends at the
 * centre): on three_shell_stats() with CGRIDMB 1 2.9 3.3 4.7 5.9 (SGRIDMB
 * 4, RCUTMB 5.6) the ranges are 1 to 1.3, 2.6 to 3.2, 3.06 to 3.54, 4.34
 * to 5.06 and 5.54 to 5.9, each holding its seed's centre, and each
 * rationale names the seed's centre and the nearest shell; the SGRIDMB
 * suggestion of an entry is 0.25/w^2 to 9/w^2 with w the FWHM of that
 * shell (4.7: w = 0.36, 1.929 to 69.44), its rationale naming the shell.
 * Before the fix: every centre got the first shell's range, 2.56 to 3.16
 * (the first peak +- its FWHM), which holds none of the seed's centres but
 * 2.9, and every width that of the first shell (2.78 to 100), each
 * rationale naming the first neighbour shell. */
TEST_CASE("Gaussian density suggestions are centred on the seed's centre",
          "[metadata][density_reach][suggest]") {
  std::unique_ptr<D_Base> d(tadah::core::factory<D_Base>("DM_EAD"));
  const auto meta = d->meta();
  REQUIRE(meta.has_value());
  const ParamMeta *centre = param_of(*meta, "CGRIDMB");
  const ParamMeta *width = param_of(*meta, "SGRIDMB");
  REQUIRE(centre != nullptr);
  REQUIRE(width != nullptr);
  Context ctx;
  ctx.add("CGRIDMB", "1 2.9 3.3 4.7 5.9");
  ctx.add("SGRIDMB", "4 4 4 4 4");
  ctx.add<double>("RCUTMB", 5.6);
  const SeedView sv(ctx);
  const DataStats data = three_shell_stats();
  struct Want { double c, lo, hi, shell; };
  const Want wants[] = {{1.0, 1.0, 1.3, 2.86},
                        {2.9, 2.6, 3.2, 2.86},
                        {3.3, 3.06, 3.54, 3.30},
                        {4.7, 4.34, 5.06, 4.67},
                        {5.9, 5.54, 5.9, 4.67}};
  for (int i = 0; i < 5; ++i) {
    const Want &w = wants[i];
    CAPTURE(i, w.c);
    const BoundSuggestion s = centre->suggest(data, sv, i);
    INFO(s.rationale);
    CHECK_FALSE(s.none);
    CHECK(s.first == Approx(w.lo).epsilon(1e-12));
    CHECK(s.second == Approx(w.hi).epsilon(1e-12));
    CHECK(s.first <= w.c);
    CHECK(s.second >= w.c);
    std::ostringstream shell;
    shell << w.shell;
    CHECK(s.rationale.find("seed's centre") != std::string::npos);
    CHECK(s.rationale.find("nearest") != std::string::npos);
    CHECK(s.rationale.find(shell.str()) != std::string::npos);
  }
  const BoundSuggestion sw = width->suggest(data, sv, 3);
  INFO(sw.rationale);
  CHECK(sw.first == Approx(0.25 / (0.36 * 0.36)).epsilon(1e-12));
  CHECK(sw.second == Approx(9.0 / (0.36 * 0.36)).epsilon(1e-12));
  CHECK(sw.rationale.find("nearest") != std::string::npos);
  CHECK(sw.rationale.find("4.67") != std::string::npos);
}
