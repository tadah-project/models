/* FEAT-02 goldens: the DM_REAM density paths against the bits of the density
 * before RHOFORM (every tadah before 1.6.0).
 *
 * ream_goldens_data.h holds, for 5 families (Blip, Blip5, Knot5, Knot3,
 * Gaussian) x 3 cutoffs (Cut_Dummy, Cut_Cos, Cut_Poly) x 2 amplitude sets
 * (signed, all non-negative), calc_rho's increment, the calc_dXijdri force
 * prefactor and eam_density's rho and rho' at 56-64 distances in (0, rcut],
 * captured as hex floats from the tree before the density forms (MODELS
 * 4c4144b) with glibc on x86-64. ream_density_before_rhoform.cpp is a
 * frozen copy of that tree's code: it gives the density before RHOFORM on
 * the platform the test runs on. RHOFORM SIGNED must reproduce it bit for
 * bit on both amplitude sets, and RHOFORM NONNEG on the non-negative set
 * (every cutoff here is >= 0, so the max(f_c, 0) clamp is inert).
 *
 * The density calls the C library's cos (Cut_Cos), sin (its derivative)
 * and exp (Gaussian), whose results the C standard does not fix to the bit:
 * musl's round some arguments one ulp away from glibc's, so that calc_rho
 * of the Blip Cut_Cos signed case at r = 0x1.e60aa64c2f838p-1 is
 * 0x1.4f8f2a2a6d567p-2 with musl against the stored 0x1.4f8f2a2a6d566p-2,
 * for the code before RHOFORM and for SIGNED alike. Hence the density
 * before RHOFORM is computed here, on the platform, and the stored goldens
 * anchor the frozen copy: bit for bit with glibc, with which they were
 * captured, and with any C library within the bound kLibmBound states.
 */
#include "tadah_test_skip.hpp"  // includes catch2/catch.hpp
#include "ream_density_before_rhoform.h"
#include "ream_goldens_data.h"
#include "ream_goldens_libm.h"

#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/core/number_format.h>
#include <tadah/core/periodic_table.h>
#include <tadah/core/registry.h>
#include <tadah/models/cutoffs.h>
#include <tadah/models/descriptors/dm/dm_ream.h>

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace tadah::core;
using namespace tadah::models;

namespace {

namespace before = ream_before_rhoform;

struct GoldenFixture {
  GoldenFixture() { PeriodicTable::initialize(); }
};

std::uint64_t bits(double v) {
  std::uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return u;
}

std::string hex(double v) {
  std::ostringstream os;
  os << std::hexfloat << v;
  return os.str();
}

int count_tokens(const char *s) {
  std::istringstream is(s);
  std::string t;
  int n = 0;
  while (is >> t) ++n;
  return n;
}

std::vector<double> numbers(const char *s) {
  std::istringstream is(s);
  std::vector<double> out;
  std::string t;
  while (is >> t) {
    double v = 0.0;
    REQUIRE(parse_number(t, v));  // the reader of DM_REAM's grids
    out.push_back(v);
  }
  return out;
}

// The context goldens/ream_goldens.cpp built (make_ctx), plus RHOFORM.
Context golden_ctx(const ream_goldens::Case &gc, const char *form) {
  const bool knot = gc.bcode == 2;
  const int n = count_tokens(gc.cgrid);
  Context c;
  c.add("VERBOSE", 0);
  c.add("ATOM", "Ta");
  c.add("WATOM", 0.8125);
  c.add("INITMB", "true");
  c.add("RCUTMB", 5.5);
  c.add("RCTYPEMB", gc.cutoff);
  c.add("TYPEMB", "DM_REAM");
  c.add("TYPEMB", gc.bcode);
  c.add("TYPEMB", n); c.add("TYPEMB", knot ? 0 : n); c.add("TYPEMB", n);
  c.add("TYPEMB", 3); c.add("TYPEMB", knot ? 0 : 3);
  c.add("TYPEMB", knot ? n : 0); c.add("TYPEMB", 0);
  c.add("TYPEMB", "Ta");
  c.add("TYPEMB", "Ta");
  c.add("CGRIDMB", gc.cgrid);
  if (!knot) c.add("SGRIDMB", gc.sgrid);
  c.add("AMPGRIDMB", gc.amp);
  c.add("CEMBFUNC", "0.0 0.5 1.0");
  if (!knot) c.add("SEMBFUNC", "1.0 1.0 1.0");
  if (knot) c.add("KDEGMB", gc.kdeg);
  c.add("RHOFORM", form);
  return c;
}

// The same density for the frozen copy of the code before RHOFORM.
before::Density golden_density(const ream_goldens::Case &gc) {
  before::Density d;
  d.bcode = gc.bcode;
  d.cutoff = gc.cutoff;
  d.rcut = 5.5;
  d.w = 0.8125;
  d.c = numbers(gc.cgrid);
  d.s = numbers(gc.sgrid);
  d.a = numbers(gc.amp);
  d.kdeg = numbers(gc.kdeg);
  return d;
}

// The density before RHOFORM at every golden point of a case, on this
// platform.
std::vector<before::Values> before_rhoform(const ream_goldens::Case &gc) {
  const before::Density d = golden_density(gc);
  const before::Libm lm = before::platform_libm();
  std::vector<before::Values> out;
  for (std::size_t k = 0; k < gc.npts; ++k)
    out.push_back(before::evaluate(d, gc.pts[k].r, lm));
  return out;
}

// Every golden point of one case through a DM_REAM in @p form against the
// density before RHOFORM on this platform; returns the number of values
// whose bits differ and describes the first one.
std::size_t mismatches(const ream_goldens::Case &gc, const char *form,
                       std::string &first) {
  Context c = golden_ctx(gc, form);
  DM_REAM dm(c);
  const std::vector<before::Values> ref = before_rhoform(gc);
  const int Z = PeriodicTable::find_by_symbol("Ta").Z;
  const Vec3d unit(1.0, 0.0, 0.0);
  std::size_t bad = 0;
  for (std::size_t k = 0; k < gc.npts; ++k) {
    const double r = gc.pts[k].r;

    rho_type rho;
    dm.init_rhoi(rho);
    rho.set_zero();
    dm.calc_rho(Z, Z, r, r * r, unit, rho, 1.0);

    rho_type rhoi;
    dm.init_rhoi(rhoi);
    rhoi.set_zero();
    for (std::size_t m = 0; m < dm.size(); ++m) rhoi(dm.get_drfidx() + m) = 1.0;
    fd_type fd(dm.size());
    fd.set_zero();
    dm.calc_dXijdri(Z, Z, r, r * r, unit, rhoi, fd, 1.0);

    double er, erp;
    dm.eam_density(Z, Z, r, er, erp);

    const double got[4] = {rho(dm.get_rfidx()), fd(0, 0), er, erp};
    const double want[4] = {ref[k].rho, ref[k].pref, ref[k].eam_rho,
                            ref[k].eam_rhop};
    const char *what[4] = {"calc_rho", "calc_dXijdri pref", "eam_density rho",
                           "eam_density rho'"};
    for (int q = 0; q < 4; ++q) {
      if (bits(got[q]) == bits(want[q])) continue;
      if (bad++ == 0)
        first = std::string(what[q]) + " at r = " + hex(r) + ": got " +
                hex(got[q]) + ", before RHOFORM " + hex(want[q]);
    }
  }
  return bad;
}

/* The anchor's bound where the C library is not the capture's: every value
 * of the frozen copy is within kLibmBound * DBL_EPSILON * scale of its
 * golden, scale from before::libm_scale() (the sum over the value's terms
 * of the magnitudes of the C library's factors times the rest of the
 * term). Two C libraries whose cos, sin and exp are within one ulp of the
 * exact result differ by at most one ulp, eps relative, per call.
 * A value adds 6 terms (each case has 6 grid entries), each a product of
 * the bare basis value or derivative with f_c or f_c', and for the
 * Gaussian with exp: at most two such factors per term, which move it by 2
 * eps of its scale each (f_c of Cut_Cos, 0.5 (cos + 1), by 2 eps of
 * (|cos| + 1)/2), and at most 4 roundings after them, 1 eps each; then 5
 * additions and the weight (and 1/r for the prefactor), 1 eps each: at
 * most 4 + 4 + 5 + 2 = 15 eps of the scale. The bound is about twice that.
 * A change of a formula, not only of its rounding, moves values far
 * further. */
constexpr double kLibmBound = 32.0;

struct AnchorCount {
  std::size_t values = 0;
  std::size_t outside_bound = 0;
  std::size_t not_bitwise = 0;
  std::string first_outside;
  std::string first_not_bitwise;
};

// The frozen copy through @p lm against every stored golden of a case.
void anchor(const ream_goldens::Case &gc, const before::Libm &lm,
            AnchorCount &n) {
  const before::Density d = golden_density(gc);
  for (std::size_t k = 0; k < gc.npts; ++k) {
    const ream_goldens::Point &p = gc.pts[k];
    const before::Values v = before::evaluate(d, p.r, lm);
    const before::Values sc = before::libm_scale(d, p.r);
    const double got[4] = {v.rho, v.pref, v.eam_rho, v.eam_rhop};
    const double want[4] = {p.rho, p.pref, p.eam_rho, p.eam_rhop};
    const double scale[4] = {sc.rho, sc.pref, sc.eam_rho, sc.eam_rhop};
    const char *what[4] = {"calc_rho", "calc_dXijdri pref", "eam_density rho",
                           "eam_density rho'"};
    for (int q = 0; q < 4; ++q) {
      ++n.values;
      const std::string at = std::string(gc.family) + " " + gc.cutoff +
                             (gc.nonneg_amps ? " non-negative" : " signed") +
                             ", " + what[q] + " at r = " + hex(p.r) +
                             ": computed " + hex(got[q]) + ", stored " +
                             hex(want[q]);
      if (bits(got[q]) != bits(want[q]) && n.not_bitwise++ == 0)
        n.first_not_bitwise = at;
      const double bound = kLibmBound * DBL_EPSILON * scale[q];
      if (!(std::fabs(got[q] - want[q]) <= bound) && n.outside_bound++ == 0)
        n.first_outside = at + ", bound " + hex(bound);
    }
  }
}

// A C library whose cos, sin and exp round results one ulp away from this
// platform's: every result up (M = 1), every result down (M = 2), or by
// the argument's two lowest bits (M = 0: 01 up, 10 down, else unmoved).
template <int M>
double one_ulp_away(double x, double y) {
  std::uint64_t u = bits(x) & 3u;
  if (M == 1) u = 1;
  if (M == 2) u = 2;
  if (u == 1) return std::nextafter(y, std::numeric_limits<double>::infinity());
  if (u == 2) return std::nextafter(y, -std::numeric_limits<double>::infinity());
  return y;
}
template <int M> double moved_cos(double x) { return one_ulp_away<M>(x, std::cos(x)); }
template <int M> double moved_sin(double x) { return one_ulp_away<M>(x, std::sin(x)); }
template <int M> double moved_exp(double x) { return one_ulp_away<M>(x, std::exp(x)); }
template <int M> before::Libm moved_libm() {
  return before::Libm{&moved_cos<M>, &moved_sin<M>, &moved_exp<M>};
}

// The results of one function of ream_goldens_libm.h against @p f, the C
// library's: the number of arguments at which they differ, and the first.
std::size_t differing(const ream_goldens_libm::Result *tab, std::size_t n,
                      double (*f)(double), const char *name,
                      std::string &first) {
  std::size_t d = 0;
  for (std::size_t i = 0; i < n; ++i) {
    volatile double x = tab[i].x;  // the call happens here, not at compile time
    const double y = f(x);
    if (bits(y) == bits(tab[i].y)) continue;
    ++d;
    if (first.empty())
      first = std::string(name) + "(" + hex(tab[i].x) + ") = " + hex(y) +
              ", recorded " + hex(tab[i].y);
  }
  return d;
}

}  // namespace

/* Pins: the C library this program runs with gives one or more of the
 * results of cos, sin and exp recorded in ream_goldens_libm.h differently
 * (the functions the frozen copy calls, before::platform_libm()). Run only
 * by the CTest entry tadah.models_ream_goldens_glibc_no_fma_libm
 * (tests/CMakeLists.txt), next to tadah.models_ream_goldens_glibc_no_fma,
 * which runs the [goldens] cases under the same setting: it shows that
 * those cases ran with other results than the recorded ones, and is
 * Skipped where the C library gives every recorded result (the setting
 * then changed nothing).
 * Before the fix: new; passes with glibc 2.43 on x86-64 under
 * GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA (1 result of 501 differs:
 * exp(-0x1.1839f3fdcf5f8p+1) = 0x1.cabf14e03b9b2p-4, recorded
 * 0x1.cabf14e03b9b3p-4), and is Skipped without it. */
TEST_CASE("DM_REAM goldens: this C library's cos, sin and exp differ from "
          "the recorded results",
          "[.ream_goldens_libm]") {
  namespace lib = ream_goldens_libm;
  const before::Libm lm = before::platform_libm();
  const std::size_t ncos = sizeof(lib::kCos) / sizeof(lib::kCos[0]);
  const std::size_t nsin = sizeof(lib::kSin) / sizeof(lib::kSin[0]);
  const std::size_t nexp = sizeof(lib::kExp) / sizeof(lib::kExp[0]);
  std::string first;
  const std::size_t n = differing(lib::kCos, ncos, lm.cos, "cos", first) +
                        differing(lib::kSin, nsin, lm.sin, "sin", first) +
                        differing(lib::kExp, nexp, lm.exp, "exp", first);
  if (n == 0) {
    TADAH_SKIP("this C library gives every recorded result of cos, sin and "
               "exp (" << ncos + nsin + nexp << "): the golden cases run "
               "with it checked no other results");
  }
  SUCCEED(n << " of " << ncos + nsin + nexp
            << " recorded results differ; the first: " << first);
}

/* Pins: the frozen copy of the code before RHOFORM, run on this platform,
 * gives every stored golden value (1770 points, 4 values each) within
 * kLibmBound * DBL_EPSILON of the value's libm scale, and with glibc, the C
 * library of the capture, bit for bit. A copy whose formulas drift from
 * that code fails here on every platform, and with glibc any changed bit
 * fails; the cases below compare RHOFORM with the copy.
 * Before the fix: new; passes (a pin of the copy: the cases below compared
 * RHOFORM with the stored goldens and failed with musl, Alpine 3.19). */
TEST_CASE_METHOD(GoldenFixture,
                 "DM_REAM goldens: the density before RHOFORM computed here "
                 "is the stored goldens (bit for bit with glibc, within the "
                 "C library's rounding elsewhere)",
                 "[ream][rhoform][goldens]") {
  AnchorCount n;
  for (const auto &gc : ream_goldens::kCases)
    anchor(gc, before::platform_libm(), n);
  REQUIRE(n.values > 4 * 1500);  // the table is present and complete
  INFO((n.outside_bound ? "first value outside the bound: " + n.first_outside
                        : std::string("every value is within the bound")));
  CHECK(n.outside_bound == 0);
#if defined(__GLIBC__)
  INFO("glibc, the C library the goldens were captured with; first value "
       "that differs: " << n.first_not_bitwise);
  CHECK(n.not_bitwise == 0);
#endif
}

/* Pins: the anchor's bound holds for a C library whose cos, sin and exp
 * round one ulp away from this platform's (every result up, every result
 * down, or by the argument): the frozen copy through each such library
 * changes golden values (at least one value's bits) and keeps every value
 * within kLibmBound * DBL_EPSILON of its scale. The bound is what the
 * anchor asserts with a C library other than glibc (musl on the GitLab
 * runner); this case checks it against C libraries that round differently
 * on any platform.
 * Before the fix: new; passes (a pin of the bound). */
TEST_CASE_METHOD(GoldenFixture,
                 "DM_REAM goldens: the anchor's bound holds for a C library "
                 "whose cos, sin and exp round one ulp away",
                 "[ream][rhoform][goldens]") {
  const before::Libm libms[3] = {moved_libm<0>(), moved_libm<1>(),
                                 moved_libm<2>()};
  const char *names[3] = {"by the argument's lowest bits", "every result up",
                          "every result down"};
  for (int i = 0; i < 3; ++i) {
    INFO("cos, sin and exp one ulp away: " << names[i]);
    AnchorCount n;
    for (const auto &gc : ream_goldens::kCases) anchor(gc, libms[i], n);
    INFO((n.outside_bound ? "first value outside the bound: " + n.first_outside
                          : std::string("every value is within the bound")));
    CHECK(n.outside_bound == 0);
    CHECK(n.not_bitwise > 0);  // the library changes the density
  }
}

/* Pins: RHOFORM SIGNED gives the density before RHOFORM, computed here on
 * this platform by the frozen copy, bit for bit at every golden point, on
 * both amplitude sets.
 * Before the fix: passes as well (a pin of DM_REAM; the defect was this
 * test's reference, the stored goldens, which differ from the density of
 * any C library other than glibc: with musl, Alpine 3.19, it failed for
 * every Cut_Cos case and the Gaussian Cut_Poly cases, 'calc_rho at r =
 * 0x1.e60aa64c2f838p-1: got 0x1.4f8f2a2a6d567p-2, before RHOFORM
 * 0x1.4f8f2a2a6d566p-2'). */
TEST_CASE_METHOD(GoldenFixture,
                 "DM_REAM goldens: RHOFORM SIGNED is bitwise the density before RHOFORM "
                 "(signed and non-negative amplitudes)",
                 "[ream][rhoform][goldens]") {
  std::size_t npts = 0;
  for (const auto &gc : ream_goldens::kCases) {
    INFO(gc.family << " " << gc.cutoff
                   << (gc.nonneg_amps ? " non-negative" : " signed")
                   << " amplitudes");
    std::string first;
    CHECK(mismatches(gc, "SIGNED", first) == 0);
    INFO(first);
    CHECK(first.empty());
    npts += gc.npts;
  }
  REQUIRE(npts > 1500);  // the table is present and complete
}

/* Pins: RHOFORM NONNEG gives the density before RHOFORM, computed here on
 * this platform by the frozen copy, bit for bit at every golden point of
 * the 15 non-negative amplitude sets.
 * Before the fix: passes as well (a pin; with musl it failed against the
 * stored goldens for every Cut_Cos case and Gaussian Cut_Poly). */
TEST_CASE_METHOD(GoldenFixture,
                 "DM_REAM goldens: RHOFORM NONNEG is bitwise the density before RHOFORM "
                 "on non-negative amplitudes",
                 "[ream][rhoform][goldens]") {
  std::size_t ncases = 0;
  for (const auto &gc : ream_goldens::kCases) {
    if (!gc.nonneg_amps) continue;
    ++ncases;
    INFO(gc.family << " " << gc.cutoff);
    std::string first;
    CHECK(mismatches(gc, "NONNEG", first) == 0);
    INFO(first);
    CHECK(first.empty());
  }
  REQUIRE(ncases == 15);
}

TEST_CASE_METHOD(GoldenFixture,
                 "DM_REAM goldens: NONNEG refuses the signed amplitude sets",
                 "[ream][rhoform][goldens]") {
  for (const auto &gc : ream_goldens::kCases) {
    if (gc.nonneg_amps) continue;
    INFO(gc.family << " " << gc.cutoff);
    Context c = golden_ctx(gc, "NONNEG");
    REQUIRE_THROWS_AS(DM_REAM(c), DensityFormError);
  }
}

/* Pins: the density a program gets from ReamDensity is the code of
 * src/ream_density.cpp, compiled with Tadah.MODELS' flags, also when the
 * program calls ReamDensity itself (DM_REAM::density()) rather than
 * through DM_REAM: for every golden case under RHOFORM SIGNED, with f_c and
 * f_c' of a cutoff of the case's RCTYPEMB (RCUTMB 5.5) built here,
 * w * value(r, f_c) is calc_rho's increment before RHOFORM, and
 * value_deriv(r, f_c, f_c') times w gives eam_density's rho and rho'
 * before RHOFORM, bit for bit (w = 0.8125, golden_ctx's WATOM; the density
 * before RHOFORM computed here on this platform by the frozen copy). In
 * tadah.models_ream_goldens_native_link this unit is compiled with
 * -march=native: a ReamDensity method moved into ream_density.h would be
 * compiled here, with those flags, and fuse its multiply-adds into FMAs on
 * an FMA machine.
 * Before the fix: passes as well (a pin). With deg, width, value, deriv
 * and value_deriv moved into ream_density.h as inline functions it fails in
 * the native-link program at every golden case (GCC 15.2, LTO), and so do
 * the other golden cases there: the link keeps one copy of an inline
 * function, this unit's -march=native one, and DM_REAM's methods call it.
 * Without this case, whose unit calls the density itself, that move
 * passed. Its reference was the stored goldens, with which it failed with
 * musl (Alpine 3.19) for every Cut_Cos case and Gaussian Cut_Poly. */
TEST_CASE_METHOD(GoldenFixture,
                 "DM_REAM goldens: ReamDensity called by the program gives the "
                 "density before RHOFORM",
                 "[ream][rhoform][goldens]") {
  const double w = 0.8125;
  std::size_t npts = 0;
  for (const auto &gc : ream_goldens::kCases) {
    INFO(gc.family << " " << gc.cutoff
                   << (gc.nonneg_amps ? " non-negative" : " signed")
                   << " amplitudes");
    Context c = golden_ctx(gc, "SIGNED");
    const DM_REAM dm(c);
    const std::vector<before::Values> ref = before_rhoform(gc);
    double rcut = 5.5;
    const std::unique_ptr<Cut_Base> cut(
        factory<Cut_Base, double &>(gc.cutoff, rcut));
    std::size_t bad = 0;
    std::string first;
    for (std::size_t k = 0; k < gc.npts; ++k) {
      const double r = gc.pts[k].r;
      const double fc = cut->calc(r), fcp = cut->calc_prime(r);
      double v = 0.0, d = 0.0;
      dm.density().value_deriv(r, fc, fcp, v, d);
      const double got[3] = {w * dm.density().value(r, fc), w * v, w * d};
      const double want[3] = {ref[k].rho, ref[k].eam_rho, ref[k].eam_rhop};
      const char *what[3] = {"w value", "w value_deriv rho",
                             "w value_deriv rho'"};
      for (int q = 0; q < 3; ++q) {
        if (bits(got[q]) == bits(want[q])) continue;
        if (bad++ == 0)
          first = std::string(what[q]) + " at r = " + hex(r) + ": got " +
                  hex(got[q]) + ", before RHOFORM " + hex(want[q]);
      }
    }
    INFO(first);
    CHECK(bad == 0);
    npts += gc.npts;
  }
  REQUIRE(npts > 1500);
}
