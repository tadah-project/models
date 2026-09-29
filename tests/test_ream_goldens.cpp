/* FEAT-02 goldens: the DM_REAM density paths against the bits of tadah 1.5.1.
 *
 * ream_goldens_data.h holds, for 5 families (Blip, Blip5, Knot5, Knot3,
 * Gaussian) x 3 cutoffs (Cut_Dummy, Cut_Cos, Cut_Poly) x 2 amplitude sets
 * (signed, all non-negative), calc_rho's increment, the calc_dXijdri force
 * prefactor and eam_density's rho and rho' at 56-64 distances in (0, rcut],
 * captured as hex floats from tadah 1.5.1 (the tree before the density
 * forms: MODELS 4c4144b). RHOFORM SIGNED must reproduce every value bit for
 * bit on both amplitude sets, and RHOFORM NONNEG on the non-negative set
 * (every cutoff here is >= 0, so the max(f_c, 0) clamp is inert).
 */
#include "catch2/catch.hpp"
#include "ream_goldens_data.h"

#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/core/periodic_table.h>
#include <tadah/models/descriptors/dm/dm_ream.h>

#include <cstdint>
#include <cstring>
#include <sstream>
#include <string>

using namespace tadah::core;
using namespace tadah::models;

namespace {

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

// Every golden value of one case through a DM_REAM in @p form; returns the
// number of values whose bits differ and describes the first one.
std::size_t mismatches(const ream_goldens::Case &gc, const char *form,
                       std::string &first) {
  Context c = golden_ctx(gc, form);
  DM_REAM dm(c);
  const int Z = PeriodicTable::find_by_symbol("Ta").Z;
  const Vec3d unit(1.0, 0.0, 0.0);
  std::size_t bad = 0;
  for (std::size_t k = 0; k < gc.npts; ++k) {
    const ream_goldens::Point &p = gc.pts[k];
    const double r = p.r;

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
    const double want[4] = {p.rho, p.pref, p.eam_rho, p.eam_rhop};
    const char *what[4] = {"calc_rho", "calc_dXijdri pref", "eam_density rho",
                           "eam_density rho'"};
    for (int q = 0; q < 4; ++q) {
      if (bits(got[q]) == bits(want[q])) continue;
      if (bad++ == 0)
        first = std::string(what[q]) + " at r = " + hex(r) + ": got " +
                hex(got[q]) + ", tadah 1.5.1 " + hex(want[q]);
    }
  }
  return bad;
}

}  // namespace

TEST_CASE_METHOD(GoldenFixture,
                 "DM_REAM goldens: RHOFORM SIGNED is bitwise the tadah 1.5.1 density "
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

TEST_CASE_METHOD(GoldenFixture,
                 "DM_REAM goldens: RHOFORM NONNEG is bitwise the tadah 1.5.1 density "
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
