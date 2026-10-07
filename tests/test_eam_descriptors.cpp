// D2_EAM and DM_EAM read their tabulated functions from one setfl file, named
// by the SETFL key (D2_EAM the pair function, DM_EAM the density and the
// embedding function). Each descriptor takes one TYPE2B / TYPEMB parameter,
// which the key formats require (TYPEMB admits only an unsigned integer
// there): it is a placeholder, written 0. Any other value has no effect and
// is reported on standard error once per process
// (core::report_ignored_key_once, the notice of every key without effect).
//
// The contexts below are the ones a configuration file gives, e.g.
//   TYPE2B D2_EAM 0 Ta Ta
//   TYPEMB DM_EAM 0 Ta Ta
//   SETFL  tests_data/eam.alloy

#include "catch2/catch.hpp"

#include <tadah/core/context.h>
#include <tadah/core/key_notice.h>
#include <tadah/core/periodic_table.h>
#include <tadah/models/cutoffs.h>
#include <tadah/models/descriptors/d2/d2_eam.h>
#include <tadah/models/descriptors/dm/dm_eam.h>
#include <tadah/models/eam/eam_setfl_io.h>

#include <iostream>
#include <memory>
#include <sstream>
#include <string>

using tadah::core::Context;
using tadah::models::Cut_Dummy;
using tadah::models::D2_EAM;
using tadah::models::DM_EAM;
using tadah::models::EAMSetflIO;

namespace {

const char *const kSetfl = "tests_data/eam.alloy";  // Ta, cutoff 5.58190742024
const double kSetflRcut = 5.58190742024;
const int kTa = 73;

/* Captures std::cerr (the notice) while it lives. */
class CerrCapture {
public:
  CerrCapture() : old_(std::cerr.rdbuf(err_.rdbuf())) {}
  ~CerrCapture() { std::cerr.rdbuf(old_); }
  CerrCapture(const CerrCapture &) = delete;
  CerrCapture &operator=(const CerrCapture &) = delete;
  std::string str() const { return err_.str(); }

private:
  std::ostringstream err_;
  std::streambuf *old_;
};

std::size_t count_of(const std::string &text, const std::string &needle) {
  std::size_t n = 0;
  for (std::size_t p = text.find(needle); p != std::string::npos;
       p = text.find(needle, p + needle.size()))
    ++n;
  return n;
}

/* A context with the D2_EAM pair term (`param` its TYPE2B parameter). */
Context d2_context(const std::string &param) {
  Context c;
  c.add("ATOM", "Ta");
  c.remove("INIT2B");
  c.add("INIT2B", "true");
  c.add("TYPE2B", "D2_EAM");
  c.add("TYPE2B", param);
  c.add("TYPE2B", "Ta");
  c.add("TYPE2B", "Ta");
  c.add("RCUT2B", kSetflRcut);
  c.add("SETFL", kSetfl);
  return c;
}

/* A context with the DM_EAM many-body term (`param` its TYPEMB parameter). */
Context dm_context(const std::string &param) {
  Context c;
  c.add("ATOM", "Ta");
  c.remove("INITMB");
  c.add("INITMB", "true");
  c.add("TYPEMB", "DM_EAM");
  c.add("TYPEMB", param);
  c.add("TYPEMB", "Ta");
  c.add("TYPEMB", "Ta");
  c.add("RCUTMB", kSetflRcut);
  c.add("SETFL", kSetfl);
  return c;
}

}  // namespace

/* Pins: one context with `TYPE2B D2_EAM 0 Ta Ta`, `TYPEMB DM_EAM 0 Ta Ta`
 * and one SETFL builds both descriptors (size 1 each), and both evaluate the
 * functions of that file: at r = 2.86 A (the bcc Ta nearest neighbour) for a
 * Ta-Ta pair, D2_EAM's descriptor is EAMSetflIO's phi(r) and DM_EAM's pair
 * density is EAMSetflIO's rho(r), for an EAMSetflIO that read the same file
 * (the same splines; 1e-12 relative).
 * Before the fix: DM_EAM read its setfl path from TYPEMB parameter 1, the
 * placeholder "0": "EAMSetflIO: setfl file does not exist: 0", so the
 * REQUIRE_NOTHROW for DM_EAM failed. */
TEST_CASE("DM_EAM reads its setfl file from SETFL, the file D2_EAM reads",
          "[DM][D2][eam]") {
  tadah::core::PeriodicTable::initialize();
  Context c = d2_context("0");
  c.remove("INITMB");
  c.add("INITMB", "true");
  c.add("TYPEMB", "DM_EAM");
  c.add("TYPEMB", "0");
  c.add("TYPEMB", "Ta");
  c.add("TYPEMB", "Ta");
  c.add("RCUTMB", kSetflRcut);

  EAMSetflIO eam;
  eam.read_setfl(kSetfl);
  eam.build_splines();
  REQUIRE(eam.rcut() == Approx(kSetflRcut).epsilon(1e-15));
  const double r = 2.86;

  std::unique_ptr<D2_EAM> d2;
  REQUIRE_NOTHROW(d2 = std::make_unique<D2_EAM>(c));
  Cut_Dummy fc2(kSetflRcut);
  d2->set_fcut(&fc2, false);
  CHECK(d2->size() == 1);
  tadah::core::aed_type aed(d2->size());
  aed.set_zero();
  d2->calc_aed(kTa, kTa, r, r * r, aed);
  CHECK(aed(0) == Approx(eam.phi_tab(r)).epsilon(1e-12));

  std::unique_ptr<DM_EAM> dm;
  REQUIRE_NOTHROW(dm = std::make_unique<DM_EAM>(c));
  Cut_Dummy fcm(kSetflRcut);
  dm->set_fcut(&fcm, false);
  CHECK(dm->size() == 1);
  tadah::core::rho_type rho;
  dm->init_rhoi(rho);
  dm->calc_rho(kTa, kTa, r, r * r, tadah::core::Vec3d(r, 0.0, 0.0), rho);
  CHECK(rho(0) != 0.0);
  CHECK(rho(0) == Approx(eam.rho_tab(r)).epsilon(1e-12));
}

/* Pins: a setfl path in DM_EAM's TYPEMB parameter, the form the old tests
 * used (`TYPEMB DM_EAM tests_data/eam.alloy Ta Ta`, which no configuration
 * file can hold: the TYPEMB format takes an unsigned integer there), is not
 * read: without SETFL the descriptor fails naming SETFL. A SETFL naming no
 * file fails naming that path.
 * Before the fix: the first context built a DM_EAM from the TYPEMB path
 * (no throw, so CHECK_THROWS_WITH failed), and the second failed naming the
 * placeholder ("setfl file does not exist: 0"), not path/to/nowhere. */
TEST_CASE("DM_EAM: a setfl path in the TYPEMB parameter is not read",
          "[DM][eam]") {
  tadah::core::PeriodicTable::initialize();
  Context no_setfl = dm_context(kSetfl);
  no_setfl.remove("SETFL");
  {
    CerrCapture cap;  // the placeholder notice for the path, if not yet given
    CHECK_THROWS_WITH(DM_EAM(no_setfl), Catch::Contains("SETFL"));
  }

  Context nowhere = dm_context("0");
  nowhere.remove("SETFL");
  nowhere.add("SETFL", "path/to/nowhere");
  CHECK_THROWS_WITH(DM_EAM(nowhere), Catch::Contains("path/to/nowhere"));
}

/* Pins: `TYPE2B D2_EAM 12345 Ta Ta` builds the descriptor and is reported
 * once per process on standard error, with the exact line
 *   WARNING: TYPE2B D2_EAM parameter '12345' is ignored: it is a
 *   placeholder (write 0); D2_EAM reads its setfl file from SETFL.
 * (one line), however often the descriptor is built; the placeholder 0 is
 * not reported. When an earlier test case in this process reported the
 * parameter already, nothing more may be written (as for EMBF0).
 * Before the fix: nothing was written (the parameter was dropped without a
 * word), so the count CHECK and ignored_key_reported failed. */
TEST_CASE("D2_EAM reports a TYPE2B parameter other than the placeholder 0 "
          "once",
          "[D2][eam][key_notice]") {
  tadah::core::PeriodicTable::initialize();
  const std::string key = "TYPE2B D2_EAM parameter";
  const bool before = tadah::core::ignored_key_reported(key);
  const std::string line =
      "WARNING: TYPE2B D2_EAM parameter '12345' is ignored: it is a "
      "placeholder (write 0); D2_EAM reads its setfl file from SETFL.\n";

  std::string zero, value;
  {
    CerrCapture cap;
    REQUIRE_NOTHROW(D2_EAM(d2_context("0")));
    zero = cap.str();
  }
  INFO("stderr with 0:\n" << zero);
  CHECK(count_of(zero, "is ignored") == 0);
  CHECK(tadah::core::ignored_key_reported(key) == before);
  {
    CerrCapture cap;
    REQUIRE_NOTHROW(D2_EAM(d2_context("12345")));
    REQUIRE_NOTHROW(D2_EAM(d2_context("12345")));
    value = cap.str();
  }
  INFO("stderr with 12345:\n" << value);
  CHECK(count_of(value, line) == (before ? 0u : 1u));
  CHECK(count_of(value, "is ignored") == (before ? 0u : 1u));
  CHECK(tadah::core::ignored_key_reported(key));
}

/* Pins: the same for DM_EAM: `TYPEMB DM_EAM 7 Ta Ta` gives, once per
 * process,
 *   WARNING: TYPEMB DM_EAM parameter '7' is ignored: it is a placeholder
 *   (write 0); DM_EAM reads its setfl file from SETFL.
 * and the placeholder 0 gives nothing.
 * Before the fix: DM_EAM read its setfl path from the parameter, so both
 * contexts failed ("setfl file does not exist: 0", then ": 7") and the
 * REQUIRE_NOTHROW failed. */
TEST_CASE("DM_EAM reports a TYPEMB parameter other than the placeholder 0 "
          "once",
          "[DM][eam][key_notice]") {
  tadah::core::PeriodicTable::initialize();
  const std::string key = "TYPEMB DM_EAM parameter";
  const bool before = tadah::core::ignored_key_reported(key);
  const std::string line =
      "WARNING: TYPEMB DM_EAM parameter '7' is ignored: it is a "
      "placeholder (write 0); DM_EAM reads its setfl file from SETFL.\n";

  std::string zero, value;
  {
    CerrCapture cap;
    REQUIRE_NOTHROW(DM_EAM(dm_context("0")));
    zero = cap.str();
  }
  INFO("stderr with 0:\n" << zero);
  CHECK(count_of(zero, "is ignored") == 0);
  CHECK(tadah::core::ignored_key_reported(key) == before);
  {
    CerrCapture cap;
    REQUIRE_NOTHROW(DM_EAM(dm_context("7")));
    REQUIRE_NOTHROW(DM_EAM(dm_context("7")));
    value = cap.str();
  }
  INFO("stderr with 7:\n" << value);
  CHECK(count_of(value, line) == (before ? 0u : 1u));
  CHECK(count_of(value, "is ignored") == (before ? 0u : 1u));
  CHECK(tadah::core::ignored_key_reported(key));
}

/* Pins: a D2_EAM built from a context with INIT2B false is empty (size 0)
 * and still takes its one TYPE2B parameter (nparams == 1), as DM_EAM with
 * INITMB false does.
 * Before the fix: D2_EAM(Context &) returned before init(), so nparams was
 * never set (an indeterminate value; the CHECK failed). */
TEST_CASE("D2_EAM built with INIT2B false takes its one parameter",
          "[D2][eam]") {
  tadah::core::PeriodicTable::initialize();
  Context off = d2_context("0");
  off.remove("INIT2B");
  off.add("INIT2B", "false");
  auto d2 = std::make_unique<D2_EAM>(off);
  CHECK(d2->size() == 0);
  CHECK(d2->nparams == 1);

  Context off_mb = dm_context("0");
  off_mb.remove("INITMB");
  off_mb.add("INITMB", "false");
  auto dm = std::make_unique<DM_EAM>(off_mb);
  CHECK(dm->size() == 0);
  CHECK(dm->nparams == 1);
}
