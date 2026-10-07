// test_descriptor_warnings_once.cpp
//
// The warnings a descriptor gives about its configuration when it is built
// (a knot past the cutoff, a bump whose support reaches past a Cut_Dummy
// cutoff, a non-compact basis with a Cut_Dummy cutoff) are logged once per
// process, as the notice of a key without effect is: tadah hpo builds the
// descriptors of every evaluation anew, and a warning per descriptor object
// was repeated at every evaluation.
//
// Each CTest entry is a process of its own, so there each warning below is
// logged exactly once; run in one process with other cases that build the
// same configuration it may already have been logged (0 here), never more
// than once.

#include "catch2/catch.hpp"
#include <tadah/core/context.h>
#include <tadah/core/logger.h>
#include <tadah/core/periodic_table.h>
#include <tadah/models/dc_selector.h>
#include <tadah/models/descriptors/d2/d2_blip.h>
#include <tadah/models/descriptors/d2/d2_knot5.h>

#include <iostream>
#include <sstream>
#include <string>

using tadah::core::Context;

namespace {

/* Captures std::cout, where the Logger writes, while it lives, with the
 * Logger at WARNING (the level of the tadah CLI). */
class CoutCapture {
public:
  CoutCapture()
      : old_(std::cout.rdbuf(out_.rdbuf())),
        level_(tadah::core::Logger::instance().getLogLevel()) {
    tadah::core::Logger::instance().setLogLevel(
        tadah::core::LogLevel::WARNING);
  }
  ~CoutCapture() {
    tadah::core::Logger::instance().setLogLevel(level_);
    std::cout.rdbuf(old_);
  }
  CoutCapture(const CoutCapture &) = delete;
  CoutCapture &operator=(const CoutCapture &) = delete;
  std::string str() const { return out_.str(); }

private:
  std::ostringstream out_;
  std::streambuf *old_;
  tadah::core::LogLevel level_;
};

std::size_t occurrences(const std::string &text, const std::string &needle) {
  std::size_t n = 0;
  for (std::size_t p = text.find(needle); p != std::string::npos;
       p = text.find(needle, p + needle.size()))
    ++n;
  return n;
}

/* A two-body H-H context: TYPE2B <label> <sizes...> H H, the grids given,
 * RCTYPE2B and RCUT2B. */
Context two_body(const std::string &label, const std::string &sizes,
                 const std::string &cgrid, const std::string &sgrid,
                 const std::string &rctype, double rcut) {
  Context c;
  c.add("ATOM", "H");
  c.add("WATOM", 1);
  c.remove("INIT2B");
  c.add("INIT2B", "true");
  c.add("VERBOSE", 0);
  c.add("CGRID2B", cgrid);
  if (!sgrid.empty()) c.add("SGRID2B", sgrid);
  c.remove("RCTYPE2B");
  c.add("RCTYPE2B", rctype);
  c.add("RCUT2B", rcut);
  c.add("TYPE2B", label + " " + sizes + " H H");
  return c;
}

}  // namespace

/* Pins: three D2_Knot5 descriptors built from one configuration whose
 * largest CGRID2B knot (4.5) lies past RCUT2B (4.0) log the warning "the
 * largest CGRID2B knot ... lies beyond RCUT2B" at most once between them.
 * Before the fix: once per descriptor object, three times (and once per
 * evaluation of a tadah hpo run that builds a model per evaluation). */
TEST_CASE("D2_Knot5 warns once per process of a knot past RCUT2B",
          "[descriptor][warn_once]") {
  tadah::core::PeriodicTable::initialize();
  const Context c =
      two_body("D2_Knot5", "3 0", "1.0 2.0 4.5", "", "Cut_Cos", 4.0);
  CoutCapture out;
  for (int k = 0; k < 3; ++k) tadah::models::D2_Knot5 d(c);
  INFO(out.str());
  CHECK(occurrences(out.str(), "lies beyond RCUT2B") <= 1);
}

/* Pins: three D2_Blip descriptors built from one configuration with a
 * Cut_Dummy cutoff and a bump whose support (centre 3.5 + 2/eta = 4.5)
 * reaches past RCUT2B (4.0) log the warning "has support reaching ...,
 * past RCUT2B" at most once between them.
 * Before the fix: three times, once per descriptor object. */
TEST_CASE("D2_Blip warns once per process of a bump past a Cut_Dummy "
          "RCUT2B",
          "[descriptor][warn_once]") {
  tadah::core::PeriodicTable::initialize();
  const Context c =
      two_body("D2_Blip", "2 2", "1.0 3.5", "2.0 2.0", "Cut_Dummy", 4.0);
  CoutCapture out;
  for (int k = 0; k < 3; ++k) tadah::models::D2_Blip d(c);
  INFO(out.str());
  CHECK(occurrences(out.str(), "has support reaching") <= 1);
}

/* Pins: three DC_Selectors (the descriptors of three models) built from one
 * configuration with the Gaussian basis of D2_BP and a Cut_Dummy cutoff
 * log the warning "uses a non-compact Gaussian basis with a discontinuous
 * (Cut_Dummy) cutoff" at most once between them.
 * Before the fix: three times, once per DC_Selector. */
TEST_CASE("A non-compact basis with Cut_Dummy is warned of once per "
          "process",
          "[descriptor][warn_once]") {
  tadah::core::PeriodicTable::initialize();
  Context c =
      two_body("D2_BP", "2 2", "1.0 2.0", "1.0 1.0", "Cut_Dummy", 4.0);
  c.remove("INITMB");
  c.add("INITMB", "false");
  CoutCapture out;
  for (int k = 0; k < 3; ++k) tadah::models::DC_Selector s(c);
  INFO(out.str());
  CHECK(occurrences(out.str(), "uses a non-compact") <= 1);
}
