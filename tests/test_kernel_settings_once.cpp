// test_kernel_settings_once.cpp
//
// The line a kernel prints when it is built at VERBOSE > 0 (its label and
// settings, "Kern_RBF | sigma: 1 | gamma: 1") is printed once per process
// for each kernel, as the descriptors' configuration warnings are logged
// (warn_once.h): tadah hpo builds the model of every evaluation anew, and
// the line was repeated at every evaluation.
//
// Each case records whether the line was already printed in this process
// (printed_once) before it builds anything: in a CTest entry, a process of
// its own, it was not, and the case requires the line exactly once; run in
// one process with other cases that build the same kernel at VERBOSE > 0,
// it requires none. Either way the line is printed in the process when the
// case ends, so a line that is never printed fails.
#include "catch2/catch.hpp"
#include <tadah/core/context.h>
#include <tadah/core/logger.h>
#include <tadah/models/descriptors/warn_once.h>
#include <tadah/models/functions/kernels/kern_lq.h>
#include <tadah/models/functions/kernels/kern_polynomial.h>
#include <tadah/models/functions/kernels/kern_quadratic.h>
#include <tadah/models/functions/kernels/kern_rbf.h>
#include <tadah/models/functions/kernels/kern_sigmoid.h>

#include <iostream>
#include <sstream>
#include <string>

using tadah::core::Context;

namespace {

/* Captures std::cout while it lives, and restores the Logger level that a
 * Context's VERBOSE sets (Context::add). */
class CoutCapture {
public:
  CoutCapture()
      : old_(std::cout.rdbuf(out_.rdbuf())),
        level_(tadah::core::Logger::instance().getLogLevel()) {}
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

/* A kernel context: VERBOSE @p verbose and, when @p mparams is not empty,
 * MPARAMS. */
Context kernel_ctx(int verbose, const std::string &mparams) {
  Context c;
  c.add("VERBOSE", verbose);
  if (!mparams.empty()) c.add("MPARAMS", mparams);
  return c;
}

/* Builds the kernel K once at VERBOSE 0, then three times at VERBOSE 1 from
 * the same MPARAMS, and checks what std::cout received: nothing from the
 * VERBOSE 0 build; @p line (the kernel's settings line) once when this
 * process had not printed the line of @p label before the case, else not
 * at all, and the label nowhere else; and the line printed in the process
 * now. */
template <class K>
void check_printed_once(const std::string &label, const std::string &mparams,
                        const std::string &line) {
  const bool before = tadah::models::printed_once(label);
  std::string quiet, verbose;
  {
    CoutCapture out;
    K k(kernel_ctx(0, mparams));
    quiet = out.str();
  }
  {
    CoutCapture out;
    for (int n = 0; n < 3; ++n) K k(kernel_ctx(1, mparams));
    verbose = out.str();
  }
  INFO("VERBOSE 0 output: '" << quiet << "'");
  INFO("VERBOSE 1 output: '" << verbose << "'");
  INFO("printed before this case: " << (before ? "yes" : "no"));
  CHECK(quiet.empty());
  CHECK(occurrences(verbose, line) == (before ? 0u : 1u));
  CHECK(occurrences(verbose, label) == (before ? 0u : 1u));
  CHECK(tadah::models::printed_once(label));
}

}  // namespace

/* Pins: a Kern_RBF built at VERBOSE 0 prints nothing, and three built at
 * VERBOSE 1 from MPARAMS 1.0 print "\nKern_RBF | sigma: 1 | gamma: 1\n"
 * exactly once between them (none when this process printed it before the
 * case), and the line is printed in the process afterwards
 * (check_printed_once).
 * Before the fix: three times, once per kernel object (tadah hpo printed it
 * once per evaluation); printed_once does not exist. */
TEST_CASE("Kern_RBF prints its settings once per process at VERBOSE > 0",
          "[kernel][print_once]") {
  check_printed_once<tadah::models::Kern_RBF>(
      "Kern_RBF", "1.0", "\nKern_RBF | sigma: 1 | gamma: 1\n");
}

/* Pins: as for Kern_RBF, Kern_Polynomial (MPARAMS 2 0.5 1.5): nothing at
 * VERBOSE 0, and "\nKern_Polynomial | degree: 2 | gamma: 0.5 | c: 1.5\n"
 * exactly once between three VERBOSE 1 kernels (none when this process
 * printed it before), printed in the process afterwards.
 * Before the fix: three times; printed_once does not exist. */
TEST_CASE("Kern_Polynomial prints its settings once per process at "
          "VERBOSE > 0",
          "[kernel][print_once]") {
  check_printed_once<tadah::models::Kern_Polynomial>(
      "Kern_Polynomial", "2 0.5 1.5",
      "\nKern_Polynomial | degree: 2 | gamma: 0.5 | c: 1.5\n");
}

/* Pins: as for Kern_RBF, Kern_Sigmoid (MPARAMS 0.25 -1): nothing at
 * VERBOSE 0, and "\nKern_Sigmoid | gamma: 0.25 | c: -1\n" exactly once
 * between three VERBOSE 1 kernels (none when this process printed it
 * before), printed in the process afterwards.
 * Before the fix: three times; printed_once does not exist. */
TEST_CASE("Kern_Sigmoid prints its settings once per process at VERBOSE > 0",
          "[kernel][print_once]") {
  check_printed_once<tadah::models::Kern_Sigmoid>(
      "Kern_Sigmoid", "0.25 -1", "\nKern_Sigmoid | gamma: 0.25 | c: -1\n");
}

/* Pins: as for Kern_RBF, Kern_Quadratic: nothing at VERBOSE 0, and
 * "\nKern_Quadratic\n" exactly once between three VERBOSE 1 kernels (none
 * when this process printed it before), printed in the process afterwards.
 * Before the fix: three times; printed_once does not exist. */
TEST_CASE("Kern_Quadratic prints its label once per process at VERBOSE > 0",
          "[kernel][print_once]") {
  check_printed_once<tadah::models::Kern_Quadratic>("Kern_Quadratic", "",
                                                     "\nKern_Quadratic\n");
}

/* Pins: as for Kern_RBF, Kern_LQ: nothing at VERBOSE 0, and "\nKern_LQ\n"
 * exactly once between three VERBOSE 1 kernels (none when this process
 * printed it before), printed in the process afterwards.
 * Before the fix: three times; printed_once does not exist. */
TEST_CASE("Kern_LQ prints its label once per process at VERBOSE > 0",
          "[kernel][print_once]") {
  check_printed_once<tadah::models::Kern_LQ>("Kern_LQ", "", "\nKern_LQ\n");
}
