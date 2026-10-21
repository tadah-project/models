// test_rank_threshold.cpp
//
// The rank threshold of OALGO 1 and 3, OALGO's second value M, as the help
// of OALGO states it:
//
//  - OALGO 1 (DGELSD) leaves out of the solution a singular direction
//    whose singular value is at most M times the largest for 0 < M < 1,
//    and at most the machine precision of LAPACK, 2^-53, times the largest
//    for M < 0;
//  - OALGO 3 (the SVD solve) leaves out one whose singular value is at
//    most M times the largest for 0 < M < 1, and only one of singular
//    value 0 for M <= 0; its WARNING names that threshold.
//
// The other values of M are refused (test_oalgo_requests.cpp).
//
// The design matrix: two data rows, (1, 0) and (0, t), with targets 1 and
// 1, and the two zero ridge-augmentation rows of the DesignMatrix layout
// (LAMBDA 0). Its singular values are 1 and t, and its weights are
// (1, 1/t) with both directions kept, (1, 0) without the second.

#include "catch2/catch.hpp"

#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/core/logger.h>
#include <tadah/models/linear_regressor.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

using namespace tadah::core;
using namespace tadah::models;

namespace {

/* The weights of a LAMBDA 0 fit with OALGO "<oalgo> <M>" of the design
 * matrix above. */
t_type fit(int oalgo, const std::string &M, double t) {
  Context c;
  c.add("LAMBDA", "0");
  c.add("OALGO", std::to_string(oalgo) + " " + M);
  phi_type Phi(4, 2);
  Phi.set_zero();
  Phi(0, 0) = 1.0;
  Phi(1, 1) = t;
  t_type T(4);
  T[0] = 1.0;
  T[1] = 1.0;
  T[2] = 0.0;
  T[3] = 0.0;
  t_type w;
  LinearRegressor::train(c, Phi, T, w);
  return w;
}

struct Case {
  double t;
  const char *M;
  bool first_kept;
  bool second_kept;
};

void check(int oalgo, const Case &k) {
  DYNAMIC_SECTION("OALGO " << oalgo << " " << k.M << ", t " << k.t) {
    const t_type w = fit(oalgo, k.M, k.t);
    REQUIRE(w.size() == 2);
    if (k.first_kept)
      CHECK(w[0] == Approx(1.0).epsilon(1e-12));
    else
      CHECK(w[0] == 0.0);
    if (k.second_kept)
      CHECK(w[1] == Approx(1.0 / k.t).epsilon(1e-12));
    else
      CHECK(w[1] == 0.0);
  }
}

/* Captures std::cout (where the Logger prints) at WARNING level and
 * restores the stream and the level afterwards. */
class LogCapture {
public:
  LogCapture()
      : old_level_(Logger::instance().getLogLevel()),
        old_out_(std::cout.rdbuf(out_.rdbuf())) {
    Logger::instance().setLogLevel(LogLevel::WARNING);
  }
  ~LogCapture() {
    std::cout.rdbuf(old_out_);
    Logger::instance().setLogLevel(old_level_);
  }
  LogCapture(const LogCapture &) = delete;
  LogCapture &operator=(const LogCapture &) = delete;
  std::string text() const { return out_.str(); }

private:
  std::ostringstream out_;
  LogLevel old_level_;
  std::streambuf *old_out_;
};

}  // namespace

/* Pins: OALGO 1 with M 1e-8 or 0.5 leaves out the singular value 1e-12;
 * with M -1 it keeps it, and 1.5e-16, and leaves out 1e-16: DGELSD's
 * threshold is then 2^-53 (1.11e-16) times the largest.
 * A guard (DGELSD's rule, which the help of OALGO states): a fit that
 * passed DGELSD an M below 0 as the smallest positive double fails its
 * case t 1e-16 with M -1. */
TEST_CASE("OALGO 1: an M below 0 is the machine precision",
          "[oalgo1][rank_threshold]") {
  for (const Case &k :
       {Case{1e-12, "1e-8", true, false}, Case{1e-12, "0.5", true, false},
        Case{1e-6, "1e-8", true, true}, Case{1e-12, "-1", true, true},
        Case{1.5e-16, "-1", true, true}, Case{1e-16, "-1", true, false}})
    check(1, k);
}

/* Pins: OALGO 3 with M -1 or 0 keeps the singular value 1e-17, below the
 * machine precision times the largest (the coefficient 1e17), and with M
 * 1e-8 or 0.5 leaves it out.
 * A guard (the rule of OLS::pseudoinverse_threshold, which the help of
 * OALGO states): a threshold of 2^-53 times the largest for an M of 0 or
 * less fails its cases t 1e-17 with M -1 and 0. */
TEST_CASE("OALGO 3: an M of 0 or less keeps every singular value above 0",
          "[oalgo3][rank_threshold]") {
  for (const Case &k :
       {Case{1e-17, "-1", true, true}, Case{1e-17, "0", true, true},
        Case{1e-12, "-1", true, true}, Case{1e-17, "1e-8", true, false},
        Case{1e-12, "1e-8", true, false}, Case{1e-6, "1e-8", true, true},
        Case{0.2, "0.5", true, false}})
    check(3, k);
}

/* Pins: OALGO 3 with M -1 or 0 and an exactly zero column (t 0) leaves out
 * its direction, and the one WARNING names the threshold 0:
 * "LinearRegressor (OALGO 3): 1 of 2 singular directions of the design
 * matrix are at or below the rank threshold (0, as OALGO's rcond -1 is not
 * above 0): the solve leaves them out of the weights."
 * Before the fix: the WARNING named the threshold "OALGO's rcond -1 times
 * the largest singular value", a negative threshold the solve does not
 * use. */
TEST_CASE("OALGO 3: the WARNING of an M of 0 or less names the threshold 0",
          "[oalgo3][rank_threshold]") {
  for (const std::string M : {"-1", "0"}) {
    DYNAMIC_SECTION("OALGO 3 " << M) {
      t_type w;
      std::string log;
      {
        LogCapture cap;
        w = fit(3, M, 0.0);
        log = cap.text();
      }
      INFO(log);
      REQUIRE(w.size() == 2);
      CHECK(w[0] == Approx(1.0).epsilon(1e-12));
      CHECK(w[1] == 0.0);
      const std::string warning =
          "WARNING: LinearRegressor (OALGO 3): 1 of 2 singular directions "
          "of the design matrix are at or below the rank threshold (0, as "
          "OALGO's rcond " + M + " is not above 0): the solve leaves them "
          "out of the weights.\n";
      REQUIRE(log.size() >= warning.size());
      CHECK(log.compare(log.size() - warning.size(), warning.size(),
                        warning) == 0);
      CHECK(log.find("WARNING") == log.rfind("WARNING"));
    }
  }
}
