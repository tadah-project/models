// test_model_sigma.cpp
//
// SIGMA, the weight covariance an OALGO 3 fit writes into the model
// Context (LinearRegressor::config_add_sigma) and get_param_file() copies
// into the potential, belongs to the fit that computed it: a fit replaces
// it, a fit that computes none (another OALGO, every weight pinned) leaves
// none, and read_sigma refuses a SIGMA whose count is not N*N+1.
//
// Also: an M_KRR model built from a Context with WEIGHTS (a potential used
// as the configuration) can be trained, as an M_BLR model can.

#include "catch2/catch.hpp"
#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/models/functions/basis_functions/bf_linear.h>
#include <tadah/models/functions/kernels/kern_linear.h>
#include <tadah/models/linear_regressor.h>
#include <tadah/models/m_blr_train.h>
#include <tadah/models/m_krr_train.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

using namespace tadah::core;
using namespace tadah::models;

namespace {

std::uint64_t bits(double v) {
  std::uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return u;
}

// Deterministic, well-conditioned 8x4 design matrix (as test_fixed_weights).
phi_type make_phi() {
  phi_type Phi(8, 4);
  double v = 0.17;
  for (size_t j = 0; j < Phi.cols(); ++j)
    for (size_t i = 0; i < Phi.rows(); ++i) {
      v = std::fmod(v * 1.618 + 0.31 * static_cast<double>(j + 1), 1.0) - 0.5;
      Phi(i, j) = v + (i == j ? 2.0 : 0.0);
    }
  return Phi;
}

t_type make_target(const phi_type &Phi) {
  t_type w_true(Phi.cols());
  for (size_t j = 0; j < Phi.cols(); ++j)
    w_true[j] = 0.5 * static_cast<double>(j + 1) * (j % 2 ? -1.0 : 1.0);
  t_type T = Phi * w_true;
  for (size_t i = 0; i < T.size(); ++i) T[i] += 0.01 * static_cast<double>(i % 3);
  return T;
}

Context make_ctx(const std::string &oalgo) {
  Context c;
  c.add("VERBOSE", "0");
  c.add("LAMBDA", "0");
  c.add("OALGO", oalgo);
  return c;
}

/* A SIGMA as an earlier OALGO 3 fit of a 4-weight model leaves it: N = 4
 * and 16 values. */
void add_earlier_sigma(Context &c) {
  c.add("SIGMA", "4");
  for (int k = 0; k < 16; ++k) c.add("SIGMA", 0.25 * (k + 1));
}

}  // namespace

/* Pins: two OALGO 3 fits (LinearRegressor::train) on one Context leave one
 * SIGMA, N*N+1 = 17 values for 4 weights: N = 4, then the 16 values of the
 * covariance the second fit returned, column-major, bit for bit.
 * Before the fix: 34 values (17 + 17): config_add_sigma appended SIGMA to
 * the earlier one (Context::add appends to a multi-value key), so a
 * potential trained from a potential, or HPO's reused trainer, carried
 * every earlier SIGMA. */
TEST_CASE("OALGO 3: a second fit on one Context replaces SIGMA", "[sigma]") {
  Context c = make_ctx("3");
  const phi_type Phi0 = make_phi();
  const t_type T0 = make_target(Phi0);
  Matrix Sigma;
  for (int fit = 0; fit < 2; ++fit) {
    phi_type Phi = Phi0;
    t_type T = T0;
    t_type w;
    LinearRegressor::train(c, Phi, T, w, Sigma);
  }
  REQUIRE(c.size("SIGMA") == 17);
  CHECK(c.get<size_t>("SIGMA") == 4);
  for (size_t col = 0; col < 4; ++col)
    for (size_t row = 0; row < 4; ++row)
      CHECK(bits(c.get<double>("SIGMA", 1 + 4 * col + row)) ==
            bits(Sigma(row, col)));
}

/* Pins: a fit that computes no covariance writes no SIGMA, whatever SIGMA
 * the Context carried from an earlier fit (the 17 values of a 4-weight
 * OALGO 3 fit): M_BLR_Train<BF_Linear> and M_KRR_Train<Kern_Linear> with
 * OALGO 1, with no FIXINDEX, with FIXINDEX 2 (one weight pinned) and with
 * FIXINDEX 1-4 (every weight pinned), and M_BLR_Train with OALGO 3 and
 * every weight pinned (the branch that solves nothing); with every weight
 * pinned the weights are the FIXWEIGHT values. An OALGO 3 fit that solves
 * writes one SIGMA of N*N+1 values. (M_KRR_Train with every weight pinned
 * and OALGO 3 solves on zeroed columns, and the covariance of such fits is
 * the subject of another change; it is not run here.)
 * Before the fix: every fit that computes no covariance kept the earlier 17
 * values, so the potential paired them with weights they were not computed
 * for (and `predict --error` used them); the OALGO 3 fits that solve
 * appended their SIGMA to the earlier one. */
TEST_CASE("A fit replaces or removes the SIGMA of an earlier fit",
          "[sigma]") {
  struct Pins {
    const char *name;
    const char *fixindex;   // nullptr: no FIXINDEX
    std::vector<double> fixweight;
  };
  const Pins pins[] = {{"no FIXINDEX", nullptr, {}},
                       {"FIXINDEX 2", "2", {0.7}},
                       {"FIXINDEX 1-4", "1-4", {0.1, 0.2, 0.3, 0.4}}};
  for (const char *oalgo : {"1", "3"}) {
    for (const Pins &p : pins) {
      for (const char *model : {"M_BLR BF_Linear", "M_KRR Kern_Linear"}) {
        const bool all_pinned = p.fixindex && std::string(p.fixindex) == "1-4";
        if (all_pinned && std::string(oalgo) == "3" &&
            std::string(model) == "M_KRR Kern_Linear")
          continue;
        DYNAMIC_SECTION(model << ", OALGO " << oalgo << ", " << p.name) {
          Context c = make_ctx(oalgo);
          if (p.fixindex) {
            c.add("FIXINDEX", p.fixindex);
            for (double w : p.fixweight) c.add("FIXWEIGHT", w);
          }
          add_earlier_sigma(c);
          phi_type Phi = make_phi();
          t_type T = make_target(Phi);
          t_type w;
          if (std::string(model) == "M_BLR BF_Linear") {
            BF_Linear bf;
            M_BLR_Train<BF_Linear> m(bf, c);
            REQUIRE_NOTHROW(m.train(Phi, T));
            w = m.get_weights();
          } else {
            Kern_Linear k;
            M_KRR_Train<Kern_Linear> m(k, c);
            REQUIRE_NOTHROW(m.train(Phi, T));
            w = m.get_weights();
          }
          REQUIRE(w.size() == 4);
          if (all_pinned)
            for (size_t i = 0; i < 4; ++i) CHECK(w[i] == p.fixweight[i]);
          const bool computed = std::string(oalgo) == "3" && !all_pinned;
          if (computed) {
            // The covariance of this fit (its N is the fitted columns', the
            // dimension under FIXINDEX is another lane's subject).
            CHECK(c.contains("SIGMA"));
            const size_t n = c.get<size_t>("SIGMA");
            CHECK(c.size("SIGMA") == n * n + 1);
          } else {
            CHECK_FALSE(c.has_stored("SIGMA"));
          }
        }
      }
    }
  }
}

/* Pins: read_sigma refuses a SIGMA whose value count is not N*N+1, naming
 * SIGMA, N, the count it needs and the count it has: N = 2 with 3 values
 * after it ("2 1 0 0", 4 of the 5 it needs) and N = 4 with 33 (two
 * appended SIGMAs of a 4-weight model, 34 values); "2 1 0 0 1" is read as
 * the 2x2 identity.
 * Before the fix: "2 1 0 0" was read with a 0 for the missing value, no
 * error; the 34 values were refused by Context::get's size check (CORE),
 * whose message names neither N nor the count SIGMA needs. */
TEST_CASE("read_sigma refuses a SIGMA whose count is not N*N+1", "[sigma]") {
  Context c;
  Matrix S;
  c.add("SIGMA", "2 1 0 0");
  std::string msg;
  try {
    LinearRegressor::read_sigma(c, S);
  } catch (const std::exception &e) {
    msg = e.what();
  }
  INFO("message: " << msg);
  CHECK(msg.find("SIGMA") != std::string::npos);
  CHECK(msg.find("N = 2") != std::string::npos);
  CHECK(msg.find("5 values") != std::string::npos);
  CHECK(msg.find("has 4") != std::string::npos);

  c.remove("SIGMA");
  c.add("SIGMA", "4");
  for (int k = 0; k < 33; ++k) c.add("SIGMA", 0.5);
  msg.clear();
  try {
    LinearRegressor::read_sigma(c, S);
  } catch (const std::exception &e) {
    msg = e.what();
  }
  INFO("message: " << msg);
  CHECK(msg.find("N = 4") != std::string::npos);
  CHECK(msg.find("17 values") != std::string::npos);
  CHECK(msg.find("has 34") != std::string::npos);

  c.remove("SIGMA");
  c.add("SIGMA", "2 1 0 0 1");
  REQUIRE_NOTHROW(LinearRegressor::read_sigma(c, S));
  REQUIRE(S.rows() == 2);
  REQUIRE(S.cols() == 2);
  CHECK(S(0, 0) == 1.0);
  CHECK(S(1, 0) == 0.0);
  CHECK(S(0, 1) == 0.0);
  CHECK(S(1, 1) == 1.0);
}

/* Pins: an M_KRR model built from a Context that has WEIGHTS (a potential
 * used as the configuration, `tadah train -c pot.tadah`, or an HPO -c
 * potential) holds those weights but is not trained, and train() fits it:
 * the weights are those of the same fit on a Context without WEIGHTS, bit
 * for bit. The same for M_BLR (which already behaved so).
 * Before the fix: M_KRR_Core::init set trained when WEIGHTS were present,
 * so train() threw "This object is already trained!". */
TEST_CASE("An M_KRR model built with WEIGHTS can be trained",
          "[krr][weights]") {
  const phi_type Phi0 = make_phi();
  const t_type T0 = make_target(Phi0);
  t_type w_ref;
  {
    Context c = make_ctx("1");
    Kern_Linear k;
    M_KRR_Train<Kern_Linear> m(k, c);
    phi_type Phi = Phi0;
    t_type T = T0;
    m.train(Phi, T);
    w_ref = m.get_weights();
  }
  REQUIRE(w_ref.size() == 4);
  for (const char *model : {"M_KRR Kern_Linear", "M_BLR BF_Linear"}) {
    DYNAMIC_SECTION(model) {
      Context c = make_ctx("1");
      c.add("WEIGHTS", "9 8 7 6");
      phi_type Phi = Phi0;
      t_type T = T0;
      t_type w;
      bool trained_before = true;
      std::string msg;
      auto fit = [&](auto &m) {
        trained_before = m.is_trained();
        CHECK(m.get_weights().size() == 4);
        try {
          m.train(Phi, T);
        } catch (const std::exception &e) {
          msg = e.what();
        }
        w = m.get_weights();
      };
      if (std::string(model) == "M_KRR Kern_Linear") {
        Kern_Linear k;
        M_KRR_Train<Kern_Linear> m(k, c);
        fit(m);
      } else {
        BF_Linear bf;
        M_BLR_Train<BF_Linear> m(bf, c);
        fit(m);
      }
      INFO("train: " << msg);
      CHECK_FALSE(trained_before);
      CHECK(msg.empty());
      REQUIRE(w.size() == 4);
      for (size_t i = 0; i < 4; ++i) CHECK(bits(w[i]) == bits(w_ref[i]));
    }
  }
}
