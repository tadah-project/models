// test_fixed_weights.cpp
//
// FIXINDEX/FIXWEIGHT coverage for M_BLR_Train / M_KRR_Train: partial fixing
// (non-prefix, unsorted index sets), the all-fixed scatter, and the
// lone-key/out-of-bounds error paths. The partial-fix reference is the
// conditional least-squares solution computed directly on the free columns
// with the fixed contribution subtracted from the target.

#include "catch2/catch.hpp"
#include <tadah/models/m_blr_train.h>
#include <tadah/models/m_krr_train.h>
#include <tadah/models/functions/basis_functions/bf_linear.h>
#include <tadah/models/functions/kernels/kern_linear.h>
#include <tadah/models/ols.h>
#include <tadah/core/context.h>
#include <tadah/core/core_types.h>

#include <cmath>
#include <string>
#include <vector>

using namespace tadah::core;
using namespace tadah::models;

namespace {

// Deterministic, well-conditioned 8x4 design matrix.
phi_type make_phi() {
  phi_type Phi(8, 4);
  double v = 0.17;
  for (size_t j = 0; j < Phi.cols(); ++j)
    for (size_t i = 0; i < Phi.rows(); ++i) {
      v = std::fmod(v * 1.618 + 0.31 * static_cast<double>(j + 1), 1.0) - 0.5;
      Phi(i, j) = v + (i == j ? 2.0 : 0.0);  // diagonal boost: full rank
    }
  return Phi;
}

t_type make_target(const phi_type &Phi) {
  t_type w_true(Phi.cols());
  for (size_t j = 0; j < Phi.cols(); ++j)
    w_true[j] = 0.5 * static_cast<double>(j + 1) * (j % 2 ? -1.0 : 1.0);
  t_type T = Phi * w_true;
  for (size_t i = 0; i < T.size(); ++i)  // perturb so fixing changes the fit
    T[i] += 0.01 * static_cast<double>(i % 3);
  return T;
}

Context make_ctx() {
  Context c;
  c.add("VERBOSE", "0");
  c.add("LAMBDA", "0");
  c.add("OALGO", "1");  // GELSD, rcond 1e-8 (same route the trainer takes)
  return c;
}

// Conditional LS reference: OLS on the free columns of Phi against
// T - Phi_fixed * w_fixed. fix_idx is 1-based, paired with w_fix.
t_type reference_free_weights(const phi_type &Phi, const t_type &T,
                              const std::vector<size_t> &fix_idx,
                              const std::vector<double> &w_fix) {
  std::vector<size_t> free_cols;
  for (size_t j = 0; j < Phi.cols(); ++j) {
    bool fixed = false;
    for (size_t k : fix_idx) fixed |= (k - 1 == j);
    if (!fixed) free_cols.push_back(j);
  }
  phi_type Phi_r(Phi.rows(), free_cols.size());
  for (size_t j = 0; j < free_cols.size(); ++j)
    for (size_t i = 0; i < Phi.rows(); ++i)
      Phi_r(i, j) = Phi(i, free_cols[j]);
  t_type T_r = T;
  for (size_t k = 0; k < fix_idx.size(); ++k)
    for (size_t i = 0; i < Phi.rows(); ++i)
      T_r[i] -= Phi(i, fix_idx[k] - 1) * w_fix[k];

  t_type w(free_cols.size());
  OLS ols;
  ols.solve(Phi_r, T_r, w, 1e-8, Algorithm::GELSD);
  return w;
}

}  // namespace

TEST_CASE("BLR partial fix: single non-prefix index", "[fixed_weights]") {
  // Fixing a column that is NOT a prefix of {1..k} exercises the
  // move_map-based reassembly; before the two-scatter fix this read w_f out
  // of bounds and produced garbage in the fixed slot.
  Context c = make_ctx();
  c.add("FIXINDEX", "2");
  c.add("FIXWEIGHT", "0.7");

  const phi_type Phi0 = make_phi();
  const t_type T0 = make_target(Phi0);

  BF_Linear bf;
  M_BLR_Train<BF_Linear> m(bf, c);
  phi_type Phi = Phi0;
  t_type T = T0;
  m.train(Phi, T);

  const t_type &w = m.get_weights();
  REQUIRE(w.size() == 4);
  REQUIRE(w[1] == 0.7);  // exact: fixed value must be placed verbatim

  const t_type wref = reference_free_weights(Phi0, T0, {2}, {0.7});
  REQUIRE(std::abs(w[0] - wref[0]) < 1e-10);
  REQUIRE(std::abs(w[2] - wref[1]) < 1e-10);
  REQUIRE(std::abs(w[3] - wref[2]) < 1e-10);
}

TEST_CASE("BLR partial fix: unsorted non-prefix pair", "[fixed_weights]") {
  // FIXINDEX order defines the FIXWEIGHT pairing: 4 -> 0.3, 1 -> -1.1.
  Context c = make_ctx();
  c.add("FIXINDEX", "4,1");
  c.add("FIXWEIGHT", "0.3");
  c.add("FIXWEIGHT", "-1.1");

  const phi_type Phi0 = make_phi();
  const t_type T0 = make_target(Phi0);

  BF_Linear bf;
  M_BLR_Train<BF_Linear> m(bf, c);
  phi_type Phi = Phi0;
  t_type T = T0;
  m.train(Phi, T);

  const t_type &w = m.get_weights();
  REQUIRE(w.size() == 4);
  REQUIRE(w[3] == 0.3);
  REQUIRE(w[0] == -1.1);

  const t_type wref = reference_free_weights(Phi0, T0, {4, 1}, {0.3, -1.1});
  REQUIRE(std::abs(w[1] - wref[0]) < 1e-10);
  REQUIRE(std::abs(w[2] - wref[1]) < 1e-10);
}

TEST_CASE("BLR all fixed: permuted scatter, no LS solve", "[fixed_weights]") {
  Context c = make_ctx();
  c.add("FIXINDEX", "3,1,2");
  c.add("FIXWEIGHT", "10.0");
  c.add("FIXWEIGHT", "20.0");
  c.add("FIXWEIGHT", "30.0");

  phi_type Phi(8, 3);
  for (size_t j = 0; j < Phi.cols(); ++j)
    for (size_t i = 0; i < Phi.rows(); ++i)
      Phi(i, j) = static_cast<double>(i + j + 1);
  t_type T(8);
  for (size_t i = 0; i < T.size(); ++i) T[i] = 1.0;

  BF_Linear bf;
  M_BLR_Train<BF_Linear> m(bf, c);
  m.train(Phi, T);

  const t_type &w = m.get_weights();
  REQUIRE(w.size() == 3);
  REQUIRE(w[2] == 10.0);  // FIXINDEX[0] = 3
  REQUIRE(w[0] == 20.0);  // FIXINDEX[1] = 1
  REQUIRE(w[1] == 30.0);  // FIXINDEX[2] = 2
  REQUIRE(m.is_trained());
}

TEST_CASE("BLR partial fix: prediction consistency", "[fixed_weights]") {
  // Phi0 * w must equal Phi_free * w_free + Phi_fixed * w_fixed; this holds
  // iff every weight landed in its original column slot.
  Context c = make_ctx();
  c.add("FIXINDEX", "3");
  c.add("FIXWEIGHT", "0.25");

  const phi_type Phi0 = make_phi();
  const t_type T0 = make_target(Phi0);

  BF_Linear bf;
  M_BLR_Train<BF_Linear> m(bf, c);
  phi_type Phi = Phi0;
  t_type T = T0;
  m.train(Phi, T);

  const t_type &w = m.get_weights();
  const t_type wref = reference_free_weights(Phi0, T0, {3}, {0.25});
  const t_type pred = Phi0 * w;
  for (size_t i = 0; i < pred.size(); ++i) {
    double expected = Phi0(i, 2) * 0.25;
    expected += Phi0(i, 0) * wref[0] + Phi0(i, 1) * wref[1] +
                Phi0(i, 3) * wref[2];
    REQUIRE(std::abs(pred[i] - expected) < 1e-10);
  }
}

TEST_CASE("BLR lone FIXINDEX/FIXWEIGHT throws", "[fixed_weights]") {
  const phi_type Phi0 = make_phi();
  const t_type T0 = make_target(Phi0);
  BF_Linear bf;

  SECTION("lone FIXINDEX") {
    Context c = make_ctx();
    c.add("FIXINDEX", "2");
    M_BLR_Train<BF_Linear> m(bf, c);
    phi_type Phi = Phi0;
    t_type T = T0;
    REQUIRE_THROWS_WITH(m.train(Phi, T), Catch::Contains("requires"));
  }
  SECTION("lone FIXWEIGHT") {
    Context c = make_ctx();
    c.add("FIXWEIGHT", "0.7");
    M_BLR_Train<BF_Linear> m(bf, c);
    phi_type Phi = Phi0;
    t_type T = T0;
    REQUIRE_THROWS_WITH(m.train(Phi, T), Catch::Contains("requires"));
  }
}

TEST_CASE("BLR FIXINDEX out of bounds throws", "[fixed_weights]") {
  Context c = make_ctx();
  c.add("FIXINDEX", "7");
  c.add("FIXWEIGHT", "0.7");

  phi_type Phi = make_phi();
  t_type T = make_target(Phi);

  BF_Linear bf;
  M_BLR_Train<BF_Linear> m(bf, c);
  REQUIRE_THROWS_WITH(m.train(Phi, T), Catch::Contains("out of bounds"));
}

TEST_CASE("KRR linear partial fix places fixed values", "[fixed_weights]") {
  // Locks current KRR behaviour with Kern_Linear (no EKM projection):
  // fixed columns are zeroed for the solve and the fixed values written
  // back by index; free weights match the conditional LS reference.
  Context c = make_ctx();
  c.add("FIXINDEX", "4,1");
  c.add("FIXWEIGHT", "0.3");
  c.add("FIXWEIGHT", "-1.1");

  const phi_type Phi0 = make_phi();
  const t_type T0 = make_target(Phi0);

  Kern_Linear kern;
  M_KRR_Train<Kern_Linear> m(kern, c);
  phi_type Phi = Phi0;
  t_type T = T0;
  m.train(Phi, T);

  const t_type &w = m.get_weights();
  REQUIRE(w.size() == 4);
  REQUIRE(w[3] == 0.3);
  REQUIRE(w[0] == -1.1);

  const t_type wref = reference_free_weights(Phi0, T0, {4, 1}, {0.3, -1.1});
  REQUIRE(std::abs(w[1] - wref[0]) < 1e-8);
  REQUIRE(std::abs(w[2] - wref[1]) < 1e-8);
}
