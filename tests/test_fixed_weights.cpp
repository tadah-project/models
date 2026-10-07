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

#include <algorithm>
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

/* Pins: with FIXINDEX 2 / FIXWEIGHT 0.7, OALGO 3 (SVD) and LAMBDA 0.1 on
 * the 8 x 4 design matrix (with the 4 ridge-augmentation rows a design
 * matrix appends, zero for OALGO 3), M_BLR_Train<BF_Linear> and
 * M_KRR_Train<Kern_Linear> give the covariance of all 4 weights, as
 * predict --error multiplies it by the full design matrix: Sigma is 4 x 4,
 * its free block (rows and columns 1, 3, 4) is the covariance of the same
 * fit on the free columns alone (the 8 rows and 3 augmentation rows: the
 * same 8 observations; 1e-12 relative), row and column 2 of the pinned
 * weight are 0 (no variance), and the Context holds no SIGMA (the model
 * holds the covariance; the potential gets it as WCOV).
 * Before the fix: M_BLR_Train left the 3 x 3 covariance of the free
 * weights (SIGMA "3 ..." with 10 values); M_KRR_Train zeroed the pinned
 * column and kept it, so its Sigma was 4 x 4 but not that of the free
 * weights: its free block differed from the reference by up to 3.3, and
 * the pinned weight's row and column were not 0 (up to 6.7).
 * Before WCOV and WNOISE: the Context held SIGMA, 4 then the 16 values
 * (this case checked it), and the matrices of 8 rows had no augmentation
 * rows (their zero rows change no SIGMA). */
TEST_CASE("BLR and KRR with FIXINDEX and OALGO 3: the covariance of every "
          "weight, zero for a pinned one",
          "[fixed_weights][sigma]") {
  const phi_type Phi8 = make_phi();
  const t_type T8 = make_target(Phi8);
  // The design matrix of a fit: the 8 rows, then 4 augmentation rows.
  phi_type Phi0(12, 4);
  Phi0.set_zero();
  t_type T0(12);
  for (size_t i = 0; i < 12; ++i) T0[i] = i < 8 ? T8[i] : 0.0;
  for (size_t j = 0; j < 4; ++j)
    for (size_t i = 0; i < 8; ++i) Phi0(i, j) = Phi8(i, j);
  auto oalgo3 = [] {
    Context c;
    c.add("VERBOSE", "0");
    c.add("LAMBDA", "0.1");
    c.add("OALGO", "3");
    return c;
  };
  // The reference: the same fit on the free columns 0, 2, 3 (8 rows and 3
  // augmentation rows).
  phi_type Phi_r(11, 3);
  Phi_r.set_zero();
  t_type T_r(11);
  const size_t free_cols[3] = {0, 2, 3};
  for (size_t j = 0; j < 3; ++j)
    for (size_t i = 0; i < 8; ++i) Phi_r(i, j) = Phi8(i, free_cols[j]);
  for (size_t i = 0; i < 11; ++i)
    T_r[i] = i < 8 ? T8[i] - Phi8(i, 1) * 0.7 : 0.0;
  Context cr = oalgo3();
  BF_Linear bfr;
  M_BLR_Train<BF_Linear> ref(bfr, cr);
  ref.train(Phi_r, T_r);
  REQUIRE(ref.Sigma.rows() == 3);
  double smax = 0;
  for (size_t a = 0; a < 3; ++a)
    for (size_t b = 0; b < 3; ++b) smax = std::max(smax, std::abs(ref.Sigma(a, b)));

  for (int krr = 0; krr < 2; ++krr) DYNAMIC_SECTION(
      (krr ? "M_KRR_Train<Kern_Linear>" : "M_BLR_Train<BF_Linear>")) {
    Context c = oalgo3();
    c.add("FIXINDEX", "2");
    c.add("FIXWEIGHT", "0.7");
    phi_type Phi = Phi0;
    t_type T = T0;
    Matrix S;
    if (krr) {
      Kern_Linear k;
      M_KRR_Train<Kern_Linear> m(k, c);
      m.train(Phi, T);
      CHECK(m.get_weights()[1] == 0.7);
      S = m.Sigma;
    } else {
      BF_Linear bf;
      M_BLR_Train<BF_Linear> m(bf, c);
      m.train(Phi, T);
      CHECK(m.get_weights()[1] == 0.7);
      S = m.Sigma;
    }
    REQUIRE(S.rows() == 4);
    REQUIRE(S.cols() == 4);
    for (size_t a = 0; a < 3; ++a)
      for (size_t b = 0; b < 3; ++b)
        CHECK(std::abs(S(free_cols[a], free_cols[b]) - ref.Sigma(a, b)) <
              1e-12 * smax);
    for (size_t k = 0; k < 4; ++k) {
      CHECK(S(1, k) == 0.0);
      CHECK(S(k, 1) == 0.0);
    }
    CHECK_FALSE(c.has_stored("SIGMA"));
  }
}
