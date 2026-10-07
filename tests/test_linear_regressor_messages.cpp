// test_linear_regressor_messages.cpp
//
// What LinearRegressor::train says when LAPACK cannot fit the weights:
// a rank-deficient design matrix under OALGO 2 (DGELS, QR) and OALGO 4
// (normal equations, Cholesky), and a design matrix or target with a value
// that is not finite, under every OALGO. The design matrices are laid out as
// DesignMatrix lays them out: the data rows, then one ridge row per column
// (zero with LAMBDA 0).

#include "catch2/catch.hpp"
#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/models/factor_handle.h>
#include <tadah/models/linear_regressor.h>
#include <tadah/models/memory/models_workspace_manager.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
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

Context make_ctx(const std::string &oalgo) {
  Context c;
  c.add("VERBOSE", "0");
  c.add("LAMBDA", "0");
  c.add("OALGO", oalgo);
  return c;
}

/* `data` rows of N columns, then N zero ridge rows; column j of data row i
 * is 1 + i + 2 j + i j % 3 for j < `nonzero`, 0 from there. */
phi_type design(size_t data, size_t n, size_t nonzero) {
  phi_type Phi(data + n, n);
  Phi.set_zero();
  for (size_t i = 0; i < data; ++i)
    for (size_t j = 0; j < nonzero; ++j)
      Phi(i, j) = 1.0 + static_cast<double>(i + 2 * j + (i * j) % 3);
  return Phi;
}

t_type target(size_t data, size_t n) {
  t_type T(data + n);
  for (size_t i = 0; i < T.size(); ++i)
    T[i] = i < data ? 0.5 * static_cast<double>(i + 1) : 0.0;
  return T;
}

/* The message LinearRegressor::train throws, or "no exception". */
std::string fit_refusal(const std::string &oalgo, phi_type Phi, t_type T) {
  Context c = make_ctx(oalgo);
  t_type w;
  Matrix Sigma;
  try {
    LinearRegressor::train(c, Phi, T, w, Sigma);
  } catch (const std::exception &e) {
    return e.what();
  }
  return "no exception";
}

const char *const kRemedies =
    " Use OALGO 1 (DGELSD) or OALGO 3 (SVD), which give the minimum-norm "
    "weights, LAMBDA > 0 (a ridge penalty makes every column it penalises "
    "full rank; it does not penalise the intercept columns of BIAS), or more "
    "data (structures, FORCE or STRESS rows).";

}  // namespace

/* Pins: with LAMBDA 0, a fit whose design matrix does not have full column
 * rank is refused by OALGO 2 and OALGO 4 with a message that says so, the
 * column where it shows, the data rows and weights, and the remedies:
 *  - OALGO 2, 2 data rows for 4 weights (DGELS: R(3,3) is exactly 0, the
 *    ridge rows being zero): "OALGO 2 (DGELS, QR) cannot fit these 4
 *    weights: the design matrix (2 data rows for 4 weights, LAMBDA 0) does
 *    not have full column rank: its column 3 is zero or a linear combination
 *    of the columns before it (fewer data rows than weights), so the
 *    least-squares weights are not unique." + the remedies;
 *  - OALGO 2 and OALGO 4, 6 data rows for 4 weights whose columns 3 and 4
 *    are zero: the same for column 3, and for OALGO 4 "the matrix Phi^T Phi
 *    of the design matrix (...) is not positive definite (its leading
 *    minor of order 3 is not)".
 * OALGO 1 and OALGO 3 fit both (the minimum-norm weights: the zero columns'
 * weights are 0 to 1e-12). (Pins as well.)
 * Before the fix: "Error in DGELS: info = 3" and "Error in DPOTRF:
 * tadah::core::Matrix is not positive definite; info = 3". */
TEST_CASE("LinearRegressor: OALGO 2 and 4 refuse a rank-deficient design "
          "matrix, saying why and what to do",
          "[linear_regressor][lapack_messages]") {
  const std::string few =
      "OALGO 2 (DGELS, QR) cannot fit these 4 weights: the design matrix (2 "
      "data rows for 4 weights, LAMBDA 0) does not have full column rank: its "
      "column 3 is zero or a linear combination of the columns before it "
      "(fewer data rows than weights), so the least-squares weights are not "
      "unique.";
  CHECK(fit_refusal("2", design(2, 4, 4), target(2, 4)) == few + kRemedies);

  const std::string zero_cols =
      "the design matrix (6 data rows for 4 weights, LAMBDA 0) does not have "
      "full column rank: its column 3 is zero or a linear combination of the "
      "columns before it, so the least-squares weights are not unique.";
  CHECK(fit_refusal("2", design(6, 4, 2), target(6, 4)) ==
        "OALGO 2 (DGELS, QR) cannot fit these 4 weights: " + zero_cols +
            kRemedies);
  CHECK(fit_refusal("4", design(6, 4, 2), target(6, 4)) ==
        std::string("OALGO 4 (normal equations, Cholesky) cannot fit these 4 "
                    "weights: the matrix Phi^T Phi of the design matrix (6 "
                    "data rows for 4 weights, LAMBDA 0) is not positive "
                    "definite (its leading minor of order 3 is not), so the "
                    "design matrix does not have full column rank or is too "
                    "ill-conditioned for the normal equations.") +
            kRemedies);

  for (const char *oalgo : {"1", "3"}) {
    INFO("OALGO " << oalgo);
    CHECK(fit_refusal(oalgo, design(2, 4, 4), target(2, 4)) == "no exception");
    Context c = make_ctx(oalgo);
    phi_type Phi = design(6, 4, 2);
    t_type T = target(6, 4);
    t_type w;
    Matrix Sigma;
    REQUIRE_NOTHROW(LinearRegressor::train(c, Phi, T, w, Sigma));
    REQUIRE(w.size() == 4);
    CHECK(std::abs(w[2]) < 1e-12);
    CHECK(std::abs(w[3]) < 1e-12);
  }
}

/* Pins: a design matrix or a target that holds a value that is not finite
 * (NaN, +-inf) is refused before LAPACK is called, under every OALGO,
 * naming the value's row and column (1-based) and that nothing was fitted:
 * "LinearRegressor::train: the design matrix holds a value that is not
 * finite (nan or inf) at row 2, column 3 (of 6 data rows, then 4 ridge
 * rows, and 4 columns): a descriptor or basis function gave it; no weights
 * were fitted", and for the target "the target vector holds a value that is
 * not finite (nan or inf) at row 5 (...)". The tests compare bits, not
 * std::isnan (the -ffast-math builds).
 * Before the fix: a NaN in the design matrix stopped in LAPACK with "Error
 * in DGELSD: info = -4" (OALGO 1), "Error in DGELS: info = -5" (OALGO 2),
 * "Error in DGESDD during SVD computation: info = 3" (OALGO 3) and "Error
 * in DPOTRF: tadah::core::Matrix is not positive definite; info = 3" (OALGO
 * 4); an infinity with "Error in DGELSD: info = 3", "Error in DGELS: info =
 * 1" and the same OALGO 3 and 4 messages; a NaN target with "Error in
 * DGELSD: info = -6" and "Error in DGELS: info = -7", and OALGO 3 and 4
 * fitted it without an error, as every OALGO fitted an infinite target. */
TEST_CASE("LinearRegressor: a design matrix or target that is not finite is "
          "refused under every OALGO, naming the row and column",
          "[linear_regressor][lapack_messages]") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  for (const char *oalgo : {"1", "2", "3", "4"}) {
    INFO("OALGO " << oalgo);
    for (const double bad : {nan, inf, -inf}) {
      INFO("value bits " << bits(bad));
      phi_type Phi = design(6, 4, 4);
      Phi(1, 2) = bad;
      CHECK(fit_refusal(oalgo, Phi, target(6, 4)) ==
            "LinearRegressor::train: the design matrix holds a value that is "
            "not finite (nan or inf) at row 2, column 3 (of 6 data rows, then "
            "4 ridge rows, and 4 columns): a descriptor or basis function "
            "gave it; no weights were fitted");
      t_type T = target(6, 4);
      T[4] = bad;
      CHECK(fit_refusal(oalgo, design(6, 4, 4), T) ==
            "LinearRegressor::train: the target vector holds a value that is "
            "not finite (nan or inf) at row 5 (of 6 data rows, then 4 ridge "
            "rows); no weights were fitted");
    }
  }
}

/* Pins: train_with_factor (OALGO 4 with the Cholesky factor cached for
 * tadah hpo's analytical gradients) refuses a non-finite design matrix as
 * train() does: "LinearRegressor::train: the design matrix holds a value
 * that is not finite (nan or inf) at row 2, column 3 (...)".
 * Before the fix: "OLS::solveAndCacheFactor: DPOTRF failed (matrix is not
 * positive definite); info = 3. This typically means a rank-deficient
 * design matrix or LAMBDA=0 with insufficient data; ...". */
TEST_CASE("LinearRegressor::train_with_factor refuses a design matrix that "
          "is not finite",
          "[linear_regressor][lapack_messages]") {
  Context c = make_ctx("4");
  phi_type Phi = design(6, 4, 4);
  Phi(1, 2) = std::numeric_limits<double>::quiet_NaN();
  t_type T = target(6, 4);
  t_type w;
  Matrix Sigma;
  FactorHandle factor;
  tadah::models::memory::ModelsWorkspaceManager wm;
  std::string msg = "no exception";
  try {
    LinearRegressor::train_with_factor(c, Phi, T, w, Sigma, factor, wm);
  } catch (const std::exception &e) {
    msg = e.what();
  }
  CHECK(msg ==
        "LinearRegressor::train: the design matrix holds a value that is not "
        "finite (nan or inf) at row 2, column 3 (of 6 data rows, then 4 ridge "
        "rows, and 4 columns): a descriptor or basis function gave it; no "
        "weights were fitted");
}
