// test_weight_covariance.cpp
//
// The weight covariance an OALGO 3 fit (LinearRegressor::train, the SVD
// solve) computes: the covariance of the weights it returns, in their own
// parameterisation,
//
//     Sigma = s^2 sum_{kept i} v_i v_i^T / (s_i^2 + lambda),
//     s^2   = RSS / (n_obs - gamma),  gamma = sum_{kept i} s_i^2/(s_i^2 + lambda),
//
// where A = U S V^T is the (weighted) design matrix, the kept directions
// are those the solve keeps (s_i > rcond * s_max), RSS is the residual sum
// of squares of the fitted weights and n_obs the number of observations:
// the rows of A without the ridge-augmentation rows the design matrix
// appends (one per column, all zero for OALGO 3). With every direction kept
// Sigma = s^2 (A^T A + lambda I)^-1, the posterior covariance of Bayesian
// linear regression with the noise variance estimated from the residuals.
//
// The references are computed here in long double from the real rows of A
// alone (Gauss-Jordan inverses, Gram-Schmidt bases), never through the
// library's SVD.

#include "catch2/catch.hpp"

#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/core/logger.h>
#include <tadah/models/functions/basis_functions/bf_linear.h>
#include <tadah/models/functions/kernels/kern_linear.h>
#include <tadah/models/functions/kernels/kern_rbf.h>
#include <tadah/models/linear_regressor.h>
#include <tadah/models/m_blr_train.h>
#include <tadah/models/m_krr_train.h>
#include <tadah/models/memory/models_workspace_manager.h>
#include <tadah/models/ols.h>
#include <tadah/models/svd.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace tadah::core;
using namespace tadah::models;

namespace {

using LD = long double;
using LMat = std::vector<std::vector<LD>>;  // row-major, [row][col]

/* True unless the bits of @p v are those of an infinity or a NaN (the
 * build may use -ffast-math, so std::isfinite is not relied on). */
bool finite_bits(double v) {
  std::uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return ((u >> 52) & 0x7ffu) != 0x7ffu;
}

/* A deterministic, well-conditioned m x n matrix (row-major, long
 * double): entries in [-0.5, 0.5) plus 2 on the diagonal. */
LMat make_A(size_t m, size_t n, double seed) {
  LMat A(m, std::vector<LD>(n));
  double v = seed;
  for (size_t j = 0; j < n; ++j)
    for (size_t i = 0; i < m; ++i) {
      v = std::fmod(v * 1.618 + 0.31 * static_cast<double>(j + 1), 1.0) - 0.5;
      A[i][j] = v + (i == j ? 2.0 : 0.0);
    }
  return A;
}

/* Targets of A: A w_true plus a perturbation, so the residuals (and s^2)
 * are not zero. */
std::vector<LD> make_T(const LMat &A) {
  const size_t m = A.size(), n = A[0].size();
  std::vector<LD> T(m, 0.0L);
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = 0; j < n; ++j)
      T[i] += A[i][j] * 0.5L * static_cast<LD>(j + 1) * (j % 2 ? -1 : 1);
    T[i] += 0.05L * std::sin(1.7L * static_cast<LD>(i) + 0.3L);
  }
  return T;
}

/* The design matrix as the DesignMatrix gives it: the m real rows of A,
 * then n ridge-augmentation rows, zero (OALGO 3 writes none). */
phi_type with_augmentation(const LMat &A) {
  const size_t m = A.size(), n = A[0].size();
  phi_type Phi(m + n, n);
  Phi.set_zero();
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j) Phi(i, j) = static_cast<double>(A[i][j]);
  return Phi;
}

t_type target_with_augmentation(const std::vector<LD> &T, size_t n) {
  t_type t(T.size() + n);
  for (size_t i = 0; i < t.size(); ++i)
    t[i] = i < T.size() ? static_cast<double>(T[i]) : 0.0;
  return t;
}

LMat transpose_times(const LMat &A, const LMat &B) {  // A^T B
  const size_t m = A.size(), n = A[0].size(), p = B[0].size();
  LMat C(n, std::vector<LD>(p, 0.0L));
  for (size_t i = 0; i < m; ++i)
    for (size_t a = 0; a < n; ++a)
      for (size_t b = 0; b < p; ++b) C[a][b] += A[i][a] * B[i][b];
  return C;
}

/* Gauss-Jordan inverse with partial pivoting (long double). */
LMat inverse(LMat M) {
  const size_t n = M.size();
  LMat I(n, std::vector<LD>(n, 0.0L));
  for (size_t i = 0; i < n; ++i) I[i][i] = 1.0L;
  for (size_t c = 0; c < n; ++c) {
    size_t p = c;
    for (size_t r = c + 1; r < n; ++r)
      if (std::fabs(M[r][c]) > std::fabs(M[p][c])) p = r;
    std::swap(M[c], M[p]);
    std::swap(I[c], I[p]);
    const LD d = M[c][c];
    for (size_t k = 0; k < n; ++k) {
      M[c][k] /= d;
      I[c][k] /= d;
    }
    for (size_t r = 0; r < n; ++r)
      if (r != c) {
        const LD f = M[r][c];
        for (size_t k = 0; k < n; ++k) {
          M[r][k] -= f * M[c][k];
          I[r][k] -= f * I[c][k];
        }
      }
  }
  return I;
}

/* The ridge fit of the real rows: w = (A^T A + lambda I)^-1 A^T T, its RSS,
 * gamma = tr(A^T A (A^T A + lambda I)^-1), s^2 = RSS/(m - gamma) and
 * Sigma = s^2 (A^T A + lambda I)^-1. */
struct RidgeRef {
  std::vector<LD> w;
  LMat Sigma;
  LD rss = 0, gamma = 0, s2 = 0;
};

RidgeRef ridge_reference(const LMat &A, const std::vector<LD> &T, LD lambda) {
  const size_t m = A.size(), n = A[0].size();
  LMat AtA = transpose_times(A, A);
  LMat M = AtA;
  for (size_t i = 0; i < n; ++i) M[i][i] += lambda;
  const LMat Minv = inverse(M);
  std::vector<LD> AtT(n, 0.0L);
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j) AtT[j] += A[i][j] * T[i];
  RidgeRef r;
  r.w.assign(n, 0.0L);
  for (size_t a = 0; a < n; ++a)
    for (size_t b = 0; b < n; ++b) r.w[a] += Minv[a][b] * AtT[b];
  for (size_t i = 0; i < m; ++i) {
    LD e = T[i];
    for (size_t j = 0; j < n; ++j) e -= A[i][j] * r.w[j];
    r.rss += e * e;
  }
  for (size_t a = 0; a < n; ++a)
    for (size_t b = 0; b < n; ++b) r.gamma += AtA[a][b] * Minv[b][a];
  r.s2 = r.rss / (static_cast<LD>(m) - r.gamma);
  r.Sigma = Minv;
  for (auto &row : r.Sigma)
    for (auto &x : row) x *= r.s2;
  return r;
}

/* max |S - R| / max |R|: S the library's covariance, R the reference. */
double rel_diff(const Matrix &S, const LMat &R) {
  LD dmax = 0, rmax = 0;
  for (size_t a = 0; a < R.size(); ++a)
    for (size_t b = 0; b < R.size(); ++b) {
      dmax = std::max(dmax, std::fabs(static_cast<LD>(S(a, b)) - R[a][b]));
      rmax = std::max(rmax, std::fabs(R[a][b]));
    }
  return static_cast<double>(dmax / rmax);
}

Context oalgo3_context(const std::string &lambda) {
  Context c;
  c.add("VERBOSE", "0");
  c.add("LAMBDA", lambda);
  c.add("OALGO", "3");
  return c;
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

/* Orthonormal columns from the columns of @p B (Gram-Schmidt twice). */
LMat orthonormal_columns(LMat B) {
  const size_t m = B.size(), n = B[0].size();
  for (size_t j = 0; j < n; ++j)
    for (int pass = 0; pass < 2; ++pass) {
      for (size_t k = 0; k < j; ++k) {
        LD d = 0;
        for (size_t i = 0; i < m; ++i) d += B[i][k] * B[i][j];
        for (size_t i = 0; i < m; ++i) B[i][j] -= d * B[i][k];
      }
      LD nrm = 0;
      for (size_t i = 0; i < m; ++i) nrm += B[i][j] * B[i][j];
      nrm = std::sqrt(nrm);
      for (size_t i = 0; i < m; ++i) B[i][j] /= nrm;
    }
  return B;
}

}  // namespace

/* Pins: with OALGO 3 the covariance LinearRegressor::train returns (its
 * Sigma out-parameter) is s^2 (A^T A + lambda I)^-1, s^2 = RSS/(m - gamma)
 * with m = 9 the real rows of a random 9 x 4 A (the 4 ridge-augmentation
 * rows of the design matrix are not observations), for LAMBDA 1e-3 and
 * LAMBDA 0, to 1e-10 relative; the weights are the ridge solution (1e-10).
 * M_BLR_Train<BF_Linear>, M_KRR_Train<Kern_Linear> and
 * M_KRR_Train<Kern_RBF> (whose design-matrix columns are its basis
 * vectors, kept as they are) hold the same covariance in their Sigma.
 * Before the fix: get_sigma summed VT[j+i*n]*VT[k+i*n]/(s_i^2+lambda),
 * V^T R V instead of V R V^T, unscaled (no s^2 = 1.53e-3): Sigma(0,0) was
 * 0.1338 against 2.087e-4 at LAMBDA 1e-3 and at LAMBDA 0, for all four. */
TEST_CASE("OALGO 3: the weight covariance is s^2 (A^T A + lambda I)^-1",
          "[weight_covariance]") {
  const LMat A = make_A(9, 4, 0.17);
  const std::vector<LD> T = make_T(A);
  for (const char *lambda : {"1e-3", "0"}) {
    const RidgeRef ref = ridge_reference(A, T, std::stold(lambda));
    for (const char *who : {"LinearRegressor", "M_BLR_Train<BF_Linear>",
                            "M_KRR_Train<Kern_Linear>",
                            "M_KRR_Train<Kern_RBF>"}) {
      DYNAMIC_SECTION(who << ", LAMBDA " << lambda) {
        Context c = oalgo3_context(lambda);
        c.add("MPARAMS", "1.0");  // Kern_RBF's sigma
        phi_type Phi = with_augmentation(A);
        t_type t = target_with_augmentation(T, 4);
        t_type w;
        Matrix S;
        const std::string name = who;
        if (name == "LinearRegressor") {
          LinearRegressor::train(c, Phi, t, w, S);
        } else if (name == "M_BLR_Train<BF_Linear>") {
          BF_Linear bf;
          M_BLR_Train<BF_Linear> m(bf, c);
          m.train(Phi, t);
          w = m.get_weights();
          S = m.Sigma;
        } else if (name == "M_KRR_Train<Kern_Linear>") {
          Kern_Linear k;
          M_KRR_Train<Kern_Linear> m(k, c);
          m.train(Phi, t);
          w = m.get_weights();
          S = m.Sigma;
        } else {
          Kern_RBF k(c);
          M_KRR_Train<Kern_RBF> m(k, c);
          m.train(Phi, t);
          w = m.get_weights();
          S = m.Sigma;
        }
        REQUIRE(w.size() == 4);
        for (size_t j = 0; j < 4; ++j)
          CHECK(std::fabs(w[j] - static_cast<double>(ref.w[j])) <
                1e-10 * std::fabs(static_cast<double>(ref.w[j])));
        REQUIRE(S.rows() == 4);
        REQUIRE(S.cols() == 4);
        INFO("s^2 " << static_cast<double>(ref.s2) << ", Sigma(0,0) "
                    << S(0, 0) << " against "
                    << static_cast<double>(ref.Sigma[0][0]));
        CHECK(rel_diff(S, ref.Sigma) < 1e-10);
      }
    }
  }
}

/* Pins: with OALGO 3, LAMBDA 0 and two equal columns (A = [a0 a1 a2 a2],
 * 9 real rows) the covariance is that of the minimum-norm weights the
 * solve returns, truncated as the solve is: finite, Sigma (0,0,1,-1) = 0
 * (no variance along the undetermined direction, 1e-10 relative), equal
 * to s^2 E+ (B^T B)^-1 E+^T with B = [a0 a1 a2], A = B E, E+ the
 * pseudo-inverse of E, and s^2 = RSS/(9 - 3) (gamma = 3, the kept
 * directions); a WARNING says that 1 of the 4 singular directions is
 * dropped from the weights and from WCOV.
 * Before the fix: every singular value got 1/(s_i^2 + lambda), so the
 * transposed covariance had Sigma(0,0) = 3.8e30 (s_4 of order 1e-16)
 * against 2.2e-4, columns 3 and 4 differed, and no WARNING was logged. */
TEST_CASE("OALGO 3, LAMBDA 0, rank deficient: the covariance of the "
          "minimum-norm weights",
          "[weight_covariance]") {
  LMat B = make_A(9, 3, 0.29);
  LMat A(9, std::vector<LD>(4));
  for (size_t i = 0; i < 9; ++i) {
    for (size_t j = 0; j < 3; ++j) A[i][j] = B[i][j];
    A[i][3] = B[i][2];
  }
  const std::vector<LD> T = make_T(A);
  // The fit of B: its weights (the solve's minimum-norm weights of A split
  // w_2 evenly over the two equal columns) and residuals.
  const RidgeRef rb = ridge_reference(B, T, 0.0L);
  REQUIRE(static_cast<double>(rb.gamma) == Approx(3.0).epsilon(1e-12));
  // E+ = E^T (E E^T)^-1 with E = [I_3 | e_3]: E E^T = I + e_3 e_3^T.
  LMat E(3, std::vector<LD>(4, 0.0L));
  for (size_t i = 0; i < 3; ++i) E[i][i] = 1.0L;
  E[2][3] = 1.0L;
  LMat EEt(3, std::vector<LD>(3, 0.0L));
  for (size_t a = 0; a < 3; ++a)
    for (size_t b = 0; b < 3; ++b)
      for (size_t k = 0; k < 4; ++k) EEt[a][b] += E[a][k] * E[b][k];
  const LMat EEt_inv = inverse(EEt);
  LMat Ep(4, std::vector<LD>(3, 0.0L));  // E^T (E E^T)^-1
  for (size_t k = 0; k < 4; ++k)
    for (size_t b = 0; b < 3; ++b)
      for (size_t a = 0; a < 3; ++a) Ep[k][b] += E[a][k] * EEt_inv[a][b];
  // B^T B ^-1 = Sigma_B / s^2_B; the reference uses A's s^2 = RSS/(9-3),
  // the same residuals as B's fit.
  LMat ref(4, std::vector<LD>(4, 0.0L));
  for (size_t p = 0; p < 4; ++p)
    for (size_t q = 0; q < 4; ++q)
      for (size_t a = 0; a < 3; ++a)
        for (size_t b = 0; b < 3; ++b)
          ref[p][q] += Ep[p][a] * rb.Sigma[a][b] * Ep[q][b];

  Context c = oalgo3_context("0");
  phi_type Phi = with_augmentation(A);
  t_type t = target_with_augmentation(T, 4);
  t_type w;
  Matrix S;
  std::string log;
  {
    LogCapture cap;
    LinearRegressor::train(c, Phi, t, w, S);
    log = cap.text();
  }
  INFO("log: " << log);
  REQUIRE(S.rows() == 4);
  REQUIRE(S.cols() == 4);
  bool finite = true;
  for (size_t a = 0; a < 4; ++a)
    for (size_t b = 0; b < 4; ++b) finite = finite && finite_bits(S(a, b));
  REQUIRE(finite);
  LD smax = 0;
  for (size_t a = 0; a < 4; ++a)
    for (size_t b = 0; b < 4; ++b) smax = std::max(smax, std::fabs(ref[a][b]));
  for (size_t a = 0; a < 4; ++a)
    CHECK(std::fabs(S(a, 2) - S(a, 3)) < 1e-10 * static_cast<double>(smax));
  INFO("Sigma(0,0) " << S(0, 0) << " against "
                     << static_cast<double>(ref[0][0]));
  CHECK(rel_diff(S, ref) < 1e-10);
  CHECK(log.find("1 of 4 singular directions") != std::string::npos);
  CHECK(log.find("WCOV") != std::string::npos);
}

/* Pins: with OALGO 3 and LAMBDA 1e-6 > 0, a direction the solve drops
 * (singular values 10, 3, 1 and 1e-10 of a 9 x 4 A = U S V^T, so s_4 is
 * below rcond 1e-8 x s_max) is dropped from the covariance too (owner
 * decision D4: one truncation rule for every lambda, the solve's):
 * Sigma = s^2 sum_{i<3} v_i v_i^T / (s_i^2 + lambda), s^2 = RSS/(9 -
 * gamma), gamma = sum_{i<3} s_i^2/(s_i^2 + lambda), to 1e-9 relative; the
 * weights have no component along v_4 (1e-12 of their norm); a WARNING
 * names 1 of the 4 directions.
 * Before the fix: the dropped direction kept 1/(s_4^2 + lambda) = 1e6 in
 * the transposed sum: Sigma(0,0) was 2.75e5 against 2.13e-4, and no
 * WARNING was logged. */
TEST_CASE("OALGO 3 with LAMBDA > 0: a direction the solve drops has no "
          "covariance",
          "[weight_covariance]") {
  const LMat U = orthonormal_columns(make_A(9, 4, 0.41));
  const LMat V = orthonormal_columns(make_A(4, 4, 0.53));
  const LD s[4] = {10.0L, 3.0L, 1.0L, 1e-10L};
  const LD lambda = 1e-6L;
  LMat A(9, std::vector<LD>(4, 0.0L));
  for (size_t i = 0; i < 9; ++i)
    for (size_t j = 0; j < 4; ++j)
      for (size_t k = 0; k < 4; ++k) A[i][j] += U[i][k] * s[k] * V[j][k];
  const std::vector<LD> T = make_T(A);
  // The solve's weights: the kept directions only.
  std::vector<LD> w_ref(4, 0.0L);
  LD gamma = 0;
  for (size_t k = 0; k < 3; ++k) {
    LD utb = 0;
    for (size_t i = 0; i < 9; ++i) utb += U[i][k] * T[i];
    for (size_t j = 0; j < 4; ++j)
      w_ref[j] += V[j][k] * s[k] / (s[k] * s[k] + lambda) * utb;
    gamma += s[k] * s[k] / (s[k] * s[k] + lambda);
  }
  LD rss = 0;
  for (size_t i = 0; i < 9; ++i) {
    LD e = T[i];
    for (size_t j = 0; j < 4; ++j) e -= A[i][j] * w_ref[j];
    rss += e * e;
  }
  const LD s2 = rss / (9.0L - gamma);
  LMat ref(4, std::vector<LD>(4, 0.0L));
  for (size_t k = 0; k < 3; ++k)
    for (size_t a = 0; a < 4; ++a)
      for (size_t b = 0; b < 4; ++b)
        ref[a][b] += s2 * V[a][k] * V[b][k] / (s[k] * s[k] + lambda);

  Context c = oalgo3_context("1e-6");
  phi_type Phi = with_augmentation(A);
  t_type t = target_with_augmentation(T, 4);
  t_type w;
  Matrix S;
  std::string log;
  {
    LogCapture cap;
    LinearRegressor::train(c, Phi, t, w, S);
    log = cap.text();
  }
  INFO("log: " << log);
  REQUIRE(w.size() == 4);
  LD wn = 0, w4 = 0;
  for (size_t j = 0; j < 4; ++j) {
    wn += static_cast<LD>(w[j]) * w[j];
    w4 += static_cast<LD>(w[j]) * V[j][3];
  }
  CHECK(static_cast<double>(std::fabs(w4)) <
        1e-12 * static_cast<double>(std::sqrt(wn)));
  REQUIRE(S.rows() == 4);
  INFO("Sigma(0,0) " << S(0, 0) << " against "
                     << static_cast<double>(ref[0][0]));
  CHECK(rel_diff(S, ref) < 1e-9);
  CHECK(log.find("1 of 4 singular directions") != std::string::npos);
}

/* Pins: LAMBDA -1 (the evidence approximation, OALGO 3) counts the
 * observations only: on the 9 real rows of a 9 x 4 A (4
 * ridge-augmentation rows appended) it converges to the lambda = alpha /
 * beta of the iteration alpha = gamma / |w|^2, beta = (N - gamma) / RSS
 * with N = 9 (gamma over the kept directions; long double, iterated to
 * 1e-15) to 1e-8 relative, and to its weights (1e-8); the covariance is
 * s^2 (A^T A + lambda I)^-1 of that lambda with s^2 = RSS/(9 - gamma)
 * (1e-7 relative).
 * Before the fix: N was max(rows, cols) = 13, which counted the
 * augmentation rows: beta was too large by about (13 - gamma)/(9 - gamma),
 * so lambda was 8.13e-3 against 1.469e-2, the weights followed, and the
 * covariance was V^T R V. */
TEST_CASE("LAMBDA -1: the evidence approximation counts the observations",
          "[weight_covariance][ea]") {
  const LMat A = make_A(9, 4, 0.71);
  std::vector<LD> T = make_T(A);
  for (size_t i = 0; i < T.size(); ++i)
    T[i] += 0.2L * std::cos(2.3L * static_cast<LD>(i));
  // Reference EA: alpha, beta from 1, 1.
  LD alpha = 1.0L, beta = 1.0L;
  RidgeRef r;
  for (int it = 0; it < 10000; ++it) {
    r = ridge_reference(A, T, alpha / beta);
    LD ww = 0;
    for (LD x : r.w) ww += x * x;
    const LD a_new = r.gamma / ww;
    const LD b_new = (9.0L - r.gamma) / r.rss;
    const bool done = std::fabs(a_new - alpha) < 1e-15L * a_new &&
                      std::fabs(b_new - beta) < 1e-15L * b_new;
    alpha = a_new;
    beta = b_new;
    if (done) break;
  }
  const LD lambda_ref = alpha / beta;
  r = ridge_reference(A, T, lambda_ref);

  Context c = oalgo3_context("-1");
  c.add("ALPHA", "1");
  c.add("BETA", "1");
  phi_type Phi = with_augmentation(A);
  t_type t = target_with_augmentation(T, 4);
  t_type w;
  Matrix S;
  LinearRegressor::train(c, Phi, t, w, S);
  const double lambda = c.get<double>("ALPHA") / c.get<double>("BETA");
  INFO("lambda " << lambda << " against " << static_cast<double>(lambda_ref));
  CHECK(std::fabs(lambda - static_cast<double>(lambda_ref)) <
        1e-8 * static_cast<double>(lambda_ref));
  REQUIRE(w.size() == 4);
  for (size_t j = 0; j < 4; ++j)
    CHECK(std::fabs(w[j] - static_cast<double>(r.w[j])) <
          1e-8 * std::fabs(static_cast<double>(r.w[j])));
  REQUIRE(S.rows() == 4);
  CHECK(rel_diff(S, r.Sigma) < 1e-7);
}

/* Pins: with FIXINDEX 2 / FIXWEIGHT 0.7 and OALGO 3, M_BLR_Train<BF_Linear>
 * and M_KRR_Train<Kern_Linear> hold the covariance of every weight (4 x
 * 4): the pinned weight's row and column are 0 and the free block (rows
 * and columns 1, 3, 4) is s^2 (A_f^T A_f + lambda I)^-1 of the fit of the
 * free columns A_f to T - 0.7 a_2, with s^2 = RSS/(9 - gamma) over the 9
 * real rows (Phi keeps all 13 rows when the pinned column leaves it, so
 * rows - free columns = 10 would count the pinned weight as an
 * observation), for LAMBDA 1e-3 and 0 (1e-10 relative); at LAMBDA 0 no
 * WARNING about dropped directions (the pinned column is removed, not
 * zeroed and kept).
 * Before the fix: the free block was V^T R V of the free columns,
 * unscaled: Sigma(0,0) 0.155 against 0.406. */
TEST_CASE("FIXINDEX with OALGO 3: one covariance row per weight, the free "
          "block of the free fit",
          "[weight_covariance][fixed_weights]") {
  const LMat A = make_A(9, 4, 0.23);
  const std::vector<LD> T = make_T(A);
  const size_t free_cols[3] = {0, 2, 3};
  LMat Af(9, std::vector<LD>(3));
  std::vector<LD> Tr = T;
  for (size_t i = 0; i < 9; ++i) {
    for (size_t j = 0; j < 3; ++j) Af[i][j] = A[i][free_cols[j]];
    Tr[i] -= 0.7L * A[i][1];
  }
  for (const char *lambda : {"1e-3", "0"}) {
    const RidgeRef ref = ridge_reference(Af, Tr, std::stold(lambda));
    for (int krr = 0; krr < 2; ++krr) {
      DYNAMIC_SECTION((krr ? "M_KRR_Train<Kern_Linear>" : "M_BLR_Train<BF_Linear>")
                      << ", LAMBDA " << lambda) {
        Context c = oalgo3_context(lambda);
        c.add("FIXINDEX", "2");
        c.add("FIXWEIGHT", "0.7");
        phi_type Phi = with_augmentation(A);
        t_type t = target_with_augmentation(T, 4);
        Matrix S;
        std::string log;
        {
          LogCapture cap;
          if (krr) {
            Kern_Linear k;
            M_KRR_Train<Kern_Linear> m(k, c);
            m.train(Phi, t);
            CHECK(m.get_weights()[1] == 0.7);
            S = m.Sigma;
          } else {
            BF_Linear bf;
            M_BLR_Train<BF_Linear> m(bf, c);
            m.train(Phi, t);
            CHECK(m.get_weights()[1] == 0.7);
            S = m.Sigma;
          }
          log = cap.text();
        }
        INFO("log: " << log);
        REQUIRE(S.rows() == 4);
        REQUIRE(S.cols() == 4);
        for (size_t k = 0; k < 4; ++k) {
          CHECK(S(1, k) == 0.0);
          CHECK(S(k, 1) == 0.0);
        }
        LD dmax = 0, rmax = 0;
        for (size_t a = 0; a < 3; ++a)
          for (size_t b = 0; b < 3; ++b) {
            dmax = std::max(dmax, std::fabs(static_cast<LD>(S(free_cols[a], free_cols[b])) -
                                            ref.Sigma[a][b]));
            rmax = std::max(rmax, std::fabs(ref.Sigma[a][b]));
          }
        INFO("Sigma(0,0) " << S(0, 0) << " against "
                           << static_cast<double>(ref.Sigma[0][0]));
        CHECK(static_cast<double>(dmax / rmax) < 1e-10);
        CHECK(log.find("singular directions") == std::string::npos);
      }
    }
  }
}

/* Pins: with OALGO 3 and FIXINDEX 1-4 (every weight pinned, nothing
 * solved) M_BLR_Train and M_KRR_Train hold a covariance of 4 x 4 zeros
 * (pinned weights have no variance) and get_weights_uncertainty() returns
 * four zeros.
 * Before the fix: Sigma was empty (0 x 0; M_KRR_Train's zeroed pinned
 * columns gave 1/lambda, inf at LAMBDA 0, before fix/multispecies-bias)
 * and get_weights_uncertainty() threw at LAMBDA 0 ("Sigma matrix is not
 * computed."). */
TEST_CASE("OALGO 3 with every weight pinned: a zero covariance",
          "[weight_covariance][fixed_weights]") {
  const LMat A = make_A(9, 4, 0.37);
  const std::vector<LD> T = make_T(A);
  for (int krr = 0; krr < 2; ++krr) {
    DYNAMIC_SECTION((krr ? "M_KRR_Train<Kern_Linear>" : "M_BLR_Train<BF_Linear>")) {
      Context c = oalgo3_context("0");
      c.add("FIXINDEX", "1-4");
      c.add("FIXWEIGHT", "0.1 0.2 0.3 0.4");
      phi_type Phi = with_augmentation(A);
      t_type t = target_with_augmentation(T, 4);
      Matrix S;
      t_type unc;
      std::string msg;
      auto run = [&](auto &m) {
        m.train(Phi, t);
        S = m.Sigma;
        try {
          unc = m.get_weights_uncertainty();
        } catch (const std::exception &e) {
          msg = e.what();
        }
      };
      if (krr) {
        Kern_Linear k;
        M_KRR_Train<Kern_Linear> m(k, c);
        run(m);
      } else {
        BF_Linear bf;
        M_BLR_Train<BF_Linear> m(bf, c);
        run(m);
      }
      INFO("get_weights_uncertainty: " << msg);
      REQUIRE(S.rows() == 4);
      REQUIRE(S.cols() == 4);
      for (size_t a = 0; a < 4; ++a)
        for (size_t b = 0; b < 4; ++b) CHECK(S(a, b) == 0.0);
      REQUIRE(unc.size() == 4);
      for (size_t a = 0; a < 4; ++a) CHECK(unc[a] == 0.0);
    }
  }
}

/* Pins: get_weights_uncertainty() returns the standard deviations of the
 * weights, sqrt(diag Sigma), when the last fit computed a covariance
 * (OALGO 3, LAMBDA 0 included: 4 values, 1e-12 relative to the square
 * roots of the reference diagonal), and throws otherwise, naming OALGO 3:
 * after an OALGO 1 fit, also when the model was fitted with OALGO 3 before
 * (reset_trained(), OALGO 1 in its Context, trained again: Sigma is empty),
 * for M_BLR_Train<BF_Linear> and M_KRR_Train<Kern_Linear>.
 * Before the fix: it tested LAMBDA, not the fit: at LAMBDA 0 it threw
 * "Sigma matrix is not computed." (M_KRR: "only computed for LAMBDA !=
 * 0") also after an OALGO 3 fit, at LAMBDA 1e-3 an OALGO 1 fit returned the
 * diagonal of the Sigma of the earlier OALGO 3 fit (a model kept its Sigma
 * when a fit computed none), and it returned variances, not standard
 * deviations. */
TEST_CASE("get_weights_uncertainty: the standard deviations of an OALGO 3 "
          "fit, refused after any other fit",
          "[weight_covariance][uncertainty]") {
  const LMat A = make_A(9, 4, 0.61);
  const std::vector<LD> T = make_T(A);
  const RidgeRef ref = ridge_reference(A, T, 0.0L);
  for (int krr = 0; krr < 2; ++krr) {
    DYNAMIC_SECTION((krr ? "M_KRR_Train<Kern_Linear>" : "M_BLR_Train<BF_Linear>")) {
      auto check_model = [&](auto &m, Context &c) {
        phi_type Phi = with_augmentation(A);
        t_type t = target_with_augmentation(T, 4);
        m.train(Phi, t);
        std::string msg;
        t_type unc;
        try {
          unc = m.get_weights_uncertainty();
        } catch (const std::exception &e) {
          msg = e.what();
        }
        INFO("OALGO 3: " << msg);
        REQUIRE(unc.size() == 4);
        for (size_t j = 0; j < 4; ++j)
          CHECK(std::fabs(unc[j] - std::sqrt(static_cast<double>(ref.Sigma[j][j]))) <
                1e-12 * std::sqrt(static_cast<double>(ref.Sigma[j][j])));

        // The same model fitted again with OALGO 1 (and LAMBDA 1e-3).
        m.reset_trained();
        c.remove("OALGO");
        c.add("OALGO", "1");
        c.remove("LAMBDA");
        c.add("LAMBDA", "1e-3");
        phi_type Phi2 = with_augmentation(A);
        t_type t2 = target_with_augmentation(T, 4);
        m.train(Phi2, t2);
        CHECK(m.Sigma.rows() == 0);
        msg.clear();
        unc = t_type();
        try {
          unc = m.get_weights_uncertainty();
        } catch (const std::exception &e) {
          msg = e.what();
        }
        INFO("OALGO 1: " << msg << " (" << unc.size() << " values)");
        CHECK(unc.size() == 0);
        CHECK(msg.find("OALGO 3") != std::string::npos);
      };
      Context c = oalgo3_context("0");
      if (krr) {
        Kern_Linear k;
        M_KRR_Train<Kern_Linear> m(k, c);
        check_model(m, c);
      } else {
        BF_Linear bf;
        M_BLR_Train<BF_Linear> m(bf, c);
        check_model(m, c);
      }
    }
  }
}

/* Pins: LAMBDA -1 (the evidence approximation, OALGO 3) refuses a fit
 * with no more observations than effective weights, before it estimates
 * the noise: a design matrix of 3 rows of weight 0 (zero rows) and its 2
 * augmentation rows, n_obs 0, throws std::runtime_error naming LAMBDA -1,
 * the 0 observations and the remedy (LAMBDA >= 0), with no covariance (Sigma
 * 0 x 0) and the Context's ALPHA and BETA as they were (1 and 1); through
 * M_BLR_Train<BF_Linear>::train(Phi, T, 0) the model is not trained and
 * has no covariance; and through train(Phi, T) on the augmentation rows
 * alone (rows - cols = 0 observations) the same.
 * Before the fix: passes as well (a pin of the refusal MODELS 9b2629f
 * added with n_obs; before that commit EA took N = max(rows, cols), which
 * counted these rows, and divided by RSS = 0: beta was infinite and there
 * was no refusal). */
TEST_CASE("LAMBDA -1 refuses a fit without observations",
          "[weight_covariance][ea]") {
  auto zero_phi = [](size_t rows) {
    phi_type Phi(rows, 2);
    Phi.set_zero();
    return Phi;
  };
  auto zero_t = [](size_t rows) {
    t_type t(rows);
    for (size_t i = 0; i < rows; ++i) t[i] = 0.0;
    return t;
  };
  {
    Context c = oalgo3_context("-1");
    c.add("ALPHA", "1");
    c.add("BETA", "1");
    phi_type Phi = zero_phi(5);
    t_type t = zero_t(5);
    t_type w;
    Matrix S;
    tadah::models::memory::ModelsWorkspaceManager wsm;
    std::string msg = "<no exception>";
    try {
      LinearRegressor::train(c, Phi, t, w, S, wsm, {}, 0,
                             LinearRegressor::Layout::of(Phi));
    } catch (const std::runtime_error &e) {
      msg = e.what();
    }
    INFO(msg);
    CHECK(msg.find("LAMBDA -1") != std::string::npos);
    CHECK(msg.find("0 observations") != std::string::npos);
    CHECK(msg.find("LAMBDA >= 0") != std::string::npos);
    CHECK(S.rows() == 0);
    CHECK(c.get<double>("ALPHA") == 1.0);
    CHECK(c.get<double>("BETA") == 1.0);
  }
  for (const bool explicit_n_obs : {true, false}) {
    DYNAMIC_SECTION("M_BLR_Train<BF_Linear>, "
                    << (explicit_n_obs ? "n_obs 0" : "rows - cols = 0")) {
      Context c = oalgo3_context("-1");
      c.add("ALPHA", "1");
      c.add("BETA", "1");
      BF_Linear bf;
      M_BLR_Train<BF_Linear> m(bf, c);
      phi_type Phi = zero_phi(explicit_n_obs ? 5 : 2);
      t_type t = zero_t(explicit_n_obs ? 5 : 2);
      std::string msg = "<no exception>";
      try {
        if (explicit_n_obs)
          m.train(Phi, t, 0);
        else
          m.train(Phi, t);
      } catch (const std::runtime_error &e) {
        msg = e.what();
      }
      INFO(msg);
      CHECK(msg.find("LAMBDA -1") != std::string::npos);
      CHECK(msg.find("0 observations") != std::string::npos);
      CHECK_FALSE(m.trained);
      CHECK_FALSE(m.has_covariance());
      CHECK(m.Sigma.rows() == 0);
    }
  }
}

/* Pins: the reason a fit leaves no weight covariance
 * (covariance_unavailable(), in get_weights_uncertainty()'s message) says
 * what to change: after an OALGO 1 fit "fit with OALGO 3"; after M_KRR
 * implementation 2 (train2) "MODEL M_KRR <kernel> 1"; after an OALGO 3 fit
 * with 2 observations for 4 weights at LAMBDA 0 (gamma = 2, the noise
 * variance undefined) "2 observations" and "add data with a non-zero
 * weight"; with every weight pinned and no observation (FIXINDEX 1-2, n_obs
 * 0) "no observations" and the same remedy. tadah train --uncertainty
 * gives the reason as its remedy.
 * Before the fix: the reasons said why, not what to change ("OALGO 1
 * computes no weight covariance (only OALGO 3, the SVD solve, does)", "...
 * so its noise variance is undefined"), and tadah train --uncertainty
 * appended "Add observations (data with a non-zero weight)" to every one,
 * M_KRR implementation 2's included. */
TEST_CASE("A fit without a weight covariance says what to change",
          "[weight_covariance][uncertainty]") {
  auto message = [](auto &m) {
    try {
      m.get_weights_uncertainty();
    } catch (const std::exception &e) {
      return std::string(e.what());
    }
    return std::string("<no exception>");
  };
  {
    Context c = oalgo3_context("1e-3");
    c.remove("OALGO");
    c.add("OALGO", "1");
    BF_Linear bf;
    M_BLR_Train<BF_Linear> m(bf, c);
    const LMat A = make_A(9, 4, 0.61);
    phi_type Phi = with_augmentation(A);
    t_type t = target_with_augmentation(make_T(A), 4);
    m.train(Phi, t);
    const std::string msg = message(m);
    INFO(msg);
    CHECK(msg.find("fit with OALGO 3") != std::string::npos);
  }
  {
    Context c = oalgo3_context("1e-3");
    Kern_Linear k;
    M_KRR_Train<Kern_Linear> m(k, c);
    Matrix B(4, 3);
    for (size_t j = 0; j < 3; ++j)
      for (size_t i = 0; i < 4; ++i) B(i, j) = (i == j ? 2.0 : 0.1 * (i + j));
    t_type T(3);
    for (size_t j = 0; j < 3; ++j) T[j] = 0.5 * static_cast<double>(j + 1);
    m.train2(B, T);
    const std::string msg = message(m);
    INFO(msg);
    CHECK(msg.find("MODEL M_KRR <kernel> 1") != std::string::npos);
  }
  {
    Context c = oalgo3_context("0");
    BF_Linear bf;
    M_BLR_Train<BF_Linear> m(bf, c);
    const LMat A = make_A(2, 4, 0.47);
    phi_type Phi = with_augmentation(A);
    t_type t = target_with_augmentation(make_T(A), 4);
    m.train(Phi, t);
    const std::string msg = message(m);
    INFO(msg);
    CHECK(msg.find("2 observations") != std::string::npos);
    CHECK(msg.find("add data with a non-zero weight") != std::string::npos);
  }
  {
    Context c = oalgo3_context("0");
    c.add("FIXINDEX", "1-2");
    c.add("FIXWEIGHT", "0.1 0.2");
    BF_Linear bf;
    M_BLR_Train<BF_Linear> m(bf, c);
    phi_type Phi(5, 2);
    Phi.set_zero();
    t_type t(5);
    for (size_t i = 0; i < 5; ++i) t[i] = 0.0;
    m.train(Phi, t, 0);
    const std::string msg = message(m);
    INFO(msg);
    CHECK(msg.find("no observations") != std::string::npos);
    CHECK(msg.find("add data with a non-zero weight") != std::string::npos);
  }
}

/* Pins: OLS::solveWithPseudoinverse of a wide 3 x 5 matrix (m < n, so
 * V^T has min(m,n) = 3 rows) returns the minimum-norm solution A^T (A
 * A^T)^-1 b at lambda 0 and the ridge solution A^T (A A^T + lambda I)^-1 b
 * at lambda 0.1 (1e-10 relative).
 * Before the fix: it multiplied by V^T with a leading dimension of n = 5
 * instead of min(m,n) = 3, which reads V^T's 15 values with the wrong
 * stride and past their end, so the weights were not the solution (lambda
 * 0: 0.103, -0.235, 3e-38, -2e-38 for weights 2 to 5 against 0.337, 0.578,
 * 0.114, -0.055; a latent defect: the design matrix always has at least as
 * many rows as columns). */
TEST_CASE("solveWithPseudoinverse of a wide matrix: the minimum-norm "
          "solution",
          "[weight_covariance][ols]") {
  const LMat A = make_A(3, 5, 0.83);
  std::vector<LD> b(3);
  for (size_t i = 0; i < 3; ++i) b[i] = 0.4L + 0.3L * static_cast<LD>(i);
  for (const LD lambda : {0.0L, 0.1L}) {
    DYNAMIC_SECTION("lambda " << static_cast<double>(lambda)) {
      LMat AAt(3, std::vector<LD>(3, 0.0L));
      for (size_t a = 0; a < 3; ++a)
        for (size_t c2 = 0; c2 < 3; ++c2)
          for (size_t k = 0; k < 5; ++k) AAt[a][c2] += A[a][k] * A[c2][k];
      for (size_t a = 0; a < 3; ++a) AAt[a][a] += lambda;
      const LMat inv = inverse(AAt);
      std::vector<LD> y(3, 0.0L), x(5, 0.0L);
      for (size_t a = 0; a < 3; ++a)
        for (size_t c2 = 0; c2 < 3; ++c2) y[a] += inv[a][c2] * b[c2];
      for (size_t k = 0; k < 5; ++k)
        for (size_t a = 0; a < 3; ++a) x[k] += A[a][k] * y[a];

      phi_type M(3, 5);
      for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 5; ++j) M(i, j) = static_cast<double>(A[i][j]);
      t_type B(3);
      for (size_t i = 0; i < 3; ++i) B[i] = static_cast<double>(b[i]);
      SVD svd;
      svd.compute(M);
      t_type w(5);
      OLS ols;
      ols.solveWithPseudoinverse(3, 5, B, w, svd.getSVDWorkspace(), 1e-8,
                                 static_cast<double>(lambda));
      for (size_t k = 0; k < 5; ++k) {
        INFO("k " << k << ": " << w[k] << " against "
                  << static_cast<double>(x[k]));
        CHECK(std::fabs(w[k] - static_cast<double>(x[k])) <
              1e-10 * std::fabs(static_cast<double>(x[k])));
      }
    }
  }
}
