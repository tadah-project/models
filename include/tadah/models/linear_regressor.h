#ifndef LINEAR_REGRESSOR_H
#define LINEAR_REGRESSOR_H

#include <cstddef>
#include <limits>
#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/core/intercept_layout.h>
#include <tadah/core/lapack.h>
#include <tadah/core/logger.h>
#include <tadah/core/utils.h>
#include <tadah/models/ea.h>
#include <tadah/models/factor_handle.h>
#include <tadah/models/ols.h>
#include <tadah/models/svd.h>
#include <tadah/models/memory/imodels_workspace_manager.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace tadah {
namespace models {

/**
 * @brief What LinearRegressor::train computes besides the weights: the
 *        weight covariance of an OALGO 3 fit, or why there is none.
 */
struct CovarianceFit {
  /** Whether the Sigma out-parameter of train() holds the covariance of
   *  this fit's weights (else it is 0 x 0). */
  bool computed = false;
  /** s^2 = RSS / (n_obs - gamma), the noise variance per unit row weight
   *  (0 when not computed). */
  double noise_variance = 0.0;
  /** The weighted residual sum of squares of the weights (OALGO 3). */
  double rss = 0.0;
  /** The effective number of weights, sum over the kept directions of
   *  s_i^2 / (s_i^2 + lambda) (OALGO 3). */
  double gamma = 0.0;
  /** The number of observations the fit was given. */
  std::size_t n_obs = 0;
  /** The directions of the weights the solve leaves out: s_i at most the
   *  rank threshold, and the n - min(m, n) directions of a matrix with
   *  fewer rows than columns (OALGO 3). */
  std::size_t n_dropped = 0;
  /** Why no covariance was computed (empty when it was). */
  std::string why_not;
};

/**
 * @class LinearRegressor
 * @brief Performs linear regression using Singular Value Decomposition (SVD).
 *
 * Utilizes SVD for training, supporting both ordinary least squares and ridge
 * regression.
 */
class LinearRegressor {

  /**
   * @brief Train the model using the provided configuration and data matrices.
   *
   * This method supports both ordinary least squares and ridge regression based
   * on the lambda value.
   *
   * @param context Context object containing training parameters.
   * @param Phi Design matrix.
   * @param T Target matrix.
   * @param weights Vector to store the resultant weights.
   * @param Sigma tadah::core::Matrix to store the covariance matrix.
   * @param workspaceManager Workspace provider.
   * @param intercepts The intercept columns of Phi, 0-based indices into
   *        its columns (Function_Base::intercept_columns() of the function
   *        that built Phi, mapped to Phi's columns when FIXINDEX removed
   *        some): the LAMBDA ridge penalty leaves them out and penalises
   *        every other column, and OALGO 3, whose ridge penalises every
   *        column alike, refuses LAMBDA != 0 when there is one.
   * @param n_obs The number of observations: the rows of Phi with a
   *        non-zero weight, without its ridge-augmentation rows (the
   *        DesignMatrix appends one per column) and without the rows of
   *        weight 0 (zero rows of the weighted matrix). OALGO 3 estimates
   *        the noise variance with it, and LAMBDA -1 the evidence.
   * @return With OALGO 3: the weight covariance in @p Sigma, s^2 times
   *         sum over the kept directions of v_i v_i^T / (s_i^2 + lambda)
   *         (get_sigma), s^2 = RSS / (n_obs - gamma), and the numbers it
   *         came from; none (Sigma 0 x 0, a WARNING) when n_obs <= gamma.
   *         Any other OALGO computes none. Nothing is written into
   *         @p context but the ALPHA and BETA of LAMBDA -1.
   */
public:
  static CovarianceFit train(tadah::core::Context &context,
                    tadah::core::phi_type &Phi,
                    tadah::core::t_type &T,
                    tadah::core::t_type &weights,
                    tadah::core::Matrix &Sigma,
                    tadah::models::memory::IModelsWorkspaceManager &workspaceManager,
                    const std::vector<size_t> &intercepts,
                    std::size_t n_obs) {

    double lambda = context.get<double>("LAMBDA", 0);
    int oalgo = context.get<int>("OALGO", 0);
    double rcond =
        context.size("OALGO") == 2 ? context.get<double>("OALGO", 1) : 1e-8;
    weights.resize(Phi.cols());

    if (std::abs(lambda + 1) < std::numeric_limits<double>::min() &&
        oalgo != 3) {
      throw std::runtime_error("LAMBDA -1 requires OALGO 3");
    }

    const size_t n_bias = intercepts.size();
    for (size_t j : intercepts)
      if (j >= Phi.cols())
        throw std::logic_error(
            "LinearRegressor::train: intercept column " + std::to_string(j) +
            " of a design matrix of " + std::to_string(Phi.cols()) +
            " columns");

    if (lambda > 0 && oalgo != 3 ) {
      // Augment the design matrix with the regularization term: row
      // aug_start + j penalises column j, except for the intercept columns,
      // which stay unpenalised (standard practice: regularising an
      // intercept biases the fit toward zero mean).
      const size_t aug_start = Phi.rows() - Phi.cols();
      std::vector<bool> intercept(Phi.cols(), false);
      for (size_t j : intercepts) intercept[j] = true;
      for (size_t j = 0; j < Phi.cols(); ++j)
        if (!intercept[j]) Phi(aug_start + j, j) = std::sqrt(lambda);
    }

    CovarianceFit fit;
    fit.n_obs = n_obs;
    fit.why_not = "OALGO " + std::to_string(oalgo) +
                  " computes no weight covariance (only OALGO 3, the SVD "
                  "solve, does)";
    Sigma = tadah::core::Matrix();

    tadah::models::OLS ols(workspaceManager);
    // tadah::models::OLS ols;
    if (oalgo == 1) {
      ols.solve(Phi, T, weights, rcond, tadah::models::Algorithm::GELSD);
    } else if (oalgo == 2) {
      ols.solve(Phi, T, weights, rcond, tadah::models::Algorithm::GELS);
    } else if (oalgo == 3) {
      // SVD path applies a uniform ridge across all singular values; there is
      // no per-column knob, so intercept-aware regularisation is not
      // supported here. Reject loudly and point at the alternatives.
      // TODO: implement Path B (residualise intercepts, run SVD on the
      //       descriptor block only) so OALGO=3 + BIAS=true becomes valid.
      if (n_bias > 0 && lambda != 0.0) {
        throw std::runtime_error(
            "OALGO=3 (SVD with uniform ridge) does not support intercept "
            "exclusion from regularisation. The model has " +
            std::to_string(n_bias) + " intercept column(s) (BIAS=true with "
            "ATOM) and LAMBDA != 0: use OALGO=1, 2 or 4 instead.");
      }
      tadah::models::SVD svd;
      svd.compute(Phi);  // Phi is overwritten: the residuals come from the SVD
      if (lambda < 0) {
        double alpha = context.get<double>("ALPHA");
        double beta = context.get<double>("BETA");
        EA ea(context, svd, T, n_obs);
        ea.run(alpha, beta);
        lambda = alpha / beta;
        context.remove("ALPHA");
        context.remove("BETA");
        context.add("ALPHA", alpha);
        context.add("BETA", beta);
      }

      core::Logger::instance().log(tadah::core::LogLevel::INFO,
          "LinearRegressor::train -> Regularization lambda: ", lambda);

      ols.solveWithPseudoinverse(Phi.rows(), Phi.cols(), T, weights,
                                 svd.getSVDWorkspace(), rcond, lambda);
      fit = covariance_of_fit(svd, T, weights, lambda, rcond, n_obs, Sigma);
      if (fit.n_dropped > 0)
        core::Logger::instance().log(tadah::core::LogLevel::WARNING,
            "LinearRegressor (OALGO 3): " + std::to_string(fit.n_dropped) +
            " of " + std::to_string(Phi.cols()) + " singular directions of "
            "the design matrix are at or below the rank threshold (OALGO's "
            "rcond " + core::to_string(rcond) + " times the largest "
            "singular value): the solve leaves them out of the weights, and "
            "the weight covariance WCOV leaves them out too (no uncertainty "
            "along them).");
      if (!fit.computed)
        core::Logger::instance().log(tadah::core::LogLevel::WARNING,
            "LinearRegressor (OALGO 3): no weight covariance: " +
            fit.why_not + ".");
    } else if (oalgo == 4) {
      ols.solve(Phi, T, weights, rcond,
                tadah::models::Algorithm::NormalEquations);
    } else {
      throw std::runtime_error("Unsupported OALGO: " + std::to_string(oalgo));
    }
    return fit;
  }
  /** @brief train() with the observations of a design matrix that has its
   *  ridge-augmentation rows and no rows of weight 0: n_obs = rows -
   *  cols (0 when Phi has no more rows than columns). */
  static CovarianceFit train(tadah::core::Context &context,
                    tadah::core::phi_type &Phi,
                    tadah::core::t_type &T,
                    tadah::core::t_type &weights,
                    tadah::core::Matrix &Sigma,
                    tadah::models::memory::IModelsWorkspaceManager &workspaceManager,
                    const std::vector<size_t> &intercepts) {
    return train(context, Phi, T, weights, Sigma, workspaceManager, intercepts,
                 augmented_observations(Phi));
  }
  /** @brief train() of a linear model: its intercept columns are the
   *  intercept columns of the descriptor, 0 .. N-1 with N =
   *  tadah::core::n_intercept_columns() (Function_Base::intercept_columns()
   *  of BF_Linear and Kern_Linear); n_obs = rows - cols. */
  static CovarianceFit train(tadah::core::Context &context,
                    tadah::core::phi_type &Phi,
                    tadah::core::t_type &T,
                    tadah::core::t_type &weights,
                    tadah::core::Matrix &Sigma,
                    tadah::models::memory::IModelsWorkspaceManager &workspaceManager) {
    std::vector<size_t> intercepts(tadah::core::n_intercept_columns(context));
    for (size_t k = 0; k < intercepts.size(); ++k) intercepts[k] = k;
    return train(context, Phi, T, weights, Sigma, workspaceManager, intercepts);
  }
  /** @brief train() of a linear model with a temporary workspace manager. */
  static CovarianceFit train(tadah::core::Context &context, tadah::core::phi_type &Phi, tadah::core::t_type &T, tadah::core::t_type &weights,
                    tadah::core::Matrix &Sigma) {
    tadah::models::memory::ModelsWorkspaceManager workspaceManager;
    return train(context, Phi, T, weights, Sigma, workspaceManager);
  }

  /** @brief The observations of a design matrix with its cols
   *  ridge-augmentation rows (the DesignMatrix contract, which the LAMBDA
   *  ridge of OALGO 1 and 2 writes into): rows - cols, or 0 when Phi has
   *  no more rows than columns. A caller that knows the rows of weight 0
   *  passes their exclusion to train() itself. */
  static std::size_t augmented_observations(const tadah::core::phi_type &Phi) {
    return Phi.rows() > Phi.cols() ? Phi.rows() - Phi.cols() : 0;
  }

  /**
   * @brief Train + cache the Cholesky factor for IFT-based gradients.
   *
   * Mirrors @ref train but additionally exposes the cached
   * factorisation through @p factor for HPO analytical/hybrid
   * gradients to consume via @ref backsolve. Restricted to OALGO=4
   * (NormalEquations) with LAMBDA >= 0 — the supported configuration
   * for Phase-2A IFT.
   *
   * The forward path of @ref train is byte-identical and unaffected by
   * this overload; HPO callers that want analytical gradients must
   * opt in by calling this method instead.
   *
   * @param context     Training parameters; must declare OALGO=4.
   * @param Phi         Design matrix (m+n rows; ridge-augmentation
   *                    block written in place if LAMBDA>0).
   * @param T           Target vector.
   * @param weights     Receives the OLS solution.
   * @param Sigma       Unused for OALGO=4 but kept for API symmetry.
   * @param factor      Output: receives the cached Cholesky factor.
   * @param workspaceManager  Workspace provider.
   */
  static void train_with_factor(
      tadah::core::Context &context,
      tadah::core::phi_type &Phi,
      tadah::core::t_type &T,
      tadah::core::t_type &weights,
      tadah::core::Matrix &Sigma,
      FactorHandle &factor,
      tadah::models::memory::IModelsWorkspaceManager &workspaceManager) {

    (void)Sigma;  // Reserved for future BLR-IFT; OALGO=4 leaves it untouched.

    const double lambda = context.get<double>("LAMBDA", 0);
    const int    oalgo  = context.get<int>("OALGO", 0);
    const double rcond  =
        context.size("OALGO") == 2 ? context.get<double>("OALGO", 1) : 1e-8;
    (void)rcond;  // NormalEquations does not consume rcond.

    if (oalgo != 4) {
      throw std::runtime_error(
          "LinearRegressor::train_with_factor requires OALGO=4 "
          "(NormalEquations); got OALGO=" + std::to_string(oalgo) +
          ". The IFT chain rule for HPO analytical gradients reuses the "
          "cached Cholesky factor, which is only produced by the "
          "NormalEquations path. Use OALGO=4 with GRAD_STRATEGY="
          "ANALYTICAL/HYBRID, or fall back to GRAD_STRATEGY=NUMERICAL "
          "with any OALGO.");
    }
    if (std::abs(lambda + 1) < std::numeric_limits<double>::min()) {
      throw std::runtime_error(
          "LinearRegressor::train_with_factor: LAMBDA=-1 (BLR/EA) is "
          "not supported under analytical gradients in Phase 2A. The "
          "α/β co-optimisation nests inside the solve and changes the "
          "IFT chain rule; revisit in Phase 2C if needed. Use a fixed "
          "LAMBDA >= 0, or GRAD_STRATEGY=NUMERICAL.");
    }
    if (lambda < 0) {
      throw std::runtime_error(
          "LinearRegressor::train_with_factor: LAMBDA must be >= 0 "
          "(got " + std::to_string(lambda) + ").");
    }

    weights.resize(Phi.cols());

    const std::size_t n_bias = tadah::core::n_intercept_columns(context);

    // Fill the ridge-augmentation block exactly as the existing
    // train() does — same byte pattern so a parity test is meaningful.
    if (lambda > 0) {
      const std::size_t aug_start = Phi.rows() - Phi.cols();
      std::size_t j = n_bias;
      for (std::size_t i = aug_start + n_bias; i < Phi.rows(); ++i) {
        Phi(i, j++) = std::sqrt(lambda);
      }
    }

    tadah::models::OLS ols(workspaceManager);
    ols.solveAndCacheFactor(Phi, T, weights, factor, n_bias);
  }

  /// Convenience overload using a temporary workspace manager.
  static void train_with_factor(
      tadah::core::Context &context,
      tadah::core::phi_type &Phi,
      tadah::core::t_type &T,
      tadah::core::t_type &weights,
      tadah::core::Matrix &Sigma,
      FactorHandle &factor) {
    tadah::models::memory::ModelsWorkspaceManager workspaceManager;
    train_with_factor(context, Phi, T, weights, Sigma, factor, workspaceManager);
  }

  /**
   * @brief Solve U^T U x = rhs using the cached Cholesky factor.
   *
   * Used by HPO analytical/hybrid gradients to compute dW_star/dtheta
   * via the implicit-function-theorem chain rule. With fixed lambda:
   *
   *     (Phi^T Phi + lambda I_no_bias) . dW_star/dtheta
   *         = (dPhi^T/dtheta)(T - Phi W_star)
   *           - Phi^T (dPhi/dtheta) W_star
   *
   * One back-sub per theta-partial; the factor is shared across all
   * calls within a single HPO gradient evaluation.
   *
   * @param f    Cached factor (must be valid and NormalEquations).
   * @param rhs  Right-hand side; size must equal f.n.
   * @param sol  Output; resized to f.n if needed.
   */
  static void backsolve(const FactorHandle &f,
                        const tadah::core::t_type &rhs,
                        tadah::core::t_type &sol) {
    if (!f.valid) {
      throw std::runtime_error(
          "LinearRegressor::backsolve: FactorHandle is invalid. The "
          "factor was either never populated or was invalidated (e.g. "
          "after a descriptor-size change in HPO).");
    }
    if (f.algo != tadah::models::Algorithm::NormalEquations) {
      throw std::runtime_error(
          "LinearRegressor::backsolve: only NormalEquations factors are "
          "supported in Phase 2A; got algo " +
          tadah::models::to_string_algorithm(f.algo));
    }
    if (rhs.size() != static_cast<std::size_t>(f.n)) {
      throw std::invalid_argument(
          "LinearRegressor::backsolve: rhs size (" +
          std::to_string(rhs.size()) + ") does not match factor "
          "dimension (" + std::to_string(f.n) + ").");
    }

    // dpotrs_ overwrites the rhs with the solution, so copy first.
    sol.resize(static_cast<std::size_t>(f.n));
    for (int i = 0; i < f.n; ++i) sol[i] = rhs[i];

    // dpotrs_ takes a non-const factor pointer for LAPACK signature
    // compatibility but does not modify it. Same const-cast pattern as
    // the existing OLS::solveImpl call site.
    char uplo = 'U';
    int  n    = f.n;
    int  nrhs = 1;
    int  info = 0;
    double *factor_ptr = const_cast<double *>(f.AtA.data());
    tadah::core::dpotrs_(&uplo, &n, &nrhs, factor_ptr, &n,
                         sol.ptr(), &n, &info);
    if (info != 0) {
      throw std::runtime_error(
          "LinearRegressor::backsolve: DPOTRS failed; info = " +
          std::to_string(info));
    }
  }

  /**
   * @brief The unscaled weight covariance of the SVD solve:
   *        sum over the kept directions i of v_i v_i^T / (s_i^2 + lambda).
   *
   * A = U S V^T is the (weighted) design matrix of @p svd (m x n, V^T
   * min(m,n) x n with leading dimension min(m,n)); v_i is row i of V^T. A
   * direction is kept when s_i is above OLS::pseudoinverse_threshold(S,
   * min(m,n), @p rcond), the rule solveWithPseudoinverse leaves the other
   * directions out of the weights with (one rule for every lambda: the
   * covariance is that of the weights the solve returns). With every
   * direction kept and m >= n it is (A^T A + lambda I)^-1. Computed as
   * W^T W (dsyrk), W the kept rows of V^T scaled by 1/sqrt(s_i^2+lambda).
   *
   * @param n_dropped Receives the number of directions left out: the
   *        directions at or below the threshold, plus n - min(m,n).
   * @return The n x n matrix, symmetric.
   */
  static tadah::core::Matrix get_sigma(tadah::models::SVD &svd, double lambda,
                                       double rcond, std::size_t &n_dropped) {
    const double *VT = svd.getVT();
    const double *S = svd.getS();
    const int m = svd.rows();
    const int n = svd.cols();
    const int k = std::min(m, n);
    const double thresh = OLS::pseudoinverse_threshold(S, k, rcond);
    std::vector<int> kept;
    for (int i = 0; i < k; ++i)
      if (S[i] > thresh) kept.push_back(i);
    n_dropped = static_cast<std::size_t>(n - static_cast<int>(kept.size()));

    tadah::core::Matrix Sigma(static_cast<size_t>(n), static_cast<size_t>(n));
    Sigma.set_zero();
    int kk = static_cast<int>(kept.size());
    if (kk == 0 || n == 0) return Sigma;
    std::vector<double> W(static_cast<size_t>(kk) * static_cast<size_t>(n));
    for (int t = 0; t < kk; ++t) {
      const int i = kept[t];
      const double f = 1.0 / std::sqrt(S[i] * S[i] + lambda);
      for (int j = 0; j < n; ++j) W[t + j * kk] = f * VT[i + j * k];
    }
    char uplo = 'U', trans = 'T';
    double one = 1.0, zero = 0.0;
    int ldc = n;
    tadah::core::dsyrk_(&uplo, &trans, &ldc, &kk, &one, W.data(), &kk, &zero,
                        Sigma.ptr(), &ldc);
    for (int c = 0; c < n; ++c)
      for (int r = c + 1; r < n; ++r) Sigma(r, c) = Sigma(c, r);
    return Sigma;
  }

  /**
   * @brief The covariance of an OALGO 3 fit's @p weights: s^2 get_sigma(),
   *        s^2 = RSS / (n_obs - gamma).
   *
   * RSS = |T - A w|^2, with A w computed from @p svd (SVDMultiplier: the
   * solve overwrote A); the augmentation rows and the rows of weight 0
   * are zero in A and in T, so they add nothing. gamma = sum over the kept
   * directions of s_i^2 / (s_i^2 + lambda). When n_obs <= gamma the noise
   * variance is undefined: @p Sigma is 0 x 0 and the result says why.
   */
  static CovarianceFit covariance_of_fit(tadah::models::SVD &svd,
                                         const tadah::core::t_type &T,
                                         tadah::core::t_type &weights,
                                         double lambda, double rcond,
                                         std::size_t n_obs,
                                         tadah::core::Matrix &Sigma) {
    CovarianceFit fit;
    fit.n_obs = n_obs;
    const int m = svd.rows();
    const int k = std::min(m, svd.cols());
    const double *S = svd.getS();
    const double thresh = OLS::pseudoinverse_threshold(S, k, rcond);
    for (int i = 0; i < k; ++i)
      if (S[i] > thresh) fit.gamma += S[i] * S[i] / (S[i] * S[i] + lambda);
    std::vector<double> Tpred(static_cast<size_t>(m));
    SVDMultiplier(svd).compute(weights.ptr(), Tpred.data());
    for (int i = 0; i < m; ++i) {
      const double r = T[static_cast<size_t>(i)] - Tpred[static_cast<size_t>(i)];
      fit.rss += r * r;
    }
    Sigma = get_sigma(svd, lambda, rcond, fit.n_dropped);
    const double dof = static_cast<double>(n_obs) - fit.gamma;
    if (!(dof > 0.0)) {
      Sigma = tadah::core::Matrix();
      fit.why_not = "the fit has " + std::to_string(n_obs) +
                    " observations (rows with a non-zero weight) for an "
                    "effective " + core::to_string(fit.gamma) +
                    " weights (gamma), so its noise variance is undefined";
      return fit;
    }
    fit.noise_variance = fit.rss / dof;
    for (size_t c = 0; c < Sigma.cols(); ++c)
      for (size_t r = 0; r < Sigma.rows(); ++r) Sigma(r, c) *= fit.noise_variance;
    fit.computed = true;
    fit.why_not.clear();
    return fit;
  }
};
} // namespace models
} // namespace tadah
#endif
