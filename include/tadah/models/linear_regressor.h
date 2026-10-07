#ifndef LINEAR_REGRESSOR_H
#define LINEAR_REGRESSOR_H

#include <cstddef>
#include <limits>
#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/core/intercept_layout.h>
#include <tadah/core/lapack.h>
#include <tadah/models/ea.h>
#include <tadah/models/factor_handle.h>
#include <tadah/models/ols.h>
#include <tadah/models/svd.h>
#include <tadah/models/memory/imodels_workspace_manager.h>

#include <vector>

namespace tadah {
namespace models {

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
   */
public:
  static void train(tadah::core::Context &context,
                    tadah::core::phi_type &Phi,
                    tadah::core::t_type &T,
                    tadah::core::t_type &weights,
                    tadah::core::Matrix &Sigma,
                    tadah::models::memory::IModelsWorkspaceManager &workspaceManager,
                    const std::vector<size_t> &intercepts) {

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
      svd.compute(Phi);
      if (lambda < 0) {
        double alpha = context.get<double>("ALPHA");
        double beta = context.get<double>("BETA");
        EA ea(context, svd, T);
        ea.run(alpha, beta);
        lambda = alpha / beta;
        context.remove("ALPHA");
        context.remove("BETA");
        context.add("ALPHA", alpha);
        context.add("BETA", beta);
      }

      Sigma = get_sigma(svd, lambda);
      config_add_sigma(context, Sigma);

      core::Logger::instance().log(tadah::core::LogLevel::INFO,
          "LinearRegressor::train -> Regularization lambda: ", lambda);

      ols.solveWithPseudoinverse(Phi.rows(), Phi.cols(), T, weights,
                                 svd.getSVDWorkspace(), rcond, lambda);
    } else if (oalgo == 4) {
      ols.solve(Phi, T, weights, rcond,
                tadah::models::Algorithm::NormalEquations);
    } else {
      throw std::runtime_error("Unsupported OALGO: " + std::to_string(oalgo));
    }

  }
  /** @brief train() of a linear model: its intercept columns are the
   *  intercept columns of the descriptor, 0 .. N-1 with N =
   *  tadah::core::n_intercept_columns() (Function_Base::intercept_columns()
   *  of BF_Linear and Kern_Linear). */
  static void train(tadah::core::Context &context,
                    tadah::core::phi_type &Phi,
                    tadah::core::t_type &T,
                    tadah::core::t_type &weights,
                    tadah::core::Matrix &Sigma,
                    tadah::models::memory::IModelsWorkspaceManager &workspaceManager) {
    std::vector<size_t> intercepts(tadah::core::n_intercept_columns(context));
    for (size_t k = 0; k < intercepts.size(); ++k) intercepts[k] = k;
    train(context, Phi, T, weights, Sigma, workspaceManager, intercepts);
  }
  /** @brief train() of a linear model with a temporary workspace manager. */
  static void train(tadah::core::Context &context, tadah::core::phi_type &Phi, tadah::core::t_type &T, tadah::core::t_type &weights,
                    tadah::core::Matrix &Sigma) {
    tadah::models::memory::ModelsWorkspaceManager workspaceManager;
    train(context, Phi, T, weights, Sigma, workspaceManager);
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
   * @brief Compute and return the covariance matrix.
   *
   * Uses the SVD to calculate the covariance matrix with regularization.
   *
   * @param svd SVD object of the design matrix.
   * @param lambda Regularization parameter.
   * @return Covariance matrix.
   */
  static tadah::core::Matrix get_sigma(tadah::models::SVD &svd, double lambda) {
    double *VT = svd.getVT();
    double *S = svd.getS();
    int n = static_cast<int>(svd.cols());

    tadah::core::Matrix Sigma((size_t)n, (size_t)n);
    Sigma.set_zero();
    double *sigma = &Sigma.data()[0];

    for (int i = 0; i < n; ++i) {
      double ridge = 1.0 / (S[i] * S[i] + lambda);
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          sigma[j + k * n] +=
              VT[j + i * n] * ridge * VT[k + i * n]; // Column-major
        }
      }
    }
    return Sigma;
  }

  /**
   * @brief Add the covariance matrix to the configuration.
   *
   * Converts the Sigma matrix to a vector and stores it in the configuration.
   *
   * @param context Context object.
   * @param Sigma Covariance matrix to be added.
   */
  static void config_add_sigma(tadah::core::Context &context, tadah::core::Matrix &Sigma) {
    tadah::core::t_type Sigma_as_vector(Sigma.data(), Sigma.cols() * Sigma.rows());
    context.add("SIGMA", Sigma.rows());
    context.add<tadah::core::t_type>("SIGMA", Sigma_as_vector);
  }

  /**
   * @brief Read and reconstruct the covariance matrix from the configuration.
   *
   * Retrieves the Sigma matrix stored in the configuration.
   *
   * @param context Context object.
   * @param Sigma tadah::core::Matrix to store the reconstructed Sigma.
   * @throws std::runtime_error if Sigma is not computed.
   */
  static void read_sigma(const tadah::core::Context &context, tadah::core::Matrix &Sigma) {
    using V = std::vector<double>;
    size_t N;
    try {
      N = context.get<size_t>("SIGMA");
    } catch (const std::runtime_error &error) {
      throw std::runtime_error("\nSigma matrix is not computed.\n");
    }

    V Sigma_as_vector(N * N + 1);
    context.get<V>("SIGMA", Sigma_as_vector);
    Sigma_as_vector.erase(Sigma_as_vector.begin());
    Sigma.resize(N, N);
    for (size_t c = 0; c < N; ++c)
      for (size_t r = 0; r < N; ++r)

        Sigma(r, c) = Sigma_as_vector.at(N * c + r);
  }
};
} // namespace models
} // namespace tadah
#endif
