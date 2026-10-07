#ifndef LINEAR_REGRESSOR_H
#define LINEAR_REGRESSOR_H

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/core/intercept_layout.h>
#include <tadah/core/lapack.h>
#include <tadah/core/number_format.h>
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

  /* A design matrix laid out as DesignMatrix lays it out: its data rows,
   * then one ridge row per column (zero with LAMBDA 0). */
  static size_t data_rows(const tadah::core::phi_type &Phi) {
    return Phi.rows() >= Phi.cols() ? Phi.rows() - Phi.cols() : Phi.rows();
  }

  /* True for a value that is not NaN or +-inf; a test on the bits, which
   * -ffast-math does not optimise away. */
  static bool finite_bits(double v) {
    std::uint64_t u = 0;
    std::memcpy(&u, &v, sizeof u);
    return (u & UINT64_C(0x7ff0000000000000)) != UINT64_C(0x7ff0000000000000);
  }

  /* Refuses a design matrix or target with a value that is not finite,
   * naming its row and column, before LAPACK is called: on NaN or inf
   * LAPACK stops with an argument or convergence error (DGELSD info -4,
   * DGELS -5, DGESDD 3), DPOTRF calls the matrix not positive definite, and
   * an infinite target is fitted without an error. */
  static void require_finite(const tadah::core::phi_type &Phi,
                             const tadah::core::t_type &T) {
    const size_t data = data_rows(Phi);
    const std::string layout =
        std::to_string(data) + " data rows, then " +
        std::to_string(Phi.rows() - data) + " ridge rows";
    for (size_t j = 0; j < Phi.cols(); ++j)
      for (size_t i = 0; i < Phi.rows(); ++i)
        if (!finite_bits(Phi(i, j)))
          throw std::runtime_error(
              "LinearRegressor::train: the design matrix holds a value that "
              "is not finite (nan or inf) at row " + std::to_string(i + 1) +
              ", column " + std::to_string(j + 1) + " (of " + layout +
              ", and " + std::to_string(Phi.cols()) +
              " columns): a descriptor or basis function gave it; no weights "
              "were fitted");
    for (size_t i = 0; i < T.size() && i < Phi.rows(); ++i)
      if (!finite_bits(T[i]))
        throw std::runtime_error(
            "LinearRegressor::train: the target vector holds a value that is "
            "not finite (nan or inf) at row " + std::to_string(i + 1) +
            " (of " + layout + "); no weights were fitted");
  }

  /* The message of a fit OALGO 2 or 4 cannot make: what `e` found, the
   * fit's data rows, weights and LAMBDA, and the remedies. */
  [[noreturn]] static void throw_rank_deficient(int oalgo,
                                                const RankDeficientError &e,
                                                const tadah::core::phi_type &Phi,
                                                double lambda) {
    const size_t n = Phi.cols();
    const size_t data = data_rows(Phi);
    const std::string fit = std::to_string(data) + " data rows for " +
                            std::to_string(n) + " weights, LAMBDA " +
                            tadah::core::format_roundtrip(lambda);
    std::string msg;
    if (oalgo == 2) {
      msg = "OALGO 2 (DGELS, QR) cannot fit these " + std::to_string(n) +
            " weights: the design matrix (" + fit +
            ") does not have full column rank: its column " +
            std::to_string(e.column()) +
            " is zero or a linear combination of the columns before it" +
            (data < n ? " (fewer data rows than weights)" : "") +
            ", so the least-squares weights are not unique.";
    } else {
      msg = "OALGO 4 (normal equations, Cholesky) cannot fit these " +
            std::to_string(n) +
            " weights: the matrix Phi^T Phi of the design matrix (" + fit +
            ") is not positive definite (its leading minor of order " +
            std::to_string(e.column()) +
            " is not), so the design matrix does not have full column rank "
            "or is too ill-conditioned for the normal equations.";
    }
    throw std::runtime_error(
        msg +
        " Use OALGO 1 (DGELSD) or OALGO 3 (SVD), which give the minimum-norm "
        "weights, LAMBDA > 0 (a ridge penalty makes every column it "
        "penalises full rank; it does not penalise the intercept columns of "
        "BIAS), or more data (structures, FORCE or STRESS rows).");
  }


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
   * @param n_intercept The intercept columns of Phi, columns 0 ..
   *        n_intercept-1 (Function_Base::intercept_columns() of the
   *        function that built Phi): the LAMBDA ridge penalty leaves them
   *        out, and OALGO 3, whose ridge penalises every column alike,
   *        refuses LAMBDA != 0 when there is one.
   */
public:
  static void train(tadah::core::Context &context,
                    tadah::core::phi_type &Phi,
                    tadah::core::t_type &T,
                    tadah::core::t_type &weights,
                    tadah::core::Matrix &Sigma,
                    tadah::models::memory::IModelsWorkspaceManager &workspaceManager,
                    const size_t n_intercept) {

    double lambda = context.get<double>("LAMBDA", 0);
    int oalgo = context.get<int>("OALGO", 0);
    double rcond =
        context.size("OALGO") == 2 ? context.get<double>("OALGO", 1) : 1e-8;
    weights.resize(Phi.cols());

    if (std::abs(lambda + 1) < std::numeric_limits<double>::min() &&
        oalgo != 3) {
      throw std::runtime_error("LAMBDA -1 requires OALGO 3");
    }

    const size_t n_bias = n_intercept;

    if (lambda > 0 && oalgo != 3 ) {
      // Augment the design matrix with the regularization term.
      // Skip the first n_bias rows of the augmentation block: intercept
      // columns are unpenalised (standard practice — regularising an intercept
      // biases the fit toward zero mean).
      const size_t aug_start = Phi.rows() - Phi.cols();
      size_t j = n_bias;
      for (size_t i = aug_start + n_bias; i < Phi.rows(); ++i) {
        Phi(i, j++) = std::sqrt(lambda);
      }
    }

    tadah::models::OLS ols(workspaceManager);
    // tadah::models::OLS ols;
    require_finite(Phi, T);
    if (oalgo == 1) {
      ols.solve(Phi, T, weights, rcond, tadah::models::Algorithm::GELSD);
    } else if (oalgo == 2) {
      try {
        ols.solve(Phi, T, weights, rcond, tadah::models::Algorithm::GELS);
      } catch (const RankDeficientError &e) {
        throw_rank_deficient(2, e, Phi, lambda);
      }
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
      try {
        ols.solve(Phi, T, weights, rcond,
                  tadah::models::Algorithm::NormalEquations);
      } catch (const RankDeficientError &e) {
        throw_rank_deficient(4, e, Phi, lambda);
      }
    } else {
      throw std::runtime_error("Unsupported OALGO: " + std::to_string(oalgo));
    }

  }
  /** @brief train() of a linear model: its intercept columns are the
   *  bias columns of the descriptor, tadah::core::n_intercept_columns()
   *  (Function_Base::intercept_columns() of BF_Linear and Kern_Linear). */
  static void train(tadah::core::Context &context,
                    tadah::core::phi_type &Phi,
                    tadah::core::t_type &T,
                    tadah::core::t_type &weights,
                    tadah::core::Matrix &Sigma,
                    tadah::models::memory::IModelsWorkspaceManager &workspaceManager) {
    train(context, Phi, T, weights, Sigma, workspaceManager,
          tadah::core::n_intercept_columns(context));
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

    require_finite(Phi, T);
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
   * @brief Store the covariance matrix in the configuration as SIGMA.
   *
   * SIGMA is N, then the N x N matrix in column-major order. It replaces a
   * SIGMA the configuration held (of an earlier fit, or of a potential used
   * as the configuration): the key holds the covariance of this fit only.
   *
   * @param context Context object.
   * @param Sigma Covariance matrix to be added.
   */
  static void config_add_sigma(tadah::core::Context &context, tadah::core::Matrix &Sigma) {
    tadah::core::t_type Sigma_as_vector(Sigma.data(), Sigma.cols() * Sigma.rows());
    context.remove("SIGMA");
    context.add("SIGMA", Sigma.rows());
    context.add<tadah::core::t_type>("SIGMA", Sigma_as_vector);
  }

  /**
   * @brief Read and reconstruct the covariance matrix from the configuration.
   *
   * Retrieves the Sigma matrix stored in the configuration: SIGMA holds N,
   * then N*N values.
   *
   * @param context Context object.
   * @param Sigma tadah::core::Matrix to store the reconstructed Sigma.
   * @throws std::runtime_error if Sigma is not computed, or if SIGMA does
   *         not hold N*N+1 values (the message names N and both counts).
   */
  static void read_sigma(const tadah::core::Context &context, tadah::core::Matrix &Sigma) {
    using V = std::vector<double>;
    size_t N;
    try {
      N = context.get<size_t>("SIGMA");
    } catch (const std::runtime_error &error) {
      throw std::runtime_error("\nSigma matrix is not computed.\n");
    }

    const size_t stored = context.size("SIGMA");
    if (stored != N * N + 1)
      throw std::runtime_error(
          "SIGMA: N = " + std::to_string(N) + " needs N*N+1 = " +
          std::to_string(N * N + 1) + " values (N and the N x N matrix), "
          "but SIGMA has " + std::to_string(stored) + "; retrain the "
          "potential with OALGO 3 to compute it.");
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
