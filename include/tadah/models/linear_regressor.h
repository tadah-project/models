#ifndef LINEAR_REGRESSOR_H
#define LINEAR_REGRESSOR_H

#include <cstddef>
#include <limits>
#include <tadah/core/config.h>
#include <tadah/core/core_types.h>
#include <tadah/models/ea.h>
#include <tadah/models/ols.h>
#include <tadah/models/svd.h>
#include <tadah/models/memory/IModelsWorkspaceManager.h>

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
   * @param config Configuration object containing training parameters.
   * @param Phi Design matrix.
   * @param T Target matrix.
   * @param weights Vector to store the resultant weights.
   * @param Sigma tadah::core::Matrix to store the covariance matrix.
   */
public:
  static void train(tadah::core::Config &config, tadah::core::phi_type &Phi, tadah::core::t_type &T, tadah::core::t_type &weights,
                    tadah::core::Matrix &Sigma, tadah::models::memory::IModelsWorkspaceManager &workspaceManager) {

    int verbose = config.get<int>("VERBOSE");
    double lambda = config.get<double>("LAMBDA", 0);
    int oalgo = config.get<int>("OALGO", 0);
    double rcond =
        config.size("OALGO") == 2 ? config.get<double>("OALGO", 1) : 1e-8;
    weights.resize(Phi.cols());

    if (std::abs(lambda + 1) < std::numeric_limits<double>::min() &&
        oalgo != 3) {
      throw std::runtime_error("LAMBDA -1 requires OALGO 3");
    }

    if (lambda > 0 && oalgo != 3 ) {
      // Augment the design matrix with the regularization term
      size_t j = 0;
      for (size_t i = Phi.rows() - Phi.cols(); i < Phi.rows(); ++i) {
        Phi(i, j++) = std::sqrt(lambda);
      }
    }

    tadah::models::OLS ols(workspaceManager);
    // tadah::models::OLS ols;
    if (oalgo == 1) {
      ols.solve(Phi, T, weights, rcond, tadah::models::Algorithm::GELSD);
    } else if (oalgo == 2) {
      ols.solve(Phi, T, weights, rcond, tadah::models::Algorithm::GELS);
    } else if (oalgo == 3) {
      tadah::models::SVD svd;
      svd.compute(Phi);
      if (lambda < 0) {
        double alpha = config.get<double>("ALPHA");
        double beta = config.get<double>("BETA");
        EA ea(config, svd, T);
        ea.run(alpha, beta);
        lambda = alpha / beta;
        config.remove("ALPHA");
        config.remove("BETA");
        config.add("ALPHA", alpha);
        config.add("BETA", beta);
      }

      Sigma = get_sigma(svd, lambda);
      config_add_sigma(config, Sigma);

      if (verbose)
        std::cout << std::endl << "REG LAMBDA: " << lambda << std::endl;

      ols.solveWithPseudoinverse(Phi.rows(), Phi.cols(), T, weights,
                                 svd.getSVDWorkspace(), rcond, lambda);
    } else if (oalgo == 4) {
      ols.solve(Phi, T, weights, rcond,
                tadah::models::Algorithm::NormalEquations);
    } else {
      throw std::runtime_error("Unsupported OALGO: " + std::to_string(oalgo));
    }

  }
  static void train(tadah::core::Config &config, tadah::core::phi_type &Phi, tadah::core::t_type &T, tadah::core::t_type &weights,
                    tadah::core::Matrix &Sigma) {
    tadah::models::memory::ModelsWorkspaceManager workspaceManager;
    train(config, Phi, T, weights, Sigma, workspaceManager);
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
   * @param config Configuration object.
   * @param Sigma Covariance matrix to be added.
   */
  static void config_add_sigma(tadah::core::Config &config, tadah::core::Matrix &Sigma) {
    tadah::core::t_type Sigma_as_vector(Sigma.data(), Sigma.cols() * Sigma.rows());
    config.add("SIGMA", Sigma.rows());
    config.add<tadah::core::t_type>("SIGMA", Sigma_as_vector);
  }

  /**
   * @brief Read and reconstruct the covariance matrix from the configuration.
   *
   * Retrieves the Sigma matrix stored in the configuration.
   *
   * @param config Configuration object.
   * @param Sigma tadah::core::Matrix to store the reconstructed Sigma.
   * @throws std::runtime_error if Sigma is not computed.
   */
  static void read_sigma(tadah::core::Config &config, tadah::core::Matrix &Sigma) {
    using V = std::vector<double>;
    size_t N;
    try {
      N = config.get<size_t>("SIGMA");
    } catch (const std::runtime_error &error) {
      throw std::runtime_error("\nSigma matrix is not computed.\nHint: It is "
                               "only computed for LAMBDA != 0.\n");
    }

    V Sigma_as_vector(N * N + 1);
    config.get<V>("SIGMA", Sigma_as_vector);
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
