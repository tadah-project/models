#ifndef M_CORE_H
#define M_CORE_H

#include <tadah/core/core_types.h>

#include <cmath>
#include <stdexcept>
#include <string>

namespace tadah {
namespace models {

/**
 * @class M_Core
 * @brief Abstract base class for model cores.
 *
 * Provides common functionality and interface for model implementations, including training status
 * and weight management.
 */
class M_Core {
public:
    int verbose = 0; ///< Verbose level for logging.
    bool trained = false; ///< Indicates if the model has been trained.
    tadah::core::t_type weights; ///< Weights vector for the model.

    /** The covariance of the weights the last fit computed (OALGO 3,
     *  LinearRegressor::train): N x N, of `weights` in their own
     *  parameterisation, s^2 (A^T A + lambda I)^-1 over the directions the
     *  solve keeps; a pinned weight (FIXINDEX) has a zero row and column.
     *  0 x 0 when the last fit computed none (has_covariance() false).
     *  The potential gets it as WCOV (MLIP get_param_file); nothing reads
     *  or writes it in the Context. */
    tadah::core::Matrix Sigma;
    /** The noise variance s^2 of the last fit per unit row weight,
     *  RSS / (n_obs - gamma) (RSS / n_obs with every weight pinned); 0
     *  when it computed no covariance. */
    double noise_variance = 0.0;

    /** Whether the last fit computed the covariance Sigma. */
    bool has_covariance() const { return covariance_computed_; }

    /** Why the last fit computed no covariance (empty when it did). */
    const std::string &covariance_unavailable() const {
      return covariance_unavailable_;
    }

    /** The covariance @p S and noise variance @p s2 of the fit that just
     *  ran. */
    void set_covariance(const tadah::core::Matrix &S, double s2) {
      Sigma = S;
      noise_variance = s2;
      covariance_computed_ = true;
      covariance_unavailable_.clear();
    }

    /** The fit that just ran (or starts) computes no covariance, for the
     *  reason @p why; an earlier fit's is removed. */
    void clear_covariance(const std::string &why) {
      Sigma = tadah::core::Matrix();
      noise_variance = 0.0;
      covariance_computed_ = false;
      covariance_unavailable_ = why;
    }

  int is_verbose() {return verbose;}

    /**
     * @brief Virtual destructor for polymorphic deletion.
     */
    virtual ~M_Core() {}

    /**
     * @brief Checks if the model has been trained.
     *
     * @return True if the model is trained, otherwise false.
     */
    bool is_trained() const {
        return trained;
    }

    /** Reset trained state so the model can be re-trained in-place. */
    void reset_trained() {
        trained = false;
    }

    /**
     * @brief Retrieves the weights of the model.
     *
     * @return Constant reference to the weights vector.
     */
    const tadah::core::t_type &get_weights() const {
        return weights;
    }

    /**
     * @brief Sets the model weights.
     *
     * @param w New weights vector to be set.
     */
    void set_weights(const tadah::core::t_type w) {
        weights = w;
    }

    /**
     * @brief Pure virtual function for making predictions.
     *
     * Must be implemented by derived classes.
     *
     * @param v Input vector for prediction.
     * @return Predicted value.
     */
    virtual double predict(const tadah::core::aed_type &v) const = 0;

    /**
     * @brief The standard deviations of the weights, sqrt(diag Sigma), of
     *        the last fit (OALGO 3): one per weight, in its units, 0 for a
     *        pinned weight.
     *
     * @throws std::runtime_error when the last fit computed no covariance,
     *         naming why (another OALGO than 3, too few observations, no
     *         fit yet).
     */
    virtual tadah::core::t_type get_weights_uncertainty() const {
      if (!covariance_computed_)
        throw std::runtime_error(
            "The weight uncertainties are not available: " +
            (covariance_unavailable_.empty()
                 ? std::string("the model has not been fitted; only an "
                               "OALGO 3 fit computes the weight covariance")
                 : covariance_unavailable_) + ".");
      tadah::core::t_type sd(Sigma.rows());
      for (size_t i = 0; i < Sigma.rows(); ++i) {
        const double v = Sigma(i, i);
        if (!(v >= 0.0))
          throw std::logic_error(
              "M_Core::get_weights_uncertainty: the variance of weight " +
              std::to_string(i + 1) + " is " + std::to_string(v) +
              ", not a variance");
        sd(i) = std::sqrt(v);
      }
      return sd;
    }

private:
    bool covariance_computed_ = false;
    std::string covariance_unavailable_;
};

} // namespace models
} // namespace tadah
#endif // M_CORE_H

