#ifndef M_BLR_TRAIN_H
#define M_BLR_TRAIN_H

#include "tadah/core/core_types.h"
#include "tadah/core/utils.h"
#include <tadah/models/m_train.h>
#include <tadah/models/m_blr_core.h>
#include <tadah/models/linear_regressor.h>
#include <tadah/models/functions/function_base.h>
#include <tadah/models/functions/basis_functions/bf_base.h>
#include <tadah/core/context.h>
#include <tadah/models/memory/imodels_workspace_manager.h>

#include <vector>

namespace tadah {
namespace models {

/**
 * @class M_BLR_Train
 * @brief Implements training for Basis Linear Regression (BLR).
 *
 * Extends BLR core functionalities with training capabilities, using a specified basis function type.
 *
 * @tparam BF Basis function type, default is Function_Base.
 */
template <class BF = Function_Base&>
class M_BLR_Train : public M_BLR_Core<BF>, public M_Train {
public:
    using M_BLR_Core<BF>::trained;
    using M_BLR_Core<BF>::context;
    using M_BLR_Core<BF>::weights;
    using M_BLR_Core<BF>::Sigma;
    using M_BLR_Core<BF>::bf;

    /**
     * @brief Constructor using a configuration object.
     *
     * Initializes BLR core with the given configuration.
     *
     * @param c Configuration object reference.
     */
    M_BLR_Train(tadah::core::Context &c) : M_BLR_Core<BF>(c) {}

    /**
     * @brief Constructor using a basis function and configuration object.
     *
     * Initializes BLR core with the given basis function and configuration.
     *
     * @param bf Basis function object.
     * @param c Configuration object reference.
     */
    M_BLR_Train(BF &bf, tadah::core::Context &c) : M_BLR_Core<BF>(bf, c) {}
    M_BLR_Train(BF &bf, tadah::core::Context &c, tadah::models::memory::IModelsWorkspaceManager &workspaceManager) :
    M_BLR_Core<BF>(bf, c),
    M_Train(workspaceManager)
  {}

    /**
     * @brief Trains the model on given feature and target matrices, whose
     *        observations are rows - cols (a design matrix with its
     *        ridge-augmentation rows and no rows of weight 0).
     */
  void train(tadah::core::phi_type &Phi, tadah::core::t_type &T) {
    train(Phi, T, LinearRegressor::augmented_observations(Phi));
  }

    /**
     * @brief Trains the model on given feature and target matrices.
     *
     * Uses a linear regressor to perform training. The model holds the
     * weight covariance of this fit (M_Core::Sigma, noise_variance) or
     * none: OALGO 3 computes it, of every weight (a pinned one with a zero
     * row and column; with every weight pinned all zeros and the noise
     * RSS / n_obs of the pinned weights), any other OALGO none. The
     * covariance of an earlier fit is removed first, whatever this fit
     * computes. Nothing is written into the Context.
     *
     * @param Phi Feature matrix.
     * @param T Target vector.
     * @param n_obs The observations: the rows of Phi with a non-zero
     *        weight, without the ridge-augmentation rows.
     * @throws std::runtime_error if the object is already trained.
     */
  void train(tadah::core::phi_type &Phi, tadah::core::t_type &T, size_t n_obs) {
    if (trained) {
      throw std::runtime_error("This object is already trained!");
    }
    this->clear_covariance("the fit computed no weight covariance");

    const bool has_fixweight = context.contains("FIXWEIGHT");
    const bool has_fixindex = context.contains("FIXINDEX");
    if (has_fixweight != has_fixindex) {
      throw std::runtime_error(
          "FIXINDEX requires FIXWEIGHT (and vice versa). Provide both keys "
          "with matching lengths; the i-th FIXWEIGHT value pins the weight "
          "named by the i-th FIXINDEX entry. There is no auto-seeding from "
          "WEIGHTS.");
    }

    if (has_fixweight && has_fixindex) {
      std::vector<std::string> indices_str(context.size("FIXINDEX"));
      context.get("FIXINDEX", indices_str);

      // Parse indices and perform checks
      std::vector<size_t> indices = tadah::core::parse_indices(indices_str);

      // Check that the min index is >= 1 and max index is <= Phi.cols()
      if (!indices.empty()) {
        size_t max_index = *std::max_element(indices.begin(), indices.end());
        if (*std::min_element(indices.begin(), indices.end()) < 1 || max_index > Phi.cols()) {
          throw std::runtime_error("FIXINDEX: Indices out of bounds: valid range is from 1 to " + std::to_string(Phi.cols()));
        }
      }

      // Adjust for 0-based indexing
      for (auto &i : indices) i--;

      tadah::core::t_type w_f(context.size("FIXWEIGHT"));
      if (w_f.size() != indices.size()) {
        throw std::runtime_error("FIXWEIGHT and FIXINDEX differ in size: " + std::to_string(w_f.size()) + " and " + std::to_string(indices.size()));
      }

      context.get("FIXWEIGHT", w_f);

      tadah::core::t_type T_r;
      auto move_map = prep_train_with_residuals(Phi, T, indices, w_f, T_r);

      // Resize Phi for training, the rightmost columns correspond to FIXINDEX
      Phi.resize(Phi.rows(), Phi.cols() - indices.size());


      CovarianceFit fit;
      if (Phi.cols() > 0) {
        // The function's intercept columns among the free ones, at their
        // new positions: the ridge leaves those out and penalises every
        // other free column. n_obs is that of the full design matrix: the
        // pinned columns leave Phi, its rows stay.
        fit = LinearRegressor::train(context, Phi, T_r, weights, Sigma, *workspaceManager_,
                               free_intercept_columns(bf.intercept_columns(context),
                                                      move_map, Phi.cols()),
                               n_obs);
      } else {
        // Every column is fixed (FIXINDEX covers all weights): there is nothing
        // left to solve, so skip the regressor entirely and let the fixed
        // values below become the model weights. This is the path exercised by
        // `tadah refit`'s coefficient-search HPO target (OPTIM FIXWEIGHT over
        // all weights). Calling LinearRegressor::train with 0 columns would
        // dereference a zero-size workspace.
        weights = tadah::core::t_type();
      }

      // Reassemble the full weight vector. {move_map[0..n_free-1]} and
      // indices partition {0..N-1}, so every slot is written exactly once:
      // the solved coefficients go back to the original positions of the
      // free columns, and each fixed value lands on its own weight index
      // (indices order == FIXWEIGHT declaration order). The all-fixed case
      // (n_free == 0) reduces to the pure scatter.
      const size_t n_free = weights.size();
      tadah::core::t_type w_full(n_free + indices.size());
      for (size_t j = 0; j < n_free; ++j) {
        w_full[move_map[j]] = weights[j];
      }
      for (size_t k = 0; k < indices.size(); ++k) {
        w_full[indices[k]] = w_f[k];
      }
      weights = w_full;
      // OALGO 3 computed the covariance of the free weights; the model's
      // is that of every weight, a pinned one with none.
      if (n_free > 0) {
        if (fit.computed) {
          sigma_of_all_weights(Sigma, move_map, n_free, weights.size());
          this->set_covariance(Sigma, fit.noise_variance);
        } else {
          this->clear_covariance(fit.why_not);
        }
      } else if (context.template get<int>("OALGO", 0) == 3) {
        double s2 = 0.0;
        std::string why;
        if (pinned_covariance(T_r, weights.size(), n_obs, Sigma, s2, why))
          this->set_covariance(Sigma, s2);
        else
          this->clear_covariance(why);
      } else {
        this->clear_covariance(fit.why_not.empty()
            ? "OALGO " + std::to_string(context.template get<int>("OALGO", 0)) +
              " computes no weight covariance (only OALGO 3, the SVD solve, does)"
            : fit.why_not);
      }
      trained = true;
    }
    else {
      const CovarianceFit fit =
          LinearRegressor::train(context, Phi, T, weights, Sigma, *workspaceManager_,
                                 bf.intercept_columns(context), n_obs);
      if (fit.computed) this->set_covariance(Sigma, fit.noise_variance);
      else this->clear_covariance(fit.why_not);
      trained = true;
    }

  }
};

} // namespace models
} // namespace tadah
#endif // M_BLR_TRAIN_H

