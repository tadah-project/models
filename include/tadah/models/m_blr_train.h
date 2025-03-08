#ifndef M_BLR_TRAIN_H
#define M_BLR_TRAIN_H

#include "tadah/core/core_types.h"
#include "tadah/core/utils.h"
#include <tadah/models/m_train.h>
#include <tadah/models/m_blr_core.h>
#include <tadah/models/linear_regressor.h>
#include <tadah/models/functions/function_base.h>
#include <tadah/models/functions/basis_functions/bf_base.h>
#include <tadah/core/config.h>
#include <tadah/models/memory/IModelsWorkspaceManager.h>

#include <iostream>
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
    using M_BLR_Core<BF>::config;
    using M_BLR_Core<BF>::weights;
    using M_BLR_Core<BF>::Sigma;

    /**
     * @brief Constructor using a configuration object.
     *
     * Initializes BLR core with the given configuration.
     *
     * @param c Configuration object reference.
     */
    M_BLR_Train(tadah::core::Config &c) : M_BLR_Core<BF>(c) {}

    /**
     * @brief Constructor using a basis function and configuration object.
     *
     * Initializes BLR core with the given basis function and configuration.
     *
     * @param bf Basis function object.
     * @param c Configuration object reference.
     */
    M_BLR_Train(BF &bf, tadah::core::Config &c) : M_BLR_Core<BF>(bf, c) {}
    M_BLR_Train(BF &bf, tadah::core::Config &c, tadah::models::memory::IModelsWorkspaceManager &workspaceManager) :
    M_BLR_Core<BF>(bf, c),
    M_Train(workspaceManager)
  {}

    /**
     * @brief Trains the model on given feature and target matrices.
     *
     * Uses a linear regressor to perform training.
     *
     * @param Phi Feature matrix.
     * @param T Target vector.
     * @throws std::runtime_error if the object is already trained.
     */
  void train(tadah::core::phi_type &Phi, tadah::core::t_type &T) {
    if (trained) {
      throw std::runtime_error("This object is already trained!");
    }

    if (config.exist("FIXWEIGHT") && config.exist("FIXINDEX")) {
      std::vector<std::string> indices_str(config.size("FIXINDEX"));
      config.get("FIXINDEX", indices_str);

      // Parse indices and perform checks
      std::vector<size_t> indices = parse_indices(indices_str);

      // Check that the min index is >= 1 and max index is <= Phi.cols()
      if (!indices.empty()) {
        size_t max_index = *std::max_element(indices.begin(), indices.end());
        if (*std::min_element(indices.begin(), indices.end()) < 1 || max_index > Phi.cols()) {
          throw std::runtime_error("FIXINDEX: Indices out of bounds: valid range is from 1 to " + std::to_string(Phi.cols()));
        }
      }

      // Adjust for 0-based indexing
      for (auto &i : indices) i--;

      tadah::core::t_type w_f(config.size("FIXWEIGHT"));
      if (w_f.size() != indices.size()) {
        throw std::runtime_error("FIXWEIGHT and FIXINDEX differ in size: " + std::to_string(w_f.size()) + " and " + std::to_string(indices.size()));
      }

      config.get("FIXWEIGHT", w_f);

      tadah::core::t_type T_r;
      auto move_map = prep_train_with_residuals(Phi, T, indices, w_f, T_r);

      // Resize Phi for training
      Phi.resize(Phi.rows(), Phi.cols() - indices.size());

      LinearRegressor::train(config, Phi, T_r, weights, Sigma,*workspaceManager_);

      tadah::core::t_type w_temp(weights.size() + indices.size());

      for (size_t i = 0; i < w_f.size(); ++i) {
        w_temp[w_temp.size() - w_f.size() + i] = w_f[i];
      }

      for (size_t i = 0; i < weights.size(); ++i) {
        w_temp[i] = weights[i];
      }

      weights = w_temp;
      reverse_vector(weights, move_map);
      trained = true;
    }
    else {
      LinearRegressor::train(config, Phi, T, weights, Sigma,*workspaceManager_);
      trained = true;
    }

  }
};

} // namespace models
} // namespace tadah
#endif // M_BLR_TRAIN_H

