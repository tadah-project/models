#ifndef M_KRR_TRAIN_H
#define M_KRR_TRAIN_H

#include <tadah/models/linear_regressor.h>
#include <tadah/models/functions/function_base.h>
#include <tadah/models/functions/kernels/kern_base.h>
#include <tadah/models/m_train.h>
#include <tadah/models/m_krr_core.h>
#include <tadah/models/ekm.h>
#include <tadah/core/maths.h>
#include <tadah/core/config.h>
#include <tadah/models/memory/IModelsWorkspaceManager.h>

#include <iostream>
#include <limits>

namespace tadah {
namespace models {

/**
 * @class M_KRR_Train
 * @brief Implements Kernel Ridge Regression (KRR) training.
 *
 * Extends KRR core and train functionalities with specific implementations.
 *
 * @tparam Kern Kernel type, default is Function_Base.
 */
template <class Kern = Function_Base&>
class M_KRR_Train : public M_KRR_Core<Kern>, public M_Train {
public:
    EKM<Kern> ekm; ///< Extended Kernel Matrix object.

    using M_KRR_Core<Kern>::trained;
    using M_KRR_Core<Kern>::config;
    using M_KRR_Core<Kern>::kernel;
    using M_KRR_Core<Kern>::weights;
    using M_KRR_Core<Kern>::Sigma;

    /**
     * @brief Constructor using a configuration object.
     *
     * Initializes KRR core and EKM with the given configuration.
     *
     * @param c Configuration object reference.
     */
    M_KRR_Train(Config &c) : M_KRR_Core<Kern>(c), ekm(c) {}

    /**
     * @brief Constructor using a kernel and configuration object.
     *
     * Initializes KRR core and EKM with the given kernel and configuration.
     *
     * @param kernel Kernel object.
     * @param c Configuration object reference.
     */
    M_KRR_Train(Kern &kernel, Config &c) : M_KRR_Core<Kern>(kernel, c), ekm(kernel) {}
    M_KRR_Train(Kern &kernel, Config &c, tadah::models::memory::IModelsWorkspaceManager &workspaceManager) :
    M_KRR_Core<Kern>(kernel, c),
    M_Train(workspaceManager),
    ekm(kernel)
  {}

    /**
     * @brief Trains the model on given feature and target matrices.
     *
     * Uses EKM and linear regressor to perform training.
     *
     * @param Phi Feature matrix.
     * @param T Target vector.
     * @throws std::runtime_error if the object is already trained.
     */
  void train(phi_type &Phi, t_type &T) {
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

      t_type w_f(config.size("FIXWEIGHT"));
      if (w_f.size() != indices.size()) {
        throw std::runtime_error("FIXWEIGHT and FIXINDEX differ in size: " + std::to_string(w_f.size()) + " and " + std::to_string(indices.size()));
      }

      config.get("FIXWEIGHT", w_f);

      t_type T_r;
      auto move_map = prep_train_with_residuals(Phi, T, indices, w_f, T_r);

      reverse_columns(Phi, move_map);

      for (const auto &i : indices) {
        for (size_t j=0; j<Phi.rows(); ++j) {
          Phi(j,i)=0;
        }
      }

      if (kernel.get_label() != "Kern_Linear") {
        ekm.project(Phi);
      }
      LinearRegressor::train(config, Phi, T_r, weights, Sigma, *workspaceManager_);
      if (kernel.get_label() != "Kern_Linear") {
        weights = ekm.EKM_mat * weights;
      }

      // reverse_vector(weights, move_map);
      for (size_t i=0; i<indices.size(); ++i) {
        weights[indices[i]]=w_f[i];
      }
      trained = true;
    }
    else {
      if (kernel.get_label() != "Kern_Linear") {
        ekm.project(Phi);
      }
      LinearRegressor::train(config, Phi, T, weights, Sigma, *workspaceManager_);
      if (kernel.get_label() != "Kern_Linear") {
        weights = ekm.EKM_mat * weights;
      }
      trained = true;
    }

  }

    /**
     * @brief Standard KRR training using covariance matrix computation.
     *
     * Trains the model on energies using the covariance matrix K.
     *
     * @param B Basis matrix.
     * @param T Target vector.
     * @throws std::runtime_error if lambda is negative for evidence approximation.
     */
    void train2(Matrix &B, t_type &T) {
        double lambda = config.template get<double>("LAMBDA");

        if (lambda < 0)
            throw std::runtime_error("This KRR implementation cannot use LAMBDA -1 (Evidence Approximation)");
        
        if (M_KRR_Core<Kern>::is_verbose()) std::cout << "Computing Kernel Matrix" << std::endl;
        Matrix K(B.cols(), B.cols());

        for (size_t i = 0; i < K.cols(); i++) {
            for (size_t j = 0; j < K.cols(); j++) {
                K(i, j) = kernel(B.col(i), B.col(j));
                if (i == j) K(i, j) += lambda;
            }
        }
        if (M_KRR_Core<Kern>::is_verbose()) std::cout << "Matrix condition number: " << condition_number(K) << std::endl;

        if (M_KRR_Core<Kern>::is_verbose()) std::cout << "Solving..." << std::flush;
        //double rcond = config.template size("LAMBDA")==2 ? config.template get<double>("LAMBDA",1) : 1e-8;
        weights = solve_posv(K, T);
        if (M_KRR_Core<Kern>::is_verbose()) std::cout << "Done" << std::endl;


        trained = true;
} // namespace models
} // namespace tadah

    }
};

#endif // M_KRR_TRAIN_H
