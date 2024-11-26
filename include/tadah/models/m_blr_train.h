#ifndef M_BLR_TRAIN_H
#define M_BLR_TRAIN_H

#include <tadah/models/m_train.h>
#include <tadah/models/m_blr_core.h>
#include <tadah/models/linear_regressor.h>
#include <tadah/models/functions/function_base.h>
#include <tadah/models/functions/basis_functions/bf_base.h>
#include <tadah/core/config.h>

#include <iostream>

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
    M_BLR_Train(Config &c) : M_BLR_Core<BF>(c) {}

    /**
     * @brief Constructor using a basis function and configuration object.
     *
     * Initializes BLR core with the given basis function and configuration.
     *
     * @param bf Basis function object.
     * @param c Configuration object reference.
     */
    M_BLR_Train(BF &bf, Config &c) : M_BLR_Core<BF>(bf, c) {}

    /**
     * @brief Trains the model on given feature and target matrices.
     *
     * Uses a linear regressor to perform training.
     *
     * @param Phi Feature matrix.
     * @param T Target vector.
     * @throws std::runtime_error if the object is already trained.
     */
    void train(phi_type &Phi, t_type &T) {
        if (trained) {
            throw std::runtime_error("This object is already trained!");
        }
        LinearRegressor::train(config, Phi, T, weights, Sigma);
        trained = true;
    }
};

#endif // M_BLR_TRAIN_H

