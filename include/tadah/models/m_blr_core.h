#ifndef M_BLR_CORE_H
#define M_BLR_CORE_H

#include <tadah/models/m_core.h>
#include <tadah/models/functions/function_base.h>
#include <tadah/models/functions/basis_functions/bf_base.h>
#include <tadah/core/config.h>

#include <iostream>

/**
 * @class M_BLR_Core
 * @brief Implements core functionality for Basis Linear Regression (BLR).
 *
 * This template class supports BLR using a specified basis function type.
 *
 * @tparam BF Basis function type, default is Function_Base.
 */
template <class BF = Function_Base&>
class M_BLR_Core : public virtual M_Core {
public:
    Matrix Sigma; ///< Covariance matrix.
    Config &config; ///< Configuration reference.
    BF bf; ///< Basis function.

    /**
     * @brief Constructor using a configuration object.
     *
     * Initializes the basis function with the provided configuration and checks its type.
     *
     * @param c Configuration object reference.
     */
    M_BLR_Core(Config &c) : config(c), bf(c) {
        static_assert(std::is_base_of<BF_Base, BF>::value,
                      "\nThe provided Basis function is not of a BasisFunction type.\n"
                      "For BLR use basis functions defined in basis_functions.h.\n"
                      "If you want to use this kernel functions use M_KRR class instead.\n"
                      "e.g. Use BF_Linear instead of Kern_Linear.\n");
        init();
    }

    /**
     * @brief Constructor using a basis function and configuration object.
     *
     * Verifies the basis function is of the correct type and initializes.
     *
     * @param bf Basis function object.
     * @param c Configuration object reference.
     */
    M_BLR_Core(BF &bf, Config &c) : config(c), bf(bf) {
        if (dynamic_cast<BF_Base*>(&bf) == nullptr)
            throw std::invalid_argument("Provided object is not of BF_Base type");
        init();
    }

    /**
     * @brief Predicts output using the basis function.
     *
     * @param v Input vector for prediction.
     * @return Predicted value.
     */
    double predict(const aed_type2 &v) const {
        return bf.epredict(get_weights(), v);
    }

    /**
     * @brief Retrieves the uncertainty of weights.
     *
     * Throws an error if LAMBDA is not zero.
     *
     * @return Vector of uncertainties based on the Sigma diagonal.
     * @throws std::runtime_error if LAMBDA is not zero.
     */
    t_type get_weights_uncertainty() const {
        double lambda = config.template get<double>("LAMBDA");
        if (lambda != 0)
            throw std::runtime_error("Sigma matrix is only computed for LAMBDA != 0");
        return Sigma.diagonal();
    }

private:
    /**
     * @brief Initializes the model based on the configuration.
     *
     * Loads weights if available in the configuration and sets the verbose level.
     */
    void init() {
        if (config.exist("WEIGHTS")) {
            weights.resize(config.size("WEIGHTS"));
            config.get<t_type>("WEIGHTS", weights);
            trained = true;
        }
        verbose = config.get<int>("VERBOSE");
    }
};

#endif // M_BLR_CORE_H
