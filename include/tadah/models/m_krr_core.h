#ifndef M_KRR_Core_H
#define M_KRR_Core_H

#include <tadah/core/config.h>
#include <tadah/models/functions/function_base.h>
#include <tadah/models/functions/kernels/kern_base.h>
#include <tadah/models/m_core.h>

#include <iostream>

namespace tadah {
namespace models {

/**
 * @class M_KRR_Core
 * @brief Implements core functionality for Kernel Ridge Regression (KRR).
 *
 * This template class supports Kernel Ridge Regression using a specified kernel type.
 * 
 * @tparam K Kernel type, default is Function_Base.
 */
template <class K = Function_Base&>
class M_KRR_Core : public virtual M_Core {
public:
    Matrix Sigma; ///< Covariance matrix.
    Config &config; ///< Configuration reference.
    K kernel; ///< Kernel function.

    /**
     * @brief Constructor using a configuration object.
     * 
     * Initializes the kernel with the provided configuration and checks its type.
     *
     * @param c Configuration object reference.
     */
    M_KRR_Core(Config &c) :
        config(c),
        kernel(c)
    {
        static_assert(std::is_base_of<Kern_Base, K>::value,
                      "\nThe provided Kernel K is not of a Kernel type.\n"
                      "For KRR use kernels defined in kernel.h.\n"
                      "e.g. Use Kern_Linear instead of BF_Linear.\n");
        init();
    }

    /**
     * @brief Constructor using a kernel and configuration object.
     * 
     * Verifies the kernel is of the correct type and initializes.
     *
     * @param kernel Kernel object.
     * @param c Configuration object reference.
     */
    M_KRR_Core(K &kernel, Config &c) :
        config(c),
        kernel(kernel)
    {
        if (dynamic_cast<Kern_Base*>(&kernel) == nullptr)
            throw std::invalid_argument("Provided object is not of Kern_Base type");
        init();
    }

    /**
     * @brief Predicts output using the kernel function.
     * 
     * @param v Input vector for prediction.
     * @return Predicted value.
     */
    double predict(const aed_type &v) const {
        return kernel.epredict(weights, v);
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
        if(lambda != 0) 
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
} // namespace models
} // namespace tadah

    }
};

#endif // M_KRR_Core_H
