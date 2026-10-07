#ifndef M_KRR_Core_H
#define M_KRR_Core_H

#include <tadah/core/context.h>
#include <tadah/models/functions/function_base.h>
#include <tadah/models/functions/kernels/kern_base.h>
#include <tadah/models/m_core.h>

#include <iostream>
#include <stdexcept>
#include <string>

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
    tadah::core::Matrix Sigma; ///< Covariance matrix.
    tadah::core::Context &context; ///< Configuration reference.
    K kernel; ///< Kernel function.

    /**
     * @brief Constructor using a configuration object.
     * 
     * Initializes the kernel with the provided configuration and checks its type.
     *
     * @param c Configuration object reference.
     */
    M_KRR_Core(tadah::core::Context &c) :
        context(c),
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
    M_KRR_Core(K &kernel, tadah::core::Context &c) :
        context(c),
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
    double predict(const tadah::core::aed_type &v) const {
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
    tadah::core::t_type get_weights_uncertainty() const {
        double lambda = context.template get<double>("LAMBDA");
        if(lambda == 0) 
            throw std::runtime_error("Sigma matrix is only computed for LAMBDA != 0");
        return Sigma.diagonal();
    }

private:
    /**
     * @brief Initializes the model based on the configuration.
     * 
     * Loads weights if available in the configuration and sets the verbose
     * level. Loaded weights do not make the model trained (is_trained()):
     * train() fits it, replacing them.
     *
     * A non-linear kernel (any kernel but Kern_Linear) evaluates its
     * weights against its basis vectors, which it read from SBASIS and
     * BASIS (Kern_Base::read_basis_from_context): with WEIGHTS, the basis
     * must be there and hold one vector per weight, else this throws
     * naming the keys (without a basis the model predicted E = 0, and with
     * fewer vectors than weights it read past the basis).
     */
    void init() {
        if (context.contains("WEIGHTS")) {
            weights.resize(context.size("WEIGHTS"));
            context.get<tadah::core::t_type>("WEIGHTS", weights);
            // Loaded weights predict, but the model is not trained: train()
            // fits it (a potential used as the configuration), as M_BLR.
            check_basis_of_weights();
        }
        verbose = context.get<int>("VERBOSE");

    }

    void check_basis_of_weights() const {
        if (kernel.get_label() == "Kern_Linear") return;
        const Kern_Base *kb = dynamic_cast<const Kern_Base*>(&kernel);
        if (kb == nullptr) return;   // the constructors refuse a non-kernel
        const size_t nb = kb->get_basis().cols();
        if (nb == 0)
            throw std::runtime_error(
                "M_KRR " + kernel.get_label() + ": the potential has WEIGHTS "
                "but no basis vectors (SBASIS and BASIS): a non-linear kernel "
                "evaluates its weights against the basis vectors that "
                "training wrote with them.");
        if (nb != weights.size())
            throw std::runtime_error(
                "M_KRR " + kernel.get_label() + ": the potential has " +
                std::to_string(weights.size()) + " WEIGHTS but SBASIS " +
                std::to_string(nb) + " basis vectors (BASIS): a non-linear "
                "kernel has one weight per basis vector.");
    }
};

} // namespace models
} // namespace tadah
#endif // M_KRR_Core_H
