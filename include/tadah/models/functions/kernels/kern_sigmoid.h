#ifndef KERN_SIGMOID_H
#define KERN_SIGMOID_H

#include <tadah/models/functions/kernels/kern_base.h>

namespace tadah {
namespace models {

/** \brief Sigmoid kernel.
 *
 * Defined for two vectors **x** and **y**:
 *
 * \f[
 * K(\mathbf{x}, \mathbf{y}) = \tanh\big( \gamma*\mathbf{x}^T \mathbf{y}
 * + \mathrm{c} \big)
 * \f]
 *
 *  Required tadah::core::Context key: \ref MPARAMS <double> \f$ \gamma \f$  c
 *
 *  @see Kern_Base
 */
class Kern_Sigmoid : public virtual Kern_Base {
    public:
        double gamma;   // scalling
        double c;   // shift
        /**
         * Label used for this class
         */
        Kern_Sigmoid ();
        Kern_Sigmoid (const tadah::core::Context &c);
        std::string label = "Kern_Sigmoid";
        double operator() (const tadah::core::aed_type& b, const tadah::core::aed_type& af) const override;
        tadah::core::aed_type derivative(const tadah::core::aed_type& b, const tadah::core::aed_type& af) const override;
        double prime(const tadah::core::aed_type& b, const tadah::core::aed_type& af, const tadah::core::aed_type& ff) const override;

        std::string get_label() const override;

};
} // namespace models
} // namespace tadah
#endif
