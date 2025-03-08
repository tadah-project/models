#ifndef KERN_POLYNOMIAL_H
#define KERN_POLYNOMIAL_H

#include <tadah/models/functions/kernels/kern_base.h>

namespace tadah {
namespace models {

/** \brief Polynomial kernel.
 *
 * Defined for two vectors **x** and **y**:
 *
 * \f[
 * K(\mathbf{x}, \mathbf{y}) = \big( \gamma*\mathbf{x}^T \mathbf{y} +
 * \mathrm{c} \big)^{\mathrm{d}}
 * \f]
 *
 *  Required tadah::core::Config key: \ref MPARAMS <int> d  <double> \f$ \gamma \f$  c
 *
 *  @see Kern_Base
 */
class Kern_Polynomial : public virtual Kern_Base {
    public:
        Kern_Polynomial ();
        Kern_Polynomial (const tadah::core::Config &c);
        int d; // power
        double gamma;   // scalling
        double c;   // shift
        /**
         * Label used for this class
         */
        std::string label = "Kern_Polynomial";
        std::string get_label() const override;
        double operator() (const tadah::core::aed_type& b, const tadah::core::aed_type& af) const override;
        tadah::core::aed_type derivative(const tadah::core::aed_type& b, const tadah::core::aed_type& af) const override;

        double prime(const tadah::core::aed_type& b, const tadah::core::aed_type& af, const tadah::core::aed_type& ff) const override;

};
} // namespace models
} // namespace tadah
#endif
