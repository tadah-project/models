#ifndef KERN_RBF_H
#define KERN_RBF_H

#include <tadah/models/functions/kernels/kern_base.h>

namespace tadah {
namespace models {

/** \brief Radial Basis Function kernel
 *
 * Defined for two vectors **x** and **y**:
 *
 * \f[
 * K(\mathbf{x}, \mathbf{y}) = \exp\Big( -\frac{||\mathbf{x}-\mathbf{y}||^2}
 * {2\sigma^2} \Big)
 * \f]
 *
 *  Required tadah::core::Context key: \ref MPARAMS <double> \f$ \sigma \f$
 *
 *  @see Kern_Base Kern_Linear Kern_Quadratic
 */
class Kern_RBF : public virtual Kern_Base {
    public:
        double sigma;
        /**
         * Label used for this class
         */
        Kern_RBF ();
        Kern_RBF (const tadah::core::Context &c);
        std::string label = "Kern_RBF";
        double operator() (const tadah::core::aed_type& b, const tadah::core::aed_type& af) const override;
        tadah::core::aed_type derivative(const tadah::core::aed_type& b, const tadah::core::aed_type& af) const override;
        double prime(const tadah::core::aed_type& b, const tadah::core::aed_type& af, const tadah::core::aed_type& ff) const override;
        std::string get_label() const override;

    private:

        double gamma;

};
} // namespace models
} // namespace tadah
#endif
