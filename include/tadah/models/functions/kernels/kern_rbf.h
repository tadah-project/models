#ifndef KERN_RBF_H
#define KERN_RBF_H

#include <tadah/models/functions/kernels/kern_base.h>

/** \brief Radial Basis Function kernel
 *
 * Defined for two vectors **x** and **y**:
 *
 * \f[
 * K(\mathbf{x}, \mathbf{y}) = \exp\Big( -\frac{||\mathbf{x}-\mathbf{y}||^2}
 * {2\sigma^2} \Big)
 * \f]
 *
 *  Required Config key: \ref MPARAMS <double> \f$ \sigma \f$
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
        Kern_RBF (const Config &c);
        std::string label = "Kern_RBF";
        double operator() (const aed_type& b, const aed_type& af) const override;
        aed_type derivative(const aed_type& b, const aed_type& af) const override;
        double prime(const aed_type& b, const aed_type& af, const aed_type& ff) const override;
        std::string get_label() const override;

    private:
        double gamma;

};
#endif
