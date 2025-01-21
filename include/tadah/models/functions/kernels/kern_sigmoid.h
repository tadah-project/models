#ifndef KERN_SIGMOID_H
#define KERN_SIGMOID_H

#include <tadah/models/functions/kernels/kern_base.h>

/** \brief Sigmoid kernel.
 *
 * Defined for two vectors **x** and **y**:
 *
 * \f[
 * K(\mathbf{x}, \mathbf{y}) = \tanh\big( \gamma*\mathbf{x}^T \mathbf{y}
 * + \mathrm{c} \big)
 * \f]
 *
 *  Required Config key: \ref MPARAMS <double> \f$ \gamma \f$  c
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
        Kern_Sigmoid (const Config &c);
        std::string label = "Kern_Sigmoid";
        double operator() (const aed_type& b, const aed_type& af) const override;
        aed_type derivative(const aed_type& b, const aed_type& af) const override;
        double prime(const aed_type& b, const aed_type& af, const aed_type& ff) const override;
        std::string get_label() const override;

};
#endif
