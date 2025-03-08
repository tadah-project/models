#ifndef KERN_LQ_H
#define KERN_LQ_H

#include <tadah/models/functions/kernels/kern_base.h>

namespace tadah {
namespace models {

/**
 * Linear + Quadratic kernel
 *
 * Defined for two vectors **x** and **y**:
 *
 * \f[
 * K(\mathbf{x}, \mathbf{y}) = \mathbf{x}^T \mathbf{y} +
 * \Big(\mathbf{x}^T \mathbf{y} \Big)^2
 * = \sum_i x^{(i)} y^{(i)} + \Big(\sum_i x^{(i)} y^{(i)} \Big)^2
 * \f]
 *
 *  @see Kern_Base Kern_Linear Kern_Quadratic
 */
class Kern_LQ : public virtual Kern_Base {
    public:
        Kern_LQ ();
        Kern_LQ (const Config &c);
        /**
         * Label used for this class
         */
        std::string label = "Kern_LQ";
        double operator() (const aed_type& b, const aed_type& af) const override;
        aed_type derivative(const aed_type& b, const aed_type& af) const override;
        double prime(const aed_type& b, const aed_type& af, const aed_type& ff) const override;
} // namespace models
} // namespace tadah

        std::string get_label() const override;

};
#endif
