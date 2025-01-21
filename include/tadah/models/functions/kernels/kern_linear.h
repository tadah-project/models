#ifndef KERN_LINEAR_H
#define KERN_LINEAR_H

#include <tadah/models/functions/kernels/kern_base.h>

/**
 * Linear kernel also knows as dot product kernel
 *
 * Defined for two vectors **x** and **y**:
 *
 * \f[
 * K(\mathbf{x}, \mathbf{y}) = \mathbf{x}^T \mathbf{y} = \sum_i x^{(i)} y^{(i)}
 * \f]
 *
 *  @see Kern_Base BF_Linear
 */
class Kern_Linear :  public virtual Kern_Base {
    public:
        Kern_Linear ();
        Kern_Linear (const Config &c);
        std::string get_label() const override;

        /**
         * Label used for this class
         */
        std::string label = "Kern_Linear";
        double operator() (const aed_type& b, const aed_type& af) const override;
        aed_type derivative(const aed_type& b, const aed_type& ) const override;
        double prime(const aed_type& b, const aed_type& af, const aed_type& ff) const override;
        void set_basis(const Matrix ) override {}

        double epredict(const t_type &weights, const aed_type& aed) const override;
        double fpredict(const t_type &weights, const fd_type &fdij,
                const aed_type& aedi, const size_t k) const override;
        force_type fpredict(const t_type &weights, const fd_type &fdij,
                const aed_type& aedi) const override;
};
#endif
