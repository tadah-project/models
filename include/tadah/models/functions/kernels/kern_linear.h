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
        std::string get_label() const;

        /**
         * Label used for this class
         */
        std::string label = "Kern_Linear";
        double operator() (const aed_type2& b, const aed_type2& af) const;
        aed_type2 derivative(const aed_type2& b, const aed_type2& ) const;
        double prime(const aed_type2& b, const aed_type2& af, const aed_type2& ff) const;
        void set_basis(const Matrix ) {}

        double epredict(const t_type &weights, const aed_type2& aed) const;
        double fpredict(const t_type &weights, const fd_type &fdij,
                const aed_type2& aedi, const size_t k) const;
        force_type fpredict(const t_type &weights, const fd_type &fdij,
                const aed_type2& aedi) const;
};
#endif
