#ifndef KERN_LINEAR_H
#define KERN_LINEAR_H
#include "kern_base.h"
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
        std::string get_label();

        /**
         * Label used for this class
         */
        std::string label = "Kern_Linear";
        double operator() (const aed_rctype& b, const aed_rctype& af);
        aed_type derivative(const aed_rctype& b, const aed_rctype& );
        double prime(const aed_rctype& b, const aed_rctype& af, const aed_rctype& ff);
        void set_basis(const aeds_type ) {}

        double epredict(const t_type &weights, const aed_rctype& aed);
        double fpredict(const t_type &weights, const fd_type &fdij,
                const aed_rctype& aedi, const size_t k);
        force_type fpredict(const t_type &weights, const fd_type &fdij,
                const aed_rctype& aedi);
};
#endif
