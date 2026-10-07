#ifndef KERN_BASE_H
#define KERN_BASE_H

#include <tadah/core/core_types.h>
#include <tadah/models/functions/function_base.h>

#include <iostream>

namespace tadah {
namespace models {

/** \brief Abstract class to be used as a base for all kernels.
 *
 *  - b = basis vector
 *  - af = atomic energy descriptor
 *  - ff = force descriptor
 *  - all derivatives are defined wrt to the second argument
 */
class Kern_Base: public virtual Function_Base {
    public:
        tadah::core::Matrix basis;
        Kern_Base();
        Kern_Base(const tadah::core::Context &c);
        virtual ~Kern_Base();

        // /**
        //  * Calculate kernel value given two vectors
        //  *
        //  * - b = basis vector
        //  * - af = atomic energy descriptor
        //  */
        virtual double operator() (const tadah::core::aed_type &b, const tadah::core::aed_type &af) const=0;
        //
        // /**
        //  * Calculate the kernel derivative wrt to the second argument
        //  *
        //  * - b = basis vector
        //  * - af = atomic energy descriptor
        //  */
        virtual tadah::core::aed_type derivative(const tadah::core::aed_type &b, const tadah::core::aed_type &af) const=0;

        /**
         * Calculate inner product of the kernel derivative
         * wrt to the second argument (b wrt af) with the force descriptor
         *
         * - b = basis vector
         * - af = atomic energy descriptor
         * - ff = force descriptor (TODO i-th dir component of it)
         *   TODO consider calculating all 3-directions at once
         */
        virtual double prime(const tadah::core::aed_type &b, const tadah::core::aed_type &af,
                const tadah::core::aed_type &ff) const=0;

        /** \brief Closed-form derivative of K(b, af) with respect to a
         *  named kernel hyperparameter.
         *
         *  Used by HPO analytical / hybrid gradients (Phase 2) — parallel
         *  to D2_Base::calc_dDdtheta on the kernel arm. The default
         *  returns false: "this kernel does not own theta_name". Concrete
         *  kernels override for the (theta_name, arg_pos) pairs they
         *  implement.
         *
         *  Calling convention:
         *   - theta_name is the OPTIM key consumed by the kernel
         *     constructor; today the only registered key is "MPARAMS".
         *   - arg_pos is the 1-based numerical-token position within the
         *     key's value list (matching the OPTIM (n) parser
         *     convention in HPO; skips non-numerical tokens). See
         *     HPOTargets::ui_offset.
         *   - dK_out receives the scalar dK / d theta on return.
         *   - Return true iff this call wrote dK_out. The aggregator
         *     uses the bool to detect "no kernel owned this theta" and
         *     route via FD fallback (HYBRID) or raise (ANALYTICAL).
         */
        virtual bool calc_dKdtheta(
            const std::string& /*theta_name*/,
            int                /*arg_pos*/,
            const tadah::core::aed_type& /*b*/,
            const tadah::core::aed_type& /*af*/,
            double& /*dK_out*/) const {
          return false;
        }

        /** \brief Set basis for calculations */
        virtual void set_basis(const tadah::core::Matrix b);

        /** \brief One weight per basis vector (basis.cols()), whatever
         *  dsize is; 0 before the kernel has a basis. */
        std::size_t weights_count(std::size_t dsize) const override;

        /** \brief The length of the basis vectors (basis.rows()). */
        std::size_t basis_vector_size() const override;

        /** \brief Return basis stored by this object */
        virtual const tadah::core::Matrix &get_basis() const;
        /** \brief Return basis stored by this object */
        virtual tadah::core::Matrix get_basis();

        virtual double epredict(const tadah::core::t_type &kweights, const tadah::core::aed_type &aed) const;
        virtual double fpredict(const tadah::core::t_type &kweights, const tadah::core::fd_type &fdij,
                const tadah::core::aed_type &aedi, const size_t k) const;
        virtual tadah::core::force_type fpredict(const tadah::core::t_type &kweights, const tadah::core::fd_type &fdij,
                const tadah::core::aed_type &aedi) const;

        /** \brief Read the basis vectors of a potential: SBASIS vectors
         *  stored column after column in BASIS.
         *
         *  Without BASIS nothing is read and \p b is left as it is (a
         *  training configuration names SBASIS, the number of basis
         *  vectors to draw, and has no BASIS). A BASIS is refused, naming
         *  the keys, when SBASIS is missing or 0, or when the number of
         *  BASIS values is not a multiple of SBASIS. Whether the vectors
         *  have the descriptor's length is checked where that length is
         *  known (M_KRR in MLIP, pair_style tadah). */
        virtual void read_basis_from_context(const tadah::core::Context &context,
                tadah::core::Matrix &b);

    protected:
        /** \brief Function_Base::intercept_columns() of a kernel
         *  K(b, x) = f(b . x): 1 when the descriptor has a bias intercept
         *  (tadah::core::n_intercept_columns() > 0) and basis column 0 is
         *  zero but for its first (bias) component, else 0. */
        size_t dot_product_intercept_columns(const tadah::core::Context &c) const;
};
} // namespace models
} // namespace tadah

#endif
