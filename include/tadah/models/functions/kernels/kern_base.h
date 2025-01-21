#ifndef KERN_BASE_H
#define KERN_BASE_H

#include <tadah/core/core_types.h>
#include <tadah/models/functions/function_base.h>

#include <iostream>

/** \brief Abstract class to be used as a base for all kernels.
 *
 *  - b = basis vector
 *  - af = atomic energy descriptor
 *  - ff = force descriptor
 *  - all derivatives are defined wrt to the second argument
 */
class Kern_Base: public virtual Function_Base {
    public:
        Matrix basis;
        Kern_Base();
        Kern_Base(const Config &c);
        virtual ~Kern_Base();

        // /**
        //  * Calculate kernel value given two vectors
        //  *
        //  * - b = basis vector
        //  * - af = atomic energy descriptor
        //  */
        virtual double operator() (const aed_type &b, const aed_type &af) const=0;
        //
        // /**
        //  * Calculate the kernel derivative wrt to the second argument
        //  *
        //  * - b = basis vector
        //  * - af = atomic energy descriptor
        //  */
        virtual aed_type derivative(const aed_type &b, const aed_type &af) const=0;

        /**
         * Calculate inner product of the kernel derivative
         * wrt to the second argument (b wrt af) with the force descriptor
         *
         * - b = basis vector
         * - af = atomic energy descriptor
         * - ff = force descriptor (TODO i-th dir component of it)
         *   TODO consider calculating all 3-directions at once
         */
        virtual double prime(const aed_type &b, const aed_type &af,
                const aed_type &ff) const=0;

        /** \brief Set basis for calculations */
        virtual void set_basis(const Matrix b);

        /** \brief Return basis stored by this object */
        virtual const Matrix &get_basis() const;
        /** \brief Return basis stored by this object */
        virtual Matrix get_basis();

        virtual double epredict(const t_type &kweights, const aed_type &aed) const;
        virtual double fpredict(const t_type &kweights, const fd_type &fdij,
                const aed_type &aedi, const size_t k) const;
        virtual force_type fpredict(const t_type &kweights, const fd_type &fdij,
                const aed_type &aedi) const;

        virtual void read_basis_from_config(const Config &config,
                Matrix &b);
};
#endif
