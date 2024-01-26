#ifndef KERN_BASE_H
#define KERN_BASE_H

#include "../../../CORE/typedefs.h"
#include "../function_base.h"

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
        aeds_type basis;
        virtual ~Kern_Base();

        /**
         * Calculate kernel value given two vectors
         *
         * - b = basis vector
         * - af = atomic energy descriptor
         */
        virtual double operator() (const aed_rctype &b, const aed_rctype &af)=0;

        /**
         * Calculate the kernel derivative wrt to the second argument
         *
         * - b = basis vector
         * - af = atomic energy descriptor
         */
        virtual aed_type derivative(const aed_rctype &b, const aed_rctype &af)=0;

        /**
         * Calculate inner product of the kernel derivative
         * wrt to the second argument (b wrt af) with the force descriptor
         *
         * - b = basis vector
         * - af = atomic energy descriptor
         * - ff = force descriptor (TODO i-th dir component of it)
         *   TODO consider calculating all 3-directions at once
         */
        virtual double prime(const aed_rctype &b, const aed_rctype &af,
                const aed_rctype &ff)=0;

        /** \brief Set basis for calculations */
        virtual void set_basis(const aeds_type b);

        virtual double epredict(const t_type &kweights, const aed_rctype &aed);
        virtual double fpredict(const t_type &kweights, const fd_type &fdij,
                const aed_rctype &aedi, const size_t k);
        virtual force_type fpredict(const t_type &kweights, const fd_type &fdij,
                const aed_rctype &aedi);
};
#endif
