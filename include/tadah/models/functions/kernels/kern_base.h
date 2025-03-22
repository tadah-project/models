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

        /** \brief Set basis for calculations */
        virtual void set_basis(const tadah::core::Matrix b);

        /** \brief Return basis stored by this object */
        virtual const tadah::core::Matrix &get_basis() const;
        /** \brief Return basis stored by this object */
        virtual tadah::core::Matrix get_basis();

        virtual double epredict(const tadah::core::t_type &kweights, const tadah::core::aed_type &aed) const;
        virtual double fpredict(const tadah::core::t_type &kweights, const tadah::core::fd_type &fdij,
                const tadah::core::aed_type &aedi, const size_t k) const;
        virtual tadah::core::force_type fpredict(const tadah::core::t_type &kweights, const tadah::core::fd_type &fdij,
                const tadah::core::aed_type &aedi) const;

        virtual void read_basis_from_config(const tadah::core::Context &config,
                tadah::core::Matrix &b);
};
} // namespace models
} // namespace tadah

#endif
