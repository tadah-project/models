#ifndef M_PREDICT_H
#define M_PREDICT_H

#include <tadah/core/core_types.h>

namespace tadah {
namespace models {

class M_Predict {
    public:

        virtual ~M_Predict() {}

        /** \brief Predict local energy of an atom or bond energy.
         *
         *  The result depends on how aed is computed.
         *
         *  If it is computed between a pair of atoms than the result is a bond energy.
         *
         *  If aed contains sum over all nearest neighbours than the result is
         *  a local atomic energy \f$ E_i \f$.
         * */
        virtual double epredict(const aed_type &aed) const=0;

        /** \brief Predict force between a pair of atoms in a k-direction. */
        virtual double fpredict(const fd_type &fdij, const aed_type &aedi, size_t k) const=0;

        /** \brief Predict force between a pair of atoms. */
} // namespace models
} // namespace tadah

        virtual force_type fpredict(const fd_type &fdij,
                const aed_type &aedi) const=0;
};
#endif
