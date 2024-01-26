#ifndef M_BASE_H
#define M_BASE_H

#include "../LIBS/Eigen/Dense"
#include "../CORE/typedefs.h"
#include "functions/function_base.h"

class M_Base {

    public:
        //Normaliser norm;    // TODO?
        virtual ~M_Base() {};

        /** \brief Predict local energy of an atom or bond energy.
         *
         *  The result depends on how aed is computed.
         *
         *  If it is computed between a pair of atoms than the result is a bond energy.
         *
         *  If aed contains sum over all nearest neighbours than the result is
         *  a local atomic energy \f$ E_i \f$.
         * */
        virtual double epredict(const aed_rctype &aed)=0;

        /** \brief Predict force between a pair of atoms in a k-direction. */
        virtual double fpredict(const fd_type &fdij, const aed_rctype &aedi, size_t k)=0;

        /** \brief Predict force between a pair of atoms. */
        virtual force_type fpredict(const fd_type &fdij,
                const aed_rctype &aedi)=0;


};
//template<> inline Registry<M_Base,Function_Base&,Config&>::Map Registry<M_Base,Function_Base&,Config&>::registry{};
#endif
