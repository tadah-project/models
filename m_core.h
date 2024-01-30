#ifndef M_Core_H
#define M_Core_H

#include "../CORE/typedefs.h"

class M_Core {
    public:
        int verbose;
        bool trained=false;
        t_type weights;

        virtual ~M_Core() {}


        bool is_trained() const {
            return trained;
        }
        const t_type &get_weights() const {
            return weights;
        }
        void set_weights(const t_type w) {
            weights=w;
        }
        virtual double predict(const rvec &v)const=0;
};
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
        virtual double epredict(const aed_rctype &aed)=0;

        /** \brief Predict force between a pair of atoms in a k-direction. */
        virtual double fpredict(const fd_type &fdij, const aed_rctype &aedi, size_t k)=0;

        /** \brief Predict force between a pair of atoms. */
        virtual force_type fpredict(const fd_type &fdij,
                const aed_rctype &aedi)=0;
};
class M_Train {
    public:
        virtual ~M_Train() {}
        virtual void train(phi_type &Phi, const t_type &T)=0;
        virtual t_type get_weights_uncertainty()const=0;
};
#endif
