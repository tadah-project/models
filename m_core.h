#ifndef M_Core_H
#define M_Core_H

#include "../CORE/typedefs.h"

class M_Core {
    public:
        int verbose;
        bool trained=false;
        t_type weights;

        virtual ~M_Core() {}

        virtual double predict(const rvec &v)const=0;
        virtual void train(phi_type &Phi, const t_type &T)=0;
        virtual t_type get_weights_uncertainty()const=0;

        bool is_trained() const {
            return trained;
        }
        const t_type &get_weights() const {
            return weights;
        }
        void set_weights(const t_type w) {
            weights=w;
        }
};
#endif
