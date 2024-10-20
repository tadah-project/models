#ifndef M_TRAIN_H
#define M_TRAIN_H

#include <tadah/core/core_types.h>

class M_Train {
    public:
        virtual ~M_Train() {}
        virtual void train(phi_type &Phi, t_type &T)=0;
};
#endif
