#ifndef M_TRAIN_H
#define M_TRAIN_H

#include "../CORE/core_types.h"
//#include "../CORE/eigen_types.h"

class M_Train {
    public:
        virtual ~M_Train() {}
        virtual void train(phi_type &Phi, const t_type &T)=0;
};
#endif
