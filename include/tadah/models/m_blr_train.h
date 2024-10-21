#ifndef M_BLR_TRAIN_H
#define M_BLR_TRAIN_H

#include <tadah/models/m_train.h>
#include <tadah/models/m_blr_core.h>
#include <tadah/models/linear_regressor.h>
#include <tadah/models/functions/function_base.h>
#include <tadah/models/functions/basis_functions/bf_base.h>
#include <tadah/core/config.h>

#include <iostream>

template
<class BF=Function_Base&>
class M_BLR_Train:
    public M_BLR_Core<BF>,    // MUST NOT BE VIRUAL!!
    public M_Train
{

    public:
        using M_BLR_Core<BF>::trained;
        using M_BLR_Core<BF>::config;
        using M_BLR_Core<BF>::weights;
        using M_BLR_Core<BF>::Sigma;

        M_BLR_Train(Config &c):
            M_BLR_Core<BF>(c)
    {}

        M_BLR_Train(BF &bf, Config &c):
            M_BLR_Core<BF>(bf, c)
    {}

        void train(phi_type &Phi, t_type &T) {
            if(trained) {
                throw std::runtime_error("This object is already trained!");
            }
            LinearRegressor::train(config,Phi, T,weights,Sigma);
            trained=true;
        }
};
#endif
