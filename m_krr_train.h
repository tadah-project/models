#ifndef M_KRR_TRAIN_H
#define M_KRR_TRAIN_H

#include "linear_regressor.h"
#include "../CORE/config/config.h"
#include "functions/function_base.h"
#include "functions/kernels/kern_base.h"
#include "m_core.h"
#include "m_krr_core.h"
#include "ekm.h"

#include <iostream>

template
<class K=Function_Base&>
class M_KRR_Train:
    public M_KRR_Core<K>,    // MUST NOT BE VIRUAL!!
    public M_Train
{

    public:
        EKM<K> ekm;
        using M_KRR_Core<K>::trained;
        using M_KRR_Core<K>::config;
        using M_KRR_Core<K>::kernel;
        using M_KRR_Core<K>::weights;
        using M_KRR_Core<K>::Sigma;

        M_KRR_Train(Config &c):
            M_KRR_Core<K>(c),
            ekm(c)
    {}

        M_KRR_Train(K &kernel, Config &c):
            M_KRR_Core<K>(kernel, c),
            ekm(kernel)
    {}

        void train(phi_type &Phi, const t_type &T) {
            if(trained) {
                throw std::runtime_error("This object is already trained!");
            }
            if (kernel.get_label()!="Kern_Linear") {
                ekm.project(Phi);
            }
            LinearRegressor::train(config,Phi, T,weights,Sigma);

            if (kernel.get_label()!="Kern_Linear") {
                //kernalize weights
                weights = ekm.KK.transpose()*weights;
            }
            trained=true;
        }
};
#endif
