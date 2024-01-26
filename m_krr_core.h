#ifndef M_KRR_Core_H
#define M_KRR_Core_H

#include "linear_regressor.h"
#include "../CORE/config/config.h"
#include "functions/function_base.h"
#include "functions/kernels/kern_base.h"
#include "m_core.h"
#include "ekm.h"

#include <iostream>

template
<class K=Function_Base&>
class M_KRR_Core: public virtual M_Core {
    public:
        Config &config;
        mat Sigma;
        K kernel;
        EKM<K> ekm;
        M_KRR_Core(Config &c):
            config(c),
            kernel(c),
            ekm(c)
    {
        static_assert(std::is_base_of<Kern_Base, K>::value,
                "\nThe provided Kernel K is not of a Kernel type.\n\
                For KRR use kernels defined in kernel.h.\n\
                e.g. Use Kern_Linear instead of BF_Linear\n.");

        init();
    }

        M_KRR_Core(K &kernel, Config &c):
            config(c),
            kernel(kernel),
            ekm(kernel)
    {
        if (dynamic_cast<Kern_Base*>(&kernel) == nullptr)
            throw std::invalid_argument("Provided object is not of Kern_Base type");

        init();
    }
        double predict(const rvec &v) const{
            return kernel.epredict(get_weights(),v);
        };

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
        t_type get_weights_uncertainty() const{
            double lambda=config.get<double>("LAMBDA");
            if(lambda >= 0) throw std::runtime_error(
                    "Sigma matrix is only computed for LAMBDA < 0");
            return Sigma.diagonal();
        }

    private:
        void init() {
            if (config.exist("WEIGHTS")) {
                weights.resize(config.size("WEIGHTS"));
                config.get<t_type>("WEIGHTS",weights);
                trained=true;
            }

            M_Core::verbose=(config.get<int>("VERBOSE"));
        }
};
#endif
