#ifndef M_BLR_Core_H
#define M_BLR_Core_H

#include "linear_regressor.h"
#include "../CORE/config/config.h"
#include "functions/function_base.h"
#include "functions/basis_functions/bf_base.h"
#include "m_core.h"

#include <iostream>

template
<class BF=Function_Base&>
class M_BLR_Core:
    public virtual M_Core
{
    public:
        mat Sigma;
        Config &config;
        BF bf;
        M_BLR_Core(Config &c):
            config(c),
            bf(c)
    {
        static_assert(std::is_base_of<BF_Base, BF>::value,
                "\nThe provided Basis function is not of a BasisFunction type.\n\
                For BLR use basis functions defined in basis_functions.h.\n\
                If you want to use this kernel functions use M_KRR class instead.\n\
                e.g. Use BF_Linear instead of Kern_Linear\n.");

        init();
    }

        M_BLR_Core(BF &bf, Config &c):
            config(c),
            bf(bf)
    {
        if (dynamic_cast<BF_Base*>(&bf) == nullptr)
            throw std::invalid_argument("Provided object is not of BF_Base type");

        init();
    }
        double predict(const rvec &v) const{
            return bf.epredict(get_weights(),v);
        };


    private:
        void init() {
            if (config.exist("WEIGHTS")) {
                weights.resize(config.size("WEIGHTS"));
                config.get<t_type>("WEIGHTS",weights);
                trained=true;
            }

            verbose=(config.get<int>("VERBOSE"));
        }
};
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

        void train(phi_type &Phi, const t_type &T) {
            if(trained) {
                throw std::runtime_error("This object is already trained!");
            }
            LinearRegressor::train(config,Phi, T,weights,Sigma);
            trained=true;
        }
        t_type get_weights_uncertainty() const{
            double lambda=config.template get<double>("LAMBDA");
            if(lambda >= 0) throw std::runtime_error(
                    "Sigma matrix is only computed for LAMBDA < 0");
            return Sigma.diagonal();
        }
};
#endif
