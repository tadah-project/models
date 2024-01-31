#include "kern_base.h"

Kern_Base::~Kern_Base() {}
void Kern_Base::set_basis(const aeds_type2 b)
{
    basis = b;
}
double Kern_Base::epredict(const t_type &kweights, const aed_type2 &aed)
{
    double energy=0;
    for (size_t b=0; b<basis.cols(); ++b) {
        energy += kweights[b]*(*this)(basis.col(b),aed);
    }
    return energy;
}
double Kern_Base::fpredict(const t_type &kweights, const fd_type &fdij,
        const aed_type2 &aedi, const size_t k)
{
    double res=0.0;
    for (size_t b=0; b<basis.cols(); ++b) {
        res -= kweights[b]*(*this).prime(basis.col(b),aedi,fdij.col(k));
    }
    return res;
}
force_type Kern_Base::fpredict(const t_type &kweights, const fd_type &fdij,
        const aed_type2 &aedi)
{
    force_type v;
    v.setZero();
    for (size_t b=0; b<basis.cols(); ++b) {
        for (size_t k=0; k<3; ++k) {
            v(k) -= kweights[b]*(*this).prime(basis.col(b),aedi,fdij.col(k));
        }
    }
    return v;
}
