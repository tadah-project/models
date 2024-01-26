#include "kern_base.h"

Kern_Base::~Kern_Base() {}
void Kern_Base::set_basis(const aeds_type b)
{
    basis = b;
}
double Kern_Base::epredict(const t_type &kweights, const aed_rctype &aed)
{
    double energy=0;
    for (long int b=0; b<basis.cols(); ++b) {
        energy += kweights[b]*(*this)(basis.col(b),aed);
    }
    return energy;
}
double Kern_Base::fpredict(const t_type &kweights, const fd_type &fdij,
        const aed_rctype &aedi, const size_t k)
{
    double res=0.0;
    for (long int b=0; b<basis.cols(); ++b) {
        res -= kweights[b]*(*this).prime(basis.col(b),aedi,fdij.col(k));
    }
    return res;
}
force_type Kern_Base::fpredict(const t_type &kweights, const fd_type &fdij,
        const aed_rctype &aedi)
{
    force_type v(3);
    v.setZero();
    for (long int b=0; b<basis.cols(); ++b) {
        for (size_t k=0; k<3; ++k) {
            v(k) -= kweights[b]*(*this).prime(basis.col(b),aedi,fdij.col(k));
        }
    }
    return v;
}
