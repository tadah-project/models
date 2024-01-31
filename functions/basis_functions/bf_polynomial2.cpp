#include "bf_polynomial2.h"

Registry<Function_Base>::Register<BF_Polynomial2> BF_Polynomial2_1( "BF_Polynomial2" );
Registry<Function_Base,Config&>::Register<BF_Polynomial2> BF_Polynomial2_2( "BF_Polynomial2" );

BF_Polynomial2::BF_Polynomial2() {}
BF_Polynomial2::BF_Polynomial2(const Config &c): Function_Base(c)
{}
std::string BF_Polynomial2::get_label()
{
    return label;
}
double BF_Polynomial2::epredict(const t_type &weights, const aed_rctype& aed)
{
    size_t b=0;
    double res=0.0;
    for (long int i=0; i<aed.size(); ++i) {
        for (long int ii=i; ii<aed.size(); ++ii) {
            res+= weights(b++)*aed(i)*aed(ii);
        }
    }
    return res;
}
double BF_Polynomial2::fpredict(const t_type &weights, const fd_type &fdij,
        const aed_rctype& aedi, const size_t k)
{
    double res=0.0;
    size_t b=0;
    for (long int i=0; i<fdij.rows(); ++i) {
        for (long int ii=i; ii<fdij.rows(); ++ii) {
            res -= weights(b)*(fdij(i,k)*aedi(ii) + fdij(ii,k)*aedi(i));
            b++;
        }
    }
    return res;
}
force_type BF_Polynomial2::fpredict(const t_type &weights, const fd_type &fdij,
        const aed_rctype& aedi)
{
    force_type v;
    v.setZero();
    size_t b=0;
    for (long int i=0; i<fdij.rows(); ++i) {
        for (long int ii=i; ii<fdij.rows(); ++ii) {
            for (size_t k=0; k<3; ++k) {
                v(k) -= weights(b)*(fdij(i,k)*aedi(ii) + fdij(ii,k)*aedi(i));
            }
            b++;
        }
    }
    return v;
}
