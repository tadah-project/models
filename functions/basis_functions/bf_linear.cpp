#include "bf_linear.h"

Registry<Function_Base>::Register<BF_Linear> BF_Linear_1( "BF_Linear" );
Registry<Function_Base,Config&>::Register<BF_Linear> BF_Linear_2( "BF_Linear" );


BF_Linear::BF_Linear() {}
BF_Linear::BF_Linear(const Config &c): Function_Base(c) {}
std::string BF_Linear::get_label() {
    return label;
}
double BF_Linear::epredict(const t_type &weights, const aed_rctype& aed)
{
    return weights.transpose()*aed;
}
double BF_Linear::fpredict(const t_type &weights, const fd_type &fdij,
        const aed_rctype& , const size_t k)
{
    return -fdij.col(k).transpose() * weights;
}
force_type BF_Linear::fpredict(const t_type &weights, const fd_type &fdij,
        const aed_rctype& )
{
    return -Vec3d(fdij.transpose() * weights);
}
