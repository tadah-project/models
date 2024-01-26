#include "kern_linear.h"

Registry<Function_Base>::Register<Kern_Linear> Kern_Linear_1( "Kern_Linear" );
Registry<Function_Base,Config&>::Register<Kern_Linear> Kern_Linear_2( "Kern_Linear" );

Kern_Linear::Kern_Linear() {}
Kern_Linear::Kern_Linear (const Config &c):
    Function_Base(c)
{}
double Kern_Linear::operator() (const aed_rctype& b, const aed_rctype& af)
{
    return b.transpose()*af;;
}
aed_type Kern_Linear::derivative(const aed_rctype& b, const aed_rctype&)
{
    return b;
}
double Kern_Linear::prime(const aed_rctype& b, const aed_rctype& , const aed_rctype& ff)
{
    // return derivative(b,af).transpose() * ff;
    return b.transpose() * ff;
}
std::string Kern_Linear::get_label()
{
    return label;
}
double Kern_Linear::epredict(const t_type &weights, const aed_rctype& aed)
{
    return weights.transpose()*aed;
}
double Kern_Linear::fpredict(const t_type &weights, const fd_type &fdij,
        const aed_rctype& , const size_t k)
{
    return -fdij.col(k).transpose() * weights;
}
force_type Kern_Linear::fpredict(const t_type &weights, const fd_type &fdij,
        const aed_rctype& )
{
    return -fdij.transpose() * weights;
}
