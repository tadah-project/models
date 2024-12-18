#include <tadah/models/functions/kernels/kern_linear.h>

Kern_Linear::Kern_Linear() {}
Kern_Linear::Kern_Linear (const Config &c):
    Function_Base(c),
    Kern_Base(c)
{}
double Kern_Linear::operator() (const aed_type& b, const aed_type& af) const
{
    return b*af;;
}
aed_type Kern_Linear::derivative(const aed_type& b, const aed_type&) const
{
    return b;
}
double Kern_Linear::prime(const aed_type& b, const aed_type& , const aed_type& ff) const
{
    // return derivative(b,af).transpose() * ff;
    return b * ff;
}
std::string Kern_Linear::get_label() const
{
    return label;
}
double Kern_Linear::epredict(const t_type &weights, const aed_type& aed) const
{
    return weights*aed;
}
double Kern_Linear::fpredict(const t_type &weights, const fd_type &fdij,
        const aed_type& , const size_t k) const
{
    return -(fdij(k) * weights);
}
force_type Kern_Linear::fpredict(const t_type &weights, const fd_type &fdij,
        const aed_type& ) const
{
    return -(fdij * weights);
}
