#include <tadah/models/functions/kernels/kern_linear.h>
namespace tadah {
namespace models {


Kern_Linear::Kern_Linear() {}
Kern_Linear::Kern_Linear (const tadah::core::Config &c):
    tadah::core::Function_Base(c),
    Kern_Base(c)
{}
double Kern_Linear::operator() (const tadah::core::aed_type& b, const tadah::core::aed_type& af) const
{
    return b*af;;
}
tadah::core::aed_type Kern_Linear::derivative(const tadah::core::aed_type& b, const tadah::core::aed_type&) const
{
    return b;
}
double Kern_Linear::prime(const tadah::core::aed_type& b, const tadah::core::aed_type& , const tadah::core::aed_type& ff) const
{
    // return derivative(b,af).transpose() * ff;
    return b * ff;
}
std::string Kern_Linear::get_label() const
{
    return label;
}
double Kern_Linear::epredict(const tadah::core::t_type &weights, const tadah::core::aed_type& aed) const
{
    return weights*aed;
}
double Kern_Linear::fpredict(const tadah::core::t_type &weights, const tadah::core::fd_type &fdij,
        const tadah::core::aed_type& , const size_t k) const
{
    return -(fdij(k) * weights);
}
tadah::core::force_type Kern_Linear::fpredict(const tadah::core::t_type &weights, const tadah::core::fd_type &fdij,
        const tadah::core::aed_type& ) const
{
    return -(fdij * weights);
}

} // namespace models
} // namespace tadah

