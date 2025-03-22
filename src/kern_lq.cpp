#include <tadah/models/functions/kernels/kern_lq.h>
namespace tadah {
namespace models {


Kern_LQ::Kern_LQ () {}
Kern_LQ::Kern_LQ (const tadah::core::Context &c):
    Function_Base(c),
    Kern_Base(c)
{
    if (get_verbose()) std::cout << std::endl << label << std::endl;
    read_basis_from_config(c,basis);
}
double Kern_LQ::operator() (const tadah::core::aed_type& b, const tadah::core::aed_type& af) const
{
    double x = b*af;
    return x*x + x;
}
tadah::core::aed_type Kern_LQ::derivative(const tadah::core::aed_type& b, const tadah::core::aed_type& af) const
{
    double temp = 2.0 * b*af;
    return temp*b + b;
}
double Kern_LQ::prime(const tadah::core::aed_type& b, const tadah::core::aed_type& af, const tadah::core::aed_type& ff) const
{
    return derivative(b, af)*ff;
}
std::string Kern_LQ::get_label() const
{
    return label;
}

} // namespace models
} // namespace tadah

