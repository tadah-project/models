#include <tadah/models/functions/kernels/kern_polynomial.h>
namespace tadah {
namespace models {


Kern_Polynomial::Kern_Polynomial () {}
Kern_Polynomial::Kern_Polynomial (const tadah::core::Context &context):
    Function_Base(context),
    Kern_Base(context),
    d(context.get<int>("MPARAMS",0)),
    gamma(context.get<double>("MPARAMS",1)),
    c(context.get<double>("MPARAMS",2))
{
    if (get_verbose()) std::cout << std::endl << label << " | degree: " << d
        <<" | gamma: " << gamma << " | c: " << c << std::endl;
    read_basis_from_context(context,basis);
}
double Kern_Polynomial::operator() (const tadah::core::aed_type& b, const tadah::core::aed_type& af) const
{
    return std::pow(gamma*(b*af)+c,d);
}
tadah::core::aed_type Kern_Polynomial::derivative(const tadah::core::aed_type& b, const tadah::core::aed_type& af) const
{
    return  d*gamma*std::pow(gamma*(b*af)+c,d-1)*b;
}
double Kern_Polynomial::prime(const tadah::core::aed_type& b, const tadah::core::aed_type& af, const tadah::core::aed_type& ff) const
{
    return derivative(b, af)*ff;
}
std::string Kern_Polynomial::get_label() const
{
    return label;
}

} // namespace models
} // namespace tadah

