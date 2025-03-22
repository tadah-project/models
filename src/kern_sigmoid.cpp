#include <tadah/models/functions/kernels/kern_sigmoid.h>
namespace tadah {
namespace models {


Kern_Sigmoid::Kern_Sigmoid () {}
Kern_Sigmoid::Kern_Sigmoid (const tadah::core::Context &config):
    Function_Base(config),
    Kern_Base(config),
    gamma(config.get<double>("MPARAMS",0)),
    c(config.get<double>("MPARAMS",1))
{
    if (get_verbose()) std::cout << std::endl << label << " | gamma: " << gamma
        << " | c: " << c << std::endl;

    read_basis_from_config(config,basis);
}
double Kern_Sigmoid::operator() (const tadah::core::aed_type& b, const tadah::core::aed_type& af) const
{
    return std::tanh(gamma*(b*af)+c);
}
tadah::core::aed_type Kern_Sigmoid::derivative(const tadah::core::aed_type& b, const tadah::core::aed_type& af) const
{
    return gamma*(1.0-std::pow((*this)(af,b),2))*b;
}
double Kern_Sigmoid::prime(const tadah::core::aed_type& b, const tadah::core::aed_type& af, const tadah::core::aed_type& ff) const
{
    return derivative(b, af)*ff;
}
std::string Kern_Sigmoid::get_label() const
{
    return label;
}

} // namespace models
} // namespace tadah

