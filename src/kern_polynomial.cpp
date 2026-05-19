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
// Closed-form derivatives of K w.r.t. gamma and c. d is integer and
// is not HPO-tunable. Common factor: d * (gamma * b.af + c)^(d-1).
bool Kern_Polynomial::calc_dKdtheta(
    const std::string& theta_name,
    int                arg_pos,
    const tadah::core::aed_type& b,
    const tadah::core::aed_type& af,
    double& dK_out) const
{
    if (theta_name != "MPARAMS") return false;
    const double dot = b * af;
    const double base = gamma * dot + c;
    const double common = static_cast<double>(d) * std::pow(base, d - 1);
    if (arg_pos == 2) {        // dK/d gamma
        dK_out = common * dot;
        return true;
    }
    if (arg_pos == 3) {        // dK/d c
        dK_out = common;
        return true;
    }
    return false;
}

std::string Kern_Polynomial::get_label() const
{
    return label;
}

} // namespace models
} // namespace tadah

