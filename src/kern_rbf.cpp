#include <tadah/models/functions/kernels/kern_rbf.h>
namespace tadah {
namespace models {


Kern_RBF::Kern_RBF () {}
Kern_RBF::Kern_RBF (const tadah::core::Config &c):
    Function_Base(c),
    Kern_Base(c),
    sigma(c.get<double>("MPARAMS")),
    gamma(1.0/(sigma*sigma))
{
    if (get_verbose()) std::cout << std::endl << label << " | sigma: " << sigma
        << " | gamma: " << gamma << std::endl;

    read_basis_from_config(c,basis);
}
double Kern_RBF::operator() (const tadah::core::aed_type& b, const tadah::core::aed_type& af) const
{
    const tadah::core::aed_type v = (af-b);
    return exp(-0.5*gamma*(v*v));
}
tadah::core::aed_type Kern_RBF::derivative(const tadah::core::aed_type& b, const tadah::core::aed_type& af) const
{
    return  gamma*(b-af)*(*this)(b,af);
}
double Kern_RBF::prime(const tadah::core::aed_type& b, const tadah::core::aed_type& af, const tadah::core::aed_type& ff) const
{
  // std::cout << "f:     " << ff << std::endl;
    return derivative(b, af)*ff;
}
std::string Kern_RBF::get_label() const
{
    return label;
}

} // namespace models
} // namespace tadah

