#include <tadah/models/functions/kernels/kern_rbf.h>

Kern_RBF::Kern_RBF () {}
Kern_RBF::Kern_RBF (const Config &c):
    Function_Base(c),
    Kern_Base(c),
    sigma(c.get<double>("MPARAMS")),
    gamma(1.0/(sigma*sigma))
{
    if (get_verbose()) std::cout << std::endl << label << " | sigma: " << sigma
        << " | gamma: " << gamma << std::endl;

    read_basis_from_config(c,basis);
}
double Kern_RBF::operator() (const aed_type& b, const aed_type& af) const
{
    const aed_type v = (af-b);
    return exp(-0.5*gamma*(v*v));
}
aed_type Kern_RBF::derivative(const aed_type& b, const aed_type& af) const
{
    return  gamma*(b-af)*(*this)(b,af);
}
double Kern_RBF::prime(const aed_type& b, const aed_type& af, const aed_type& ff) const
{
  // std::cout << "f:     " << ff << std::endl;
    return derivative(b, af)*ff;
}
std::string Kern_RBF::get_label() const
{
    return label;
}
