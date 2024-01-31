#include "kern_rbf.h"

Registry<Function_Base>::Register<Kern_RBF> Kern_RBF_1( "Kern_RBF" );
Registry<Function_Base,Config&>::Register<Kern_RBF> Kern_RBF_2( "Kern_RBF" );

Kern_RBF::Kern_RBF () {}
Kern_RBF::Kern_RBF (const Config &c):
    Function_Base(c),
    sigma(c.get<double>("MPARAMS")),
    gamma(1.0/(2.0*sigma*sigma))
{
    if (get_verbose()) std::cout << std::endl << label << " | sigma: " << sigma
        << " | gamma: " << gamma << std::endl;
}
double Kern_RBF::operator() (const aed_type2& b, const aed_type2& af)
{
    const aed_type2 v = af-b;
    return exp(-gamma*(v*v));
}
aed_type2 Kern_RBF::derivative(const aed_type2& b, const aed_type2& af)
{
    return  -2.0*gamma*(af-b)*(*this)(af,b);
}
double Kern_RBF::prime(const aed_type2& b, const aed_type2& af, const aed_type2& ff)
{
    return derivative(b, af).transpose()*ff;
}
std::string Kern_RBF::get_label()
{
    return label;
}
