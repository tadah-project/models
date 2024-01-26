#include "kern_polynomial.h"

Registry<Function_Base>::Register<Kern_Polynomial> Kern_Polynomial_1( "Kern_Polynomial" );
Registry<Function_Base,Config&>::Register<Kern_Polynomial> Kern_Polynomial_2( "Kern_Polynomial" );

Kern_Polynomial::Kern_Polynomial () {}
Kern_Polynomial::Kern_Polynomial (const Config &c):
    Function_Base(c),
    d(c.get<int>("MPARAMS",0)),
    gamma(c.get<double>("MPARAMS",1)),
    c(c.get<double>("MPARAMS",2))
{
    if (get_verbose()) std::cout << std::endl << label << " | degree: " << d
        <<" | gamma: " << gamma << " | c: " << c << std::endl;
}
double Kern_Polynomial::operator() (const aed_rctype& b, const aed_rctype& af)
{
    return std::pow(gamma*(b.dot(af))+c,d);
}
aed_type Kern_Polynomial::derivative(const aed_rctype& b, const aed_rctype& af)
{
    return  d*gamma*std::pow(gamma*(b.dot(af))+c,d-1)*b;
}
double Kern_Polynomial::prime(const aed_rctype& b, const aed_rctype& af, const aed_rctype& ff)
{
    return derivative(b, af).transpose()*ff;
}
std::string Kern_Polynomial::get_label()
{
    return label;
}
