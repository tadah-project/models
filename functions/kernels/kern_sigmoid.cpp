#include "kern_sigmoid.h"

Registry<Function_Base>::Register<Kern_Sigmoid> Kern_Sigmoid_1( "Kern_Sigmoid" );
Registry<Function_Base,Config&>::Register<Kern_Sigmoid> Kern_Sigmoid_2( "Kern_Sigmoid" );

Kern_Sigmoid::Kern_Sigmoid () {}
Kern_Sigmoid::Kern_Sigmoid (const Config &c):
    Function_Base(c),
    gamma(c.get<double>("MPARAMS",0)),
    c(c.get<double>("MPARAMS",1))
{
    if (get_verbose()) std::cout << std::endl << label << " | gamma: " << gamma 
        << " | c: " << c << std::endl;
}
double Kern_Sigmoid::operator() (const aed_rctype& b, const aed_rctype& af)
{
    return std::tanh(gamma*(b.dot(af))+c);
}
aed_type Kern_Sigmoid::derivative(const aed_rctype& b, const aed_rctype& af)
{
    return gamma*(1.0-std::pow((*this)(af,b),2))*b;
}
double Kern_Sigmoid::prime(const aed_rctype& b, const aed_rctype& af, const aed_rctype& ff)
{
    return derivative(b, af).transpose()*ff;
}
std::string Kern_Sigmoid::get_label()
{
    return label;
}
