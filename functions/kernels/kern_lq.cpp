#include "kern_lq.h"


Kern_LQ::Kern_LQ () {}
Kern_LQ::Kern_LQ (const Config &c):
    Function_Base(c)
{
    if (get_verbose()) std::cout << std::endl << label << std::endl;
}
double Kern_LQ::operator() (const aed_type2& b, const aed_type2& af)
{
    double x = b.transpose()*af;
    return x*x + x;
}
aed_type2 Kern_LQ::derivative(const aed_type2& b, const aed_type2& af)
{
    double temp = 2.0 * b.transpose()*af;
    return temp*b + b;
}
double Kern_LQ::prime(const aed_type2& b, const aed_type2& af, const aed_type2& ff)
{
    return derivative(b, af).transpose()*ff;
}
std::string Kern_LQ::get_label()
{
    return label;
}
