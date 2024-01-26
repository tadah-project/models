#include "kern_quadratic.h"

Registry<Function_Base>::Register<Kern_Quadratic> Kern_Quadratic_1( "Kern_Quadratic" );
Registry<Function_Base,Config&>::Register<Kern_Quadratic> Kern_Quadratic_2( "Kern_Quadratic" );

Kern_Quadratic::Kern_Quadratic () {}
Kern_Quadratic::Kern_Quadratic (const Config &c): 
    Function_Base(c)
{
    if (get_verbose()) std::cout << std::endl << label << std::endl;
}
double Kern_Quadratic::operator() (const aed_rctype& b, const aed_rctype& af)
{
    double x = b.transpose()*af;
    return  x*x;
}
aed_type Kern_Quadratic::derivative(const aed_rctype& b, const aed_rctype& af)
{
    double temp = 2.0 * b.transpose()*af;
    return b*temp;;
}
double Kern_Quadratic::prime(const aed_rctype& b, const aed_rctype& af, const aed_rctype& ff)
{
    //return dot(derivative(b, af), ff);
    // it is faster this way
    double scale = 2.0*b.transpose()*af;
    return scale*ff.transpose()*b;
}
std::string Kern_Quadratic::get_label()
{
    return label;
}
