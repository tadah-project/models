#include <tadah/models/functions/function_base.h>

Function_Base::Function_Base(): verbose(0) {
}
Function_Base::Function_Base(const Config &c) {
  verbose = (c.exist("VERBOSE")) ? c.get<int>("VERBOSE") : 0;
}
void Function_Base::set_verbose(int v) { verbose=v; }
int Function_Base::get_verbose() { return verbose; }

Function_Base::~Function_Base() {}

double Function_Base::operator() (const aed_type& , const aed_type& ) const { return 0.0; }
aed_type Function_Base::derivative(const aed_type& , const aed_type& ) const { return aed_type(); }
double Function_Base::prime(const aed_type& , const aed_type& ,
                            const aed_type& ) const { return 0.0; }
void Function_Base::set_basis(const Matrix ) {}
