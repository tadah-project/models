#include <tadah/models/functions/function_base.h>
namespace tadah {
namespace models {


tadah::core::Function_Base::tadah::core::Function_Base(): verbose(0) {
}
tadah::core::Function_Base::tadah::core::Function_Base(const tadah::core::Config &c) {
  verbose = (c.exist("VERBOSE")) ? c.get<int>("VERBOSE") : 0;
}
void tadah::core::Function_Base::set_verbose(int v) { verbose=v; }
int tadah::core::Function_Base::get_verbose() { return verbose; }

tadah::core::Function_Base::~tadah::core::Function_Base() {}

double tadah::core::Function_Base::operator() (const tadah::core::aed_type& , const tadah::core::aed_type& ) const { return 0.0; }
tadah::core::aed_type tadah::core::Function_Base::derivative(const tadah::core::aed_type& , const tadah::core::aed_type& ) const { return tadah::core::aed_type(); }
double tadah::core::Function_Base::prime(const tadah::core::aed_type& , const tadah::core::aed_type& ,
                            const tadah::core::aed_type& ) const { return 0.0; }
void tadah::core::Function_Base::set_basis(const tadah::core::Matrix ) {}

} // namespace models
} // namespace tadah

