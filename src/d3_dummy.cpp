#include <tadah/models/descriptors/d3/d3_dummy.h>
namespace tadah {
namespace models {


D3_Dummy::D3_Dummy() {
  init();
}
D3_Dummy::D3_Dummy(tadah::core::Context &) {
  init();
}

void D3_Dummy::calc_aed(
        const size_t,
        const double,
        const double,
        const double,
        const double,
        tadah::core::aed_type& ) {}
void D3_Dummy::calc_fd(
        const size_t,
        const double,
        const double,
        const double,
        const double,
        const double,
        const double,
        tadah::core::fd_type &) {}
void D3_Dummy::calc_all(
        const size_t,
        const double,
        const double,
        const double,
        const double,
        const double,
        const double,
        tadah::core::aed_type& ,
        tadah::core::fd_type &) {}
size_t D3_Dummy::size() { return s;}
std::string D3_Dummy::label() {return lab;}
void D3_Dummy::init() {
  nparams=0;
}

} // namespace models
} // namespace tadah

