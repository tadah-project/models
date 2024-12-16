#include <tadah/models/descriptors/d3/d3_dummy.h>

D3_Dummy::D3_Dummy() {
  init();
}
D3_Dummy::D3_Dummy(Config &) {
  init();
}

void D3_Dummy::calc_aed(
        const size_t,
        const double,
        const double,
        const double,
        const double,
        aed_type2& ) {}
void D3_Dummy::calc_fd(
        const size_t,
        const double,
        const double,
        const double,
        const double,
        const double,
        const double,
        fd_type &) {}
void D3_Dummy::calc_all(
        const size_t,
        const double,
        const double,
        const double,
        const double,
        const double,
        const double,
        aed_type2& ,
        fd_type &) {}
size_t D3_Dummy::size() { return s;}
std::string D3_Dummy::label() {return lab;}
void D3_Dummy::init() {
  nparams=0;
}
