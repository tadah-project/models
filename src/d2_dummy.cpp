#include "tadah/models/descriptors/d2/d2_base.h"
#include <tadah/models/descriptors/d2/d2_dummy.h>
namespace tadah {
namespace models {


D2_Dummy::D2_Dummy() {
  init();
}
D2_Dummy::D2_Dummy(tadah::core::Config &c): D2_Base(c)
{
  init();
}

void D2_Dummy::calc_aed(
  const int,
  const int,
  const double,
  const double,
  tadah::core::aed_type&,
  const double) {}

void D2_Dummy::calc_dXijdri(
  const int,
  const int,
  const double,
  const double,
  tadah::core::fd_type &,
  const double) {}
void D2_Dummy::calc_all(
  const int,
  const int,
  const double,
  const double,
  tadah::core::aed_type & ,
  tadah::core::fd_type &,
  const double) {}
std::string D2_Dummy::label() { return lab; }
void D2_Dummy::init() {
  nparams=0;
}

} // namespace models
} // namespace tadah

