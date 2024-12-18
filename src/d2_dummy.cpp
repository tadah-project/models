#include "tadah/models/descriptors/d2/d2_base.h"
#include <tadah/models/descriptors/d2/d2_dummy.h>

D2_Dummy::D2_Dummy() {
  init();
}
D2_Dummy::D2_Dummy(Config &c): D2_Base(c)
{
  init();
}

void D2_Dummy::calc_aed(
  const int,
  const int,
  const double,
  const double,
  aed_type&,
  const double) {}

void D2_Dummy::calc_dXijdri(
  const int,
  const int,
  const double,
  const double,
  fd_type &,
  const double) {}
void D2_Dummy::calc_all(
  const int,
  const int,
  const double,
  const double,
  aed_type & ,
  fd_type &,
  const double) {}
std::string D2_Dummy::label() { return lab; }
void D2_Dummy::init() {
  nparams=0;
}
