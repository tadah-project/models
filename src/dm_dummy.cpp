#include "tadah/models/descriptors/dm/dm_base.h"
#include <tadah/models/descriptors/dm/dm_dummy.h>

DM_Dummy::DM_Dummy() {
  init();
}
DM_Dummy::DM_Dummy(Config&c): DM_Base(c) {
  init();
}
void DM_Dummy::calc_aed(
  rho_type&,
  aed_type2 & ) {}
void DM_Dummy::calc_dXijdri_dXjidri(
  const int,
  const int,
  const double,
  const double,
  const Vec3d &,
  rho_type&,
  rho_type&,
  fd_type &,
  const double) {}
void DM_Dummy::calc_dXijdri(
  const int,
  const double,
  const double,
  const Vec3d &,
  rho_type&,
  fd_type &,
  const double) {}
std::string DM_Dummy::label() {return lab;}
void DM_Dummy::init_rhoi(rho_type&) {}
void DM_Dummy::calc_rho(
  const int,
  const double,
  const double,
  const Vec3d &,
  rho_type&,
  const double) {}

void DM_Dummy::init() {
  nparams=0;
}
