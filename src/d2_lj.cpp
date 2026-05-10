#include "tadah/models/descriptors/d2/d2_base.h"
#include <tadah/models/descriptors/d2/d2_lj.h>
#include <tadah/models/descriptors/closed_form_meta.h>
namespace tadah {
namespace models {


D2_LJ::D2_LJ() {
  init();
}
D2_LJ::D2_LJ(const tadah::core::Context &c): D2_Base(c)
{
  if (!c.get<bool>("INIT2B")) return;
  init();
  s=2;

  auto init_atoms = get_init_atoms(c);
  init_for_atoms(init_atoms);
}

void D2_LJ::calc_aed(
  const int Zi,
  const int Zj,
  const double rij,
  const double rij_sq,
  tadah::core::aed_type &aed,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  double r2_inv = 1.0/rij_sq;
  double r6_inv = r2_inv*r2_inv*r2_inv;
  double fc_ij = fcut->calc(rij);

  aed(fidx) -= scale*weights[Zj]*(r6_inv*fc_ij);
  aed(fidx+1) += scale*weights[Zj]*(r6_inv*r6_inv*fc_ij);
}
void D2_LJ::calc_dXijdri(
  const int Zi,
  const int Zj,
  const double rij,
  const double rij_sq,
  tadah::core::fd_type &fd_ij,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;

  double r2_inv = 1.0/rij_sq;
  double r6_inv = r2_inv*r2_inv*r2_inv;

  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);
  fd_ij(fidx,0) =scale*weights[Zj]*( 6.0*r6_inv*r2_inv*rij*fc_ij - fcp_ij*r6_inv);
  fd_ij(fidx+1,0) =scale*weights[Zj]*( -12.0*r6_inv*r6_inv*r2_inv*rij*fc_ij + fcp_ij*r6_inv*r6_inv);
}
void D2_LJ::calc_all(
  const int Zi,
  const int Zj,
  const double rij,
  const double rij_sq,
  tadah::core::aed_type &aed,
  tadah::core::fd_type &fd_ij,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  double r2_inv = 1.0/rij_sq;
  double r6_inv = r2_inv*r2_inv*r2_inv;

  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);
  aed(fidx) -= scale*weights[Zj]*r6_inv*fc_ij;
  aed(fidx+1) += scale*weights[Zj]*r6_inv*r6_inv*fc_ij;

  fd_ij(fidx,0) = scale*weights[Zj]*(6.0*r6_inv*r2_inv*rij*fc_ij - fcp_ij*r6_inv);
  fd_ij(fidx+1,0) = scale*weights[Zj]*(-12.0*r6_inv*r6_inv*r2_inv*rij*fc_ij + fcp_ij*r6_inv*r6_inv);
}
std::string D2_LJ::label() const {
  return lab;
}
void D2_LJ::init() {
  nparams=0;
}

std::optional<DescriptorMeta> D2_LJ::meta() const {
  return make_closed_form_meta("D2_LJ", "RCUT2B");
}

} // namespace models
} // namespace tadah

