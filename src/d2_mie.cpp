#include "tadah/models/descriptors/d2/d2_base.h"
#include <tadah/models/descriptors/d2/d2_mie.h>
namespace tadah {
namespace models {


D2_MIE::D2_MIE() {
  init();
}
D2_MIE::D2_MIE(const tadah::core::Context &c): D2_Base(c)
{
  if (!c.get<bool>("INIT2B")) return;
  init();
  s=2;
  n = c.get<double>("TYPE2B",1);
  m = c.get<double>("TYPE2B",2);
  if (n<0 || m<0) {
    throw std::runtime_error("Both Mie exponents must by positive\n");
  }
  auto init_atoms = get_init_atoms(c);
  init_for_atoms(init_atoms);
}

void D2_MIE::calc_aed(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  tadah::core::aed_type &aed,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  double rn_inv = 1.0/pow(rij,n);
  double rm_inv = 1.0/pow(rij,m);
  double fc_ij = fcut->calc(rij);

  aed(fidx) -= scale*weights[Zj]*rn_inv*fc_ij;
  aed(fidx+1) += scale*weights[Zj]*rm_inv*fc_ij;
}
void D2_MIE::calc_dXijdri(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  tadah::core::fd_type &fd_ij,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);

  fd_ij(fidx,0)   = scale*weights[Zj]*(n/pow(rij,n+1)*fc_ij - fcp_ij/pow(rij,n));
  fd_ij(fidx+1,0) = scale*weights[Zj]*(-m/pow(rij,m+1)*fc_ij + fcp_ij/pow(rij,m));
}
void D2_MIE::calc_all(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  tadah::core::aed_type &aed,
  tadah::core::fd_type &fd_ij,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;

  double rn_inv = 1.0/pow(rij,n);
  double rm_inv = 1.0/pow(rij,m);

  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);

  aed(fidx) -= scale*weights[Zj]*rn_inv*fc_ij;
  aed(fidx+1) += scale*weights[Zj]*rm_inv*fc_ij;

  fd_ij(fidx,0) = scale*weights[Zj]*(n*rn_inv/rij*fc_ij - fcp_ij*rn_inv);
  fd_ij(fidx+1,0) = scale*weights[Zj]*(-m*rm_inv/rij*fc_ij + fcp_ij*rm_inv);

}
std::string D2_MIE::label() const {
  return lab;
}
void D2_MIE::init() {
  nparams=2;
}
bool D2_MIE::reinit(const tadah::core::Context& ctx) {
  double temp;
  temp = ctx.get<double>("TYPE2B", 1); if (temp != -1) n = temp;
  temp = ctx.get<double>("TYPE2B", 2); if (temp != -1) m = temp;
  if (fcut) fcut->set_rcut(ctx.get<double>("RCUT2BMAX"));
  return false;  // s = 2 always
}

} // namespace models
} // namespace tadah

