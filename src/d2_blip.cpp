#include "tadah/models/descriptors/d_base.h"
#include <tadah/models/descriptors/d2/d2_blip.h>
#include <tadah/models/descriptors/d_basis_functions.h>
#include <cstdlib>
#include <stdexcept>
namespace tadah {
namespace models {


D2_Blip::D2_Blip() {
  init();
}
D2_Blip::D2_Blip(const tadah::core::Context &c): D2_Base(c)
{
  if (!c.get<bool>("INIT2B")) return;
  init();

  get_grid(c,"CGRID2B",mius);
  get_grid(c,"SGRID2B",etas);

  if (mius.size()!=etas.size()) {
    throw std::runtime_error("SGRID2B and CGRID2B arrays differ in size.\n");
  }

  for (const auto &e:etas) {
    if(std::abs(e)<std::numeric_limits<double>::min())
      throw std::runtime_error("At least one of SGRID2B values is zero.\n");
  }

  s=mius.size();

  auto init_atoms = get_init_atoms(c);
  init_for_atoms(init_atoms);
}

void D2_Blip::calc_aed(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  tadah::core::aed_type &aed,
  const double scale)
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  size_t i=fidx;
  double fc_ij = fcut->calc(rij);
  for (size_t c=0; c<mius.size(); c++) {
    double t = B(etas[c]*(rij-mius[c]),fc_ij);
    aed(i++) += scale*weights[Zj]*t;
  }

}
void D2_Blip::calc_dXijdri(
  const int Zi,
  const int Zj,
  const double rij,
  const double ,
  tadah::core::fd_type &fd_ij,
  const double scale)
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  size_t i=fidx;
  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);
  for (size_t c=0; c<mius.size(); c++) {
    fd_ij(i++,0) = scale*weights[Zj]*dB(etas[c]*(rij-mius[c]),etas[c],fc_ij,fcp_ij);
  }
}
void D2_Blip::calc_all(
  const int Zi,
  const int Zj,
  const double rij,
  const double ,
  tadah::core::aed_type &aed,
  tadah::core::fd_type &fd_ij,
  const double scale)
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;

  size_t i=fidx;
  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);
  for (size_t c=0; c<mius.size(); c++) {
    aed(i) += scale*weights[Zj]*B(etas[c]*(rij-mius[c]),fc_ij);
    fd_ij(i,0) = scale*weights[Zj]*dB(etas[c]*(rij-mius[c]),etas[c],fc_ij,fcp_ij);
    ++i;
  }
}

std::string D2_Blip::label() {
  return lab;
}
void D2_Blip::init() {
  keys.push_back("CGRID2B");
  keys.push_back("SGRID2B");
  nparams=2;
}

} // namespace models
} // namespace tadah

