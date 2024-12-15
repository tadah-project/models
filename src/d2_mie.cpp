#include "tadah/models/descriptors/d2/d2_base.h"
#include <tadah/models/descriptors/d2/d2_mie.h>

D2_MIE::D2_MIE() {}
D2_MIE::D2_MIE(Config &c): D2_Base(c)
{
  if (!c.get<bool>("INIT2B")) return;
  s=2;
  get_grid(c,"SGRID2B",etas);
  if (etas.size()!=2) {
    throw std::runtime_error("Number of elements in SGRID2B is != 2\n\
Mie descriptor requires SGRID2B with two positive elements.\n");
  }
  n=etas[1];
  m=etas[0];
  if (n<0 || m<0) {
    throw std::runtime_error("Both Mie exponents must by positive\n");
  }
}

void D2_MIE::calc_aed(
  const int Zj,
  const double rij,
  const double,
  aed_type2 &aed,
  const double scale)
{
  double rn_inv = 1.0/pow(rij,n);
  double rm_inv = 1.0/pow(rij,m);
  double fc_ij = fcut->calc(rij);

  aed(fidx) -= scale*weights[Zj]*rn_inv*fc_ij;
  aed(fidx+1) += scale*weights[Zj]*rm_inv*fc_ij;
}
void D2_MIE::calc_dXijdri(
  const int Zj,
  const double rij,
  const double,
  fd_type &fd_ij,
  const double scale)
{
  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);

  fd_ij(fidx,0)   = scale*weights[Zj]*(n/pow(rij,n+1)*fc_ij - fcp_ij/pow(rij,n));
  fd_ij(fidx+1,0) = scale*weights[Zj]*(-m/pow(rij,m+1)*fc_ij + fcp_ij/pow(rij,m));
}
void D2_MIE::calc_all(
  const int Zj,
  const double rij,
  const double,
  aed_type2 &aed,
  fd_type &fd_ij,
  const double scale)
{

  double rn_inv = 1.0/pow(rij,n);
  double rm_inv = 1.0/pow(rij,m);

  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);

  aed(fidx) -= scale*weights[Zj]*rn_inv*fc_ij;
  aed(fidx+1) += scale*weights[Zj]*rm_inv*fc_ij;

  fd_ij(fidx,0) = scale*weights[Zj]*(-n*rn_inv/rij*fc_ij - fcp_ij*rn_inv);
  fd_ij(fidx+1,0) = scale*weights[Zj]*(-m*rm_inv/rij*fc_ij + fcp_ij*rm_inv);

}
std::string D2_MIE::label() {
  return lab;
}
