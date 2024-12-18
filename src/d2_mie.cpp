#include "tadah/models/descriptors/d2/d2_base.h"
#include <tadah/models/descriptors/d2/d2_mie.h>

D2_MIE::D2_MIE() {
  init();
}
D2_MIE::D2_MIE(Config &c): D2_Base(c)
{
  if (!c.get<bool>("INIT2B")) return;
  init();
  s=2;
  n = c.get<double>("TYPE2B",1);
  m = c.get<double>("TYPE2B",2);
  if (n<0 || m<0) {
    throw std::runtime_error("Both Mie exponents must by positive\n");
  }
}

void D2_MIE::calc_aed(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  aed_type &aed,
  const double scale)
{
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
  fd_type &fd_ij,
  const double scale)
{
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
  aed_type &aed,
  fd_type &fd_ij,
  const double scale)
{

  double rn_inv = 1.0/pow(rij,n);
  double rm_inv = 1.0/pow(rij,m);

  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);

  aed(fidx) -= scale*weights[Zj]*rn_inv*fc_ij;
  aed(fidx+1) += scale*weights[Zj]*rm_inv*fc_ij;

  fd_ij(fidx,0) = scale*weights[Zj]*(n*rn_inv/rij*fc_ij - fcp_ij*rn_inv);
  fd_ij(fidx+1,0) = scale*weights[Zj]*(-m*rm_inv/rij*fc_ij + fcp_ij*rm_inv);

}
std::string D2_MIE::label() {
  return lab;
}
void D2_MIE::init() {
  nparams=2;
}
