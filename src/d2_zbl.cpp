#include <tadah/models/descriptors/d2/d2_zbl.h>

D2_ZBL::D2_ZBL() {
  init();
}
D2_ZBL::D2_ZBL(Config &c): D2_Base(c)
{
  if (!c.get<bool>("INIT2B")) return;
  init();
  s=1;
  double temp;
  temp = c.get<double>("TYPE2B",1);
  if (temp!=-1) a0=temp;
  temp = c.get<double>("TYPE2B",2);
  if (temp!=-1) s0=temp;
  temp = c.get<double>("TYPE2B",3);
  if (temp!=-1) p0=temp;
  temp = c.get<double>("TYPE2B",4);
  if (temp!=-1) p1=temp;

  ntypes = c.size("ATOMS");
  types = new int[ntypes];
  a = new double*[ntypes];
  for (int i = 0; i < ntypes; ++i) {
    a[i] = new double[ntypes];
  }
  for (int i=0; i<ntypes; ++i) {
    std::string symbol = c.get<std::string>("ATOMS",i);
    int Zi = PeriodicTable::find_by_symbol(symbol).Z;
    types[i] = Zi;
  }
  for (int i=0; i<ntypes; ++i) {
    for (int j=0; j<ntypes; ++j) {
      a[i][j] = screening_length(types[i],types[j]);
    }
  }
}
D2_ZBL::~D2_ZBL() {
  if (a) {
    for (int i = 0; i < ntypes; ++i) {
      delete[] a[i];  
    }
    delete[] a;        
    delete[] types;        
    a = nullptr;       
    types = nullptr;       
  }
}

void D2_ZBL::calc_aed(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  aed_type &aed,
  const double scale)
{
  double fc_ij = fcut->calc(rij);
  double a_ = a[types[Zi]][types[Zj]];
  double x = rij / a_;
  double p = phi(x);
  aed(fidx) += scale*weights[Zj]*(Zi * Zj / rij) * p * fc_ij;
}
void D2_ZBL::calc_dXijdri(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  fd_type &fd_ij,
  const double scale)
{
  double r_inv = 1.0/rij;
  double a_ = a[types[Zi]][types[Zj]];
  double x = rij / a_;
  double p = phi(x);
  double dp = dphi(x);

  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);
  double E_ZBL = (Zi * Zj / rij) * p;
  double dE_ZBL_r = -((Zi * Zj *r_inv*r_inv)*p + (Zi * Zj * r_inv) / a_ * dp);

  fd_ij(fidx,0) = -scale*weights[Zj]*(E_ZBL*fcp_ij + dE_ZBL_r*fc_ij);
}
void D2_ZBL::calc_all(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  aed_type &aed,
  fd_type &fd_ij,
  const double scale)
{
  double r_inv = 1.0/rij;
  double a_ = a[types[Zi]][types[Zj]];
  double x = rij / a_;
  double p = phi(x);
  double dp = dphi(x);

  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);
  double E_ZBL = (Zi * Zj / rij) * p;
  double dE_ZBL_r = -((Zi * Zj *r_inv*r_inv)*p + (Zi * Zj * r_inv) / a_ * dp);

  aed(fidx) += scale*weights[Zj]*(Zi * Zj * r_inv) * p * fc_ij;
  fd_ij(fidx,0) = -scale*weights[Zj]*(E_ZBL*fcp_ij + dE_ZBL_r*fc_ij);
}
std::string D2_ZBL::label() {
  return lab;
}
void D2_ZBL::init() {
  nparams=4;
}
double D2_ZBL::phi(double x) {
  return 0.1818 * std::exp(-3.2 * x) +
  0.5099 * std::exp(-0.9423 * x) +
  0.2802 * std::exp(-0.4029 * x) +
  0.02817 * std::exp(-0.2016 * x);
}
double D2_ZBL::dphi(double x) {
  return -3.2 * 0.1818 * std::exp(-3.2 * x) -
  0.9423 * 0.5099 * std::exp(-0.9423 * x) -
  0.4029 * 0.2802 * std::exp(-0.4029 * x) -
  0.2016 * 0.02817 * std::exp(-0.2016 * x);
}
double D2_ZBL::screening_length(int Zi, int Zj) {
  return (s0 * a0) / (std::pow(Zi, p0) + std::pow(Zj, p1));
}
