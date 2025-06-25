#include <tadah/models/descriptors/d2/d2_zbl.h>
namespace tadah {
namespace models {


D2_ZBL::D2_ZBL() {
  init();
}
D2_ZBL::D2_ZBL(const tadah::core::Context &c): D2_Base(c)
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
  
  if (verbose) 
    std::cout << label() << ": " << a0 << " "
      << s0 << " " << p0 << " " << p1 << std::endl;

  ntypes = c.size("ATOM");
  a = new double*[ntypes];
  for (int i = 0; i < ntypes; ++i) {
    a[i] = new double[ntypes];
  }
  for (int i=0; i<ntypes; ++i) {
    std::string symboli = c.get<std::string>("ATOM",i);
    int Zi = tadah::core::PeriodicTable::find_by_symbol(symboli).Z;
    types[Zi] = i;
    for (int j=0; j<ntypes; ++j) {
      std::string symbolj = c.get<std::string>("ATOM",j);
      int Zj = tadah::core::PeriodicTable::find_by_symbol(symbolj).Z;
      a[i][j] = screening_length(Zi,Zj);
    }
  }
  auto init_atoms = get_init_atoms(c);
  init_for_atoms(init_atoms);
}
D2_ZBL::~D2_ZBL() {
  if (a) {
    for (int i = 0; i < ntypes; ++i) {
      delete[] a[i];  
    }
    delete[] a;        
    a = nullptr;       
  }
}

void D2_ZBL::calc_aed(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  tadah::core::aed_type &aed,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
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
  tadah::core::fd_type &fd_ij,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  double r_inv = 1.0/rij;
  double a_ = a[types[Zi]][types[Zj]];
  double x = rij / a_;
  double p = phi(x);
  double dp = dphi(x);
  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);

  double E_ZBL = (Zi * Zj * r_inv) * p;
  double dE_ZBL_r = -(Zi * Zj * r_inv*r_inv)*p + dp*(Zi * Zj * r_inv) / a_ ;
  fd_ij(fidx,0) = scale*weights[Zj]*(E_ZBL*fcp_ij + dE_ZBL_r*fc_ij);
}
void D2_ZBL::calc_all(
  const int Zi,
  const int Zj,
  const double rij,
  const double,
  tadah::core::aed_type &aed,
  tadah::core::fd_type &fd_ij,
  const double scale) const
{
  if (!is_init_for_atoms(Zi,Zj) || rij > get_rcut()) return;
  double r_inv = 1.0/rij;
  double a_ = a[types[Zi]][types[Zj]];
  double x = rij / a_;
  double p = phi(x);
  double dp = dphi(x);

  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);

  aed(fidx) += scale*weights[Zj]*(Zi * Zj * r_inv) * p * fc_ij;
  double E_ZBL = (Zi * Zj * r_inv) * p;
  double dE_ZBL_r = -(Zi * Zj * r_inv*r_inv)*p + dp*(Zi * Zj * r_inv) / a_ ;
  fd_ij(fidx,0) = scale*weights[Zj]*(E_ZBL*fcp_ij + dE_ZBL_r*fc_ij);
}
std::string D2_ZBL::label() const {
  return lab;
}
void D2_ZBL::init() {
  nparams=4;
}
double D2_ZBL::phi(double x) const {
  return 0.1818 * std::exp(-3.2 * x) +
  0.5099 * std::exp(-0.9423 * x) +
  0.2802 * std::exp(-0.4029 * x) +
  0.02817 * std::exp(-0.2016 * x);
}
double D2_ZBL::dphi(double x) const {
  return -3.2 * 0.1818 * std::exp(-3.2 * x) -
  0.9423 * 0.5099 * std::exp(-0.9423 * x) -
  0.4029 * 0.2802 * std::exp(-0.4029 * x) -
  0.2016 * 0.02817 * std::exp(-0.2016 * x);
}
double D2_ZBL::screening_length(int Zi, int Zj) const {
  return (s0 * a0) / (std::pow(Zi, p0) + std::pow(Zj, p1));
}

} // namespace models
} // namespace tadah

