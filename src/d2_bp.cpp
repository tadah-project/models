#include <tadah/models/descriptors/d2/d2_bp.h>
#include <tadah/models/descriptors/d_basis_functions.h>

D2_BP::D2_BP() {
  init();
}
D2_BP::D2_BP(Config &c): D2_Base(c) {
  init();

  if (!c.get<bool>("INIT2B")) return;

  get_grid(c,"CGRID2B",mius);
  get_grid(c,"SGRID2B",etas);
  // update config with generated grid
  c.remove("CGRID2B");
  for (const auto &m: mius) c.add("CGRID2B",m);
  c.remove("SGRID2B");
  for (const auto &e: etas) c.add("SGRID2B",e);

  if (verbose) {
    std::cout << std::endl;
    std::cout << "SGRID2B: " << etas.size() << std::endl;
    for (auto e:etas) std::cout << e << "  ";
    std::cout << std::endl;

    std::cout << "CGRID2B: " << mius.size() << std::endl;
    for (auto m:mius) std::cout << m << "  ";
    std::cout << std::endl;
  }

  if (mius.size()!=etas.size()) {
    throw std::runtime_error("SGRID2B and CGRID2B arrays differ in size.\n");
  }

  for (const auto &e:etas) {
    if(std::abs(e)<std::numeric_limits<double>::min())
      throw std::runtime_error("At least one of SGRID2B values is zero.\n");
  }

  s=mius.size();
}

void D2_BP::calc_aed(
  const int Zj,
  const double rij,
  const double ,
  aed_type2 &aed,
  const double scale)
{

  size_t i=fidx;
  double fc_ij = fcut->calc(rij);
  for (size_t c=0; c<mius.size(); c++) {         
    aed(i++) += scale*weights[Zj]*G(rij,etas[c],mius[c],fc_ij);
  }

}
void D2_BP::calc_dXijdri(
  const int Zj,
  const double rij,
  const double ,
  fd_type &fd_ij,
  const double scale)
{
  size_t i=fidx;
  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);
  for (size_t c=0; c<mius.size(); c++) {         
    fd_ij(i++,0) = scale*weights[Zj]*dG(rij,etas[c],mius[c],fc_ij,fcp_ij);
  }
}
void D2_BP::calc_all(
  const int Zj,
  const double rij,
  const double ,
  aed_type2 &aed,
  fd_type &fd_ij,
  const double scale)
{
  size_t i=fidx;
  double fc_ij = fcut->calc(rij);
  double fcp_ij = fcut->calc_prime(rij);
  for (size_t c=0; c<mius.size(); c++) {         
    aed(i) += scale*weights[Zj]*G(rij,etas[c],mius[c],fc_ij);
    fd_ij(i,0) = scale*weights[Zj]*dG(rij,etas[c],mius[c],fc_ij,fcp_ij);
    ++i;
  }
}
std::string D2_BP::label() {
  return lab;
}
void D2_BP::init() {
  keys.push_back("CGRID2B");
  keys.push_back("SGRID2B");
  nparams=2;
}
