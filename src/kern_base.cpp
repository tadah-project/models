#include <tadah/models/functions/kernels/kern_base.h>
#include <tadah/core/intercept_layout.h>

#include <stdexcept>
#include <string>
#include <vector>
namespace tadah {
namespace models {


Kern_Base::Kern_Base() {}
Kern_Base::Kern_Base(const tadah::core::Context &c): Function_Base(c) {}
Kern_Base::~Kern_Base() {}
void Kern_Base::set_basis(const tadah::core::Matrix b)
{
  basis = b;
}
double Kern_Base::epredict(const tadah::core::t_type &kweights, const tadah::core::aed_type &aed) const
{
  double energy=0;
  for (size_t b=0; b<basis.cols(); ++b) {
    energy += kweights[b]*(*this)(basis.col(b),aed);
  }
  return energy;
}
double Kern_Base::fpredict(const tadah::core::t_type &kweights, const tadah::core::fd_type &fdij,
                           const tadah::core::aed_type &aedi, const size_t k) const
{
  double res=0.0;
  for (size_t b=0; b<basis.cols(); ++b) {
    res -= kweights[b]*(*this).prime(basis.col(b),aedi,fdij(k));
  }
  return res;
}
tadah::core::force_type Kern_Base::fpredict(const tadah::core::t_type &kweights, const tadah::core::fd_type &fdij,
                               const tadah::core::aed_type &aedi) const
{
  tadah::core::force_type v;
  v.set_zero();
  for (size_t b=0; b<basis.cols(); ++b) {
    for (size_t k=0; k<3; ++k) {
      v(k) -= kweights[b]*(*this).prime(basis.col(b),aedi,fdij(k));
    }
  }
  return v;
}
void Kern_Base::read_basis_from_context(const tadah::core::Context &context, tadah::core::Matrix &basis) {
  // storage is in column-major order
  if (!context.contains("BASIS")) return;
  if (!context.contains("SBASIS"))
    throw std::runtime_error(
        "BASIS without SBASIS: SBASIS names the number of basis vectors "
        "stored in BASIS.");
  const size_t sbasis = context.get<size_t>("SBASIS");
  const size_t nvalues = context.size("BASIS");
  if (sbasis == 0)
    throw std::runtime_error(
        "SBASIS 0 with " + std::to_string(nvalues) + " BASIS values: "
        "SBASIS must be the number of basis vectors stored in BASIS.");
  if (nvalues == 0 || nvalues % sbasis != 0)
    throw std::runtime_error(
        "BASIS has " + std::to_string(nvalues) + " values, which is not a "
        "multiple of SBASIS " + std::to_string(sbasis) + ": BASIS stores "
        "SBASIS basis vectors of the descriptor's length, one after another.");
  size_t vlen = nvalues/sbasis;
  std::vector<double> temp(sbasis*vlen);
  context.get<std::vector<double>>("BASIS",temp);
  basis.resize(vlen,sbasis);  // aeds are in column-major order
  size_t ii=0;
  for (size_t i=0;i<sbasis; ++i) {
    for (size_t j=0;j<vlen; ++j) {
      basis(j,i) = temp[ii++];
    }
  }
}
size_t Kern_Base::dot_product_intercept_columns(const tadah::core::Context &c) const {
  if (tadah::core::n_intercept_columns(c) == 0 || basis.cols() == 0) return 0;
  for (size_t j=1; j<basis.rows(); ++j)
    if (basis(j,0) != 0.0) return 0;
  return 1;
}
const tadah::core::Matrix &Kern_Base::get_basis() const {
  return basis;
}
tadah::core::Matrix Kern_Base::get_basis() {
  return basis;
}

} // namespace models
} // namespace tadah

