#include <tadah/models/functions/basis_functions/bf_polynomial2.h>
namespace tadah {
namespace models {


BF_Polynomial2::BF_Polynomial2() {}
BF_Polynomial2::BF_Polynomial2(const tadah::core::Context &c): 
  Function_Base(c),
  BF_Base(c)
{}
std::string BF_Polynomial2::get_label() const
{
  return label;
}
std::size_t BF_Polynomial2::weights_count(std::size_t dsize) const
{
  return dsize*(dsize+1)/2;
}
double BF_Polynomial2::epredict(const tadah::core::t_type &weights, const tadah::core::aed_type& aed) const
{
  size_t b=0;
  double res=0.0;
  for (size_t i=0; i<aed.size(); ++i) {
    for (size_t ii=i; ii<aed.size(); ++ii) {
      res+= weights(b++)*aed(i)*aed(ii);
    }
  }
  return res;
}
double BF_Polynomial2::fpredict(const tadah::core::t_type &weights, const tadah::core::fd_type &fdij,
                                const tadah::core::aed_type& aedi, const size_t k) const
{
  double res=0.0;
  size_t b=0;
  for (size_t i=0; i<fdij.rows(); ++i) {
    for (size_t ii=i; ii<fdij.rows(); ++ii) {
      res -= weights(b++)*(fdij(i,k)*aedi(ii) + fdij(ii,k)*aedi(i));
    }
  }
  return res;
}
tadah::core::force_type BF_Polynomial2::fpredict(const tadah::core::t_type &weights, const tadah::core::fd_type &fdij,
                                    const tadah::core::aed_type& aedi) const
{
  tadah::core::force_type v;
  for (size_t k=0; k<3; ++k) {
    size_t b=0;
    for (size_t i=0; i<fdij.rows(); ++i) {
      for (size_t ii=i; ii<fdij.rows(); ++ii) {
        v(k) -= weights(b++)*(fdij(i,k)*aedi(ii) + fdij(ii,k)*aedi(i));
      }
    }
  }
  return v;
}

} // namespace models
} // namespace tadah

