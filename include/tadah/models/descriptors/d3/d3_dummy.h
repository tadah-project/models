#ifndef D3_DUMMY_H
#define D3_DUMMY_H
#include <tadah/models/descriptors/d3/d3_base.h>

namespace tadah {
namespace models {
class D3_Dummy: public D3_Base {

private:
  size_t s=0;
  std::string lab="D3_Dummy";
public:
  D3_Dummy();
  D3_Dummy(tadah::core::Config&);
  void calc_aed(
    const size_t fidx,
    const double rij,
    const double rik,
    const double fc_ij,
    const double fc_ik,
    tadah::core::aed_type& aed);
  void calc_fd(
    const size_t fidx,
    const double rij,
    const double rik,
    const double fc_ij,
    const double fc_ik,
    const double fcp_ij,
    const double fcp_ik,
    tadah::core::fd_type &fd_ij);
  void calc_all(
    const size_t fidx,
    const double rij,
    const double rik,
    const double fc_ij,
    const double fc_ik,
    const double fcp_ij,
    const double fcp_ik,
    tadah::core::aed_type& aed,
    tadah::core::fd_type &fd_ij);
  size_t size();

  std::string label();
  void init() override;
};
} // namespace models
} // namespace tadah
#endif
