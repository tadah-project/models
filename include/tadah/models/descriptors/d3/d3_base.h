#ifndef D3_BASE_H
#define D3_BASE_H

#include <tadah/models/descriptors/d_base.h>

namespace tadah {
namespace models {

class D3_Base: public D_Base {
public:
  virtual ~D3_Base() {};

  virtual void calc_aed(
    const size_t fidx,
    const double rij,
    const double rik,
    const double fc_ij,
    const double fc_ik,
    tadah::core::aed_type& aed)=0;
  virtual void calc_fd(
    const size_t fidx,
    const double rij,
    const double rik,
    const double fc_ij,
    const double fc_ik,
    const double fcp_ij,
    const double fcp_ik,
    tadah::core::fd_type &fd_ij)=0;
  virtual void calc_all(
    const size_t fidx,
    const double rij,
    const double rik,
    const double fc_ij,
    const double fc_ik,
    const double fcp_ij,
    const double fcp_ik,
    tadah::core::aed_type& aed,
    tadah::core::fd_type &fd_ij)=0;
  virtual size_t size()override=0;

  virtual std::string label()override=0;
  virtual std::vector<std::string> get_init_atoms(tadah::core::Context &c) override;
};
} // namespace models
} // namespace tadah
#endif
