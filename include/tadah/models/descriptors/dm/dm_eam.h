#ifndef DM_EAM_H
#define DM_EAM_H
#include <tadah/models/descriptors/dm/dm_base.h>
#include <tadah/models/eam/eam_setfl_io.h>

namespace tadah {
namespace models {

/** \brief many-body part for the Embedded Atom Method descriptor.
 *
 * \f[
 * V_i = F\Bigg(\sum_{j \neq i} \rho(r_{ij}) \Bigg)
 * \f]
 *
 * This descriptor will load tabulated values for the density \f$ \rho \f$
 * and embedded energy \f$ F \f$ from the provided \ref SETFL file.
 *
 * This descriptor is usually used together with the two-body descriptor \ref D2_EAM
 * although this is not required and user can mix it with any other descriptors
 * or use it on its own.
 *
 * This descriptor will enforce cutoff distance as specified in a \ref SETFL file.
 * Set \ref RCUTMB to the same value to suppress the warning message.
 *
 * The setfl parsing, cubic splines and extrapolation are provided by the shared
 * \ref EAMSetflIO class (also used by \ref D2_EAM and the `tadah refit` machinery).
 *
 * Required tadah::core::Context keys:
 * \ref INITMB \ref SETFL
 */
class DM_EAM: public DM_Base {
private:
  EAMSetflIO eam_;
  std::string lab="DM_EAM";

public:
  DM_EAM();
  DM_EAM(tadah::core::Context&);
  DM_EAM(const tadah::core::Context&);
  void calc_aed(
    tadah::core::rho_type& rho,
    tadah::core::aed_type &aed) const override;
  void calc_dXijdri_dXjidri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d &vec_ij,
    tadah::core::rho_type& rhoi,
    tadah::core::rho_type& rhoj,
    tadah::core::fd_type &fd_ij,
    const double scale=1) const override;
  void calc_dXijdri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d &vec_ij,
    tadah::core::rho_type& rhoi,
    tadah::core::fd_type &fd_ij,
    const double scale=1) const override;
  std::string label() const override;
  void init_rhoi(tadah::core::rho_type& rhoi) override;
  void calc_rho(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d &vec_ij,
    tadah::core::rho_type& rho,
    const double scale=1) const override;
  void init() override;
  std::optional<DescriptorMeta> meta() const override;
};
} // namespace models
} // namespace tadah
#endif
