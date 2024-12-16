#ifndef DM_DUMMY_H
#define DM_DUMMY_H
#include <tadah/models/descriptors/dm/dm_base.h>

/**
 * Dummy many-body descriptor.
 *
 * Use it to satisfy DescriptorsCalc requirements in case
 * when many-body descriptor is not required.
 *
 */
class DM_Dummy: public DM_Base {
private:
  std::string lab="DM_Dummy";
public:
  DM_Dummy();
  DM_Dummy(Config&);
  void calc_aed(
    rho_type& rho,
    aed_type2 &aed) override;
  void calc_dXijdri_dXjidri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const Vec3d &vec_ij,
    rho_type& rhoi,
    rho_type& rhoj,
    fd_type &fd_ij,
    const double scale=1) override;
  void calc_dXijdri(
    const int Zj,
    const double rij,
    const double rij_sq,
    const Vec3d &vec_ij,
    rho_type& rhoi,
    fd_type &fd_ij, 
    const double scale=1) override;
  std::string label() override;
  void init_rhoi(rho_type& rhoi) override;
  void calc_rho(
    const int Zj,
    const double rij,
    const double rij_sq,
    const Vec3d &vec_ij,
    rho_type& rho,
    const double scale=1) override;
  void init() override;
};
#endif
