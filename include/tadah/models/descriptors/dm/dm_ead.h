#ifndef DM_EAD_H
#define DM_EAD_H
#include <tadah/models/descriptors/dm/dm_base.h>

namespace tadah {
namespace models {

/**  Embedded Atom Descriptor
 *
 * \f[
 * V_i^{L,\eta,r_s} = \sum_{l_x,l_y,l_z}^{l_x+l_y+l_z=L} \frac{L!}{l_x!l_y!l_z!}
 * \Big( \rho_i^{\eta,r_s,l_x,l_y,l_z} \Big)^2
 * \f]
 *
 * where density \f$ \rho \f$ is calculated using Gaussian Type Orbitals:
 * \f[
 * \rho_i^{\eta,r_s,l_x,l_y,l_z} = \sum_{j \neq i} x_{ij}^{l_x}y_{ij}^{l_y}z_{ij}^{l_z}
 * \exp{\Big(-\eta(r_{ij}-r_s)^2\Big)}f_c(r_{ij})
 * \f]
 *
 * \ref CGRIDMB parameters control position \f$ r_s \f$ of the gaussian basis function.
 *
 * \ref SGRIDMB parameters control width \f$ \eta \f$ of the gaussian basis function.
 *
 * e.g. \f$L_{max}=2\f$ will calculate descriptors with \f$ L=0,1,2 \f$ (s,p,d orbitals).
 *
 * More information about this descriptor:
 *
 * <div class="csl-entry">Zhang, Y., Hu, C.,Jiang, B.
 * (2019). Embedded atom neural network potentials: efficient and accurate
 * machine learning with a physically inspired representation.
 * <i>Journal of Physical Chemistry Letters</i>, <i>10</i>(17),
 * 4962–4967. https://doi.org/10.1021/acs.jpclett.9b02037</div>
 *
 * Required tadah::core::Config keys:
 * \ref INITMB \ref CGRIDMB \ref SGRIDMB
 */
class DM_EAD: public DM_Base {
private:
  int Lmax=0;
  std::string lab="DM_EAD";
  tadah::core::v_type sgrid;
  tadah::core::v_type cgrid;
  std::vector<std::vector<std::vector<int>>> orbitals;
  std::vector<std::vector<double>> fac;
  size_t gen_atomic_orbitals(int L);
  double fact(int n);
  double my_pow(double,int);
  double rc;

public:
  DM_EAD();
  DM_EAD(tadah::core::Config&);
  void calc_aed(
    tadah::core::rho_type& rho,
    tadah::core::aed_type &aed) override;
  void calc_dXijdri_dXjidri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d &vec_ij,
    tadah::core::rho_type& rhoi,
    tadah::core::rho_type& rhoj,
    tadah::core::fd_type &fd_ij,
    const double scale=1) override;
  void calc_dXijdri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d &vec_ij,
    tadah::core::rho_type& rhoi,
    tadah::core::fd_type &fd_ij, 
    const double scale=1) override;
  std::string label() override;
  void init_rhoi(tadah::core::rho_type& rhoi) override;
  void calc_rho(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d &vec_ij,
    tadah::core::rho_type& rho,

    const double scale=1) override;
  void init() override;
};
} // namespace models
} // namespace tadah
#endif
