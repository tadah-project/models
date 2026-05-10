#ifndef DM_BLIP_H
#define DM_BLIP_H
#include <tadah/models/descriptors/dm/dm_base.h>

namespace tadah {
namespace models {

/** \brief Blip Many-Body Descriptor.
 *
 * Angular-momentum-projected descriptor built from blip (cubic B-spline) basis
 * functions.  For each angular-momentum channel \f$ L \f$ and each grid point
 * \f$(\eta,r_s)\f$ the descriptor component is:
 *
 * \f[
 *   V_i^{L,\eta,r_s} = \sum_{\substack{l_x,l_y,l_z \ge 0 \\ l_x+l_y+l_z=L}}
 *   \frac{L!}{l_x!\,l_y!\,l_z!}
 *   \Bigl(\rho_i^{\eta,r_s,l_x,l_y,l_z}\Bigr)^2
 * \f]
 *
 * where the projected density is
 *
 * \f[
 *   \rho_i^{\eta,r_s,l_x,l_y,l_z}
 *   = \sum_{j \neq i} x_{ij}^{l_x}\,y_{ij}^{l_y}\,z_{ij}^{l_z}\,
 *     \mathcal{B}\!\Bigl(\eta(r_{ij}-r_s)\Bigr)\,f_c(r_{ij})
 * \f]
 *
 * Here \f$ \mathcal{B} \f$ is the blip (cubic B-spline) basis function with
 * compact support \f$ |\eta(r-r_s)| < 2 \f$ (see tadah::models::B), and
 * \f$ f_c \f$ is the cutoff function.
 *
 * The multinomial prefactor \f$ L!/({l_x!\,l_y!\,l_z!}) \f$ is normalised by
 * \f$(4r_c)^L\f$ so that descriptor components at different \f$ L \f$ have
 * comparable magnitudes.
 *
 * \ref CGRIDMB parameters control the centre \f$ r_s \f$ of each blip basis
 * function.
 *
 * \ref SGRIDMB parameters control the inverse width \f$ \eta \f$; the full
 * support of each blip in \f$ r \f$-space is \f$ 4/\eta \f$.
 *
 * Setting \f$ L_{\max}=2 \f$ generates descriptors for \f$ L=0,1,2 \f$
 * (s-, p-, d-orbital channels). Maximum supported \f$ L_{\max} = 6 \f$.
 *
 * More information about this descriptor family:
 *
 * <div class="csl-entry">Zhang, Y., Hu, C., Jiang, B.
 * (2019). Embedded atom neural network potentials: efficient and accurate
 * machine learning with a physically inspired representation.
 * <i>Journal of Physical Chemistry Letters</i>, <i>10</i>(17),
 * 4962–4967. https://doi.org/10.1021/acs.jpclett.9b02037</div>
 *
 * Required tadah::core::Context keys:
 * \ref INITMB \ref CGRIDMB \ref SGRIDMB
 */
class DM_Blip: public DM_Base {
private:
  int Lmax=0;
  std::string lab="DM_Blip";
  tadah::core::v_type sgrid;
  tadah::core::v_type cgrid;
  std::vector<std::vector<std::vector<int>>> orbitals;
  std::vector<std::vector<double>> fac;
  size_t gen_atomic_orbitals(int L);
  double fact(int n);
  double my_pow(double,int);
  double rc;

public:
  DM_Blip();
  DM_Blip(const tadah::core::Context&);
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
  bool reinit(const tadah::core::Context& ctx) override;
  std::optional<DescriptorMeta> meta() const override;
};
} // namespace models
} // namespace tadah
#endif
