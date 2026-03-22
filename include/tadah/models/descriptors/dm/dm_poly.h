#ifndef DM_POLY_H
#define DM_POLY_H
#include <tadah/models/descriptors/dm/dm_base.h>
#include <tadah/models/descriptors/dm/radial_policies.h>

namespace tadah {
namespace models {

/** \brief Power-Spectrum Many-Body Descriptor with swappable radial basis.
 *
 * Generalisation of DM_Blip in which:
 * 1. The radial basis is a compile-time policy (`RadialBasis`).
 * 2. A channel function \f$\mathcal{F}\f$ is applied to the rotationally
 *    invariant power-spectrum scalar \f$\varrho_i^{L,c}\f$.
 *
 * For any choice of \f$\mathcal{F}\f$ **and** any radial basis the descriptor
 * is rotationally invariant.
 *
 * \f[
 *   V_i^{L,c} = \mathcal{F}\!\left(
 *     \sum_{\substack{l_x,l_y,l_z \ge 0 \\ l_x+l_y+l_z=L}}
 *     \frac{L!}{l_x!\,l_y!\,l_z!}\,
 *     \Bigl(\rho_i^{L,c,l_x,l_y,l_z}\Bigr)^2
 *   \right)
 * \f]
 *
 * where the projected density is
 *
 * \f[
 *   \rho_i^{L,c,l_x,l_y,l_z}
 *   = \sum_{j \neq i} x_{ij}^{l_x}\,y_{ij}^{l_y}\,z_{ij}^{l_z}\,
 *     \phi_c(r_{ij})\,f_c(r_{ij})
 * \f]
 *
 * with \f$\phi_c\f$ provided by the `RadialBasis` policy.
 *
 * The multinomial prefactor is normalised by \f$(4r_c)^L\f$.
 *
 * **Compare to DM_Blip / DM_mEAD:**
 * - DM_Blip: \f$\mathcal{F}=\mathrm{id}\f$, BlipRadial — raw power spectrum.
 * - DM_mEAD: \f$\mathcal{F}\f$ applied *before* the quadratic sum
 *   (orbital-wise); only rotationally invariant when \f$\mathcal{F}=\rho^2\f$.
 * - DM_Poly: \f$\mathcal{F}\f$ applied *after* the quadratic sum;
 *   invariant for any \f$\mathcal{F}\f$ and any `RadialBasis`.
 *
 * \tparam ChannelFn  Channel function (F_SQ, F_SQRT, F_RLR, F_Blip, …).
 * \tparam RadialBasis  Radial basis policy (BlipRadial, GaussianRadial, …).
 *                      Default: BlipRadial.
 *
 * Available instantiations registered in dm_all.cpp:
 * | Config name         | ChannelFn | RadialBasis      |
 * |---------------------|-----------|------------------|
 * | `DM_PS_SQ`          | F_SQ      | BlipRadial       |
 * | `DM_PS_SQRT`        | F_SQRT    | BlipRadial       |
 * | `DM_PS_RLR`         | F_RLR     | BlipRadial       |
 * | `DM_PS_Blip`        | F_Blip    | BlipRadial       |
 * | `DM_GPS_SQ`         | F_SQ      | GaussianRadial   |
 * | `DM_GPS_SQRT`       | F_SQRT    | GaussianRadial   |
 * | `DM_GPS_RLR`        | F_RLR     | GaussianRadial   |
 * | `DM_GPS_Blip`       | F_Blip    | GaussianRadial   |
 *
 * Required Config keys:
 * \ref INITMB \ref CGRIDMB \ref SGRIDMB (plus channel-function keys)
 */
template <typename ChannelFn, typename RadialBasis = BlipRadial>
class DM_Poly : public DM_Base {
private:
  int Lmax = 0;
  std::string lab = "DM_Poly";
  tadah::core::v_type sgrid;
  tadah::core::v_type cgrid;
  std::vector<std::vector<std::vector<int>>> orbitals;
  std::vector<std::vector<double>> fac;
  size_t gen_atomic_orbitals(int L);
  double rc;
  ChannelFn emb;

public:
  DM_Poly();
  DM_Poly(const tadah::core::Context&);
  void calc_aed(
    tadah::core::rho_type& rho,
    tadah::core::aed_type& aed) const override;
  void calc_dXijdri_dXjidri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d& vec_ij,
    tadah::core::rho_type& rhoi,
    tadah::core::rho_type& rhoj,
    tadah::core::fd_type& fd_ij,
    const double scale = 1) const override;
  void calc_dXijdri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d& vec_ij,
    tadah::core::rho_type& rhoi,
    tadah::core::fd_type& fd_ij,
    const double scale = 1) const override;
  std::string label() const override;
  void init_rhoi(tadah::core::rho_type& rhoi) override;
  void calc_rho(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d& vec_ij,
    tadah::core::rho_type& rho,
    const double scale = 1) const override;
  void init() override;
};

} // namespace models
} // namespace tadah

#include "dm_poly.hpp"

#endif
