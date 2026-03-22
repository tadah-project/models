#ifndef DM_POLY_H
#define DM_POLY_H
#include <tadah/models/descriptors/dm/dm_base.h>

namespace tadah {
namespace models {

/** \brief Power-Spectrum Many-Body Descriptor (Blip radial basis).
 *
 * Generalisation of DM_Blip in which a channel function
 * \f$\mathcal{F}\f$ is applied to the rotationally invariant power-spectrum
 * scalar \f$\varrho_i^{L,c}\f$ rather than to individual orbital densities.
 * For any choice of \f$\mathcal{F}\f$ the descriptor is rotationally invariant.
 *
 * \f[
 *   V_i^{L,c} = \mathcal{F}\!\left(
 *     \sum_{\substack{l_x,l_y,l_z \ge 0 \\ l_x+l_y+l_z=L}}
 *     \frac{L!}{l_x!\,l_y!\,l_z!}\,
 *     \Bigl(\rho_i^{L,c,l_x,l_y,l_z}\Bigr)^2
 *   \right)
 * \f]
 *
 * where the projected density is built from blip (cubic B-spline) basis
 * functions:
 *
 * \f[
 *   \rho_i^{L,c,l_x,l_y,l_z}
 *   = \sum_{j \neq i} x_{ij}^{l_x}\,y_{ij}^{l_y}\,z_{ij}^{l_z}\,
 *     \mathcal{B}\!\Bigl(\eta_c(r_{ij}-r_{s,c})\Bigr)\,f_c(r_{ij})
 * \f]
 *
 * The multinomial prefactor is normalised by \f$(4r_c)^L\f$.
 *
 * **Compare to DM_Blip / DM_mEAD:**
 * - DM_Blip: \f$\mathcal{F}=\mathrm{id}\f$ applied after the quadratic sum
 *   (i.e.\ \f$V = \varrho_L\f$, the raw power spectrum).
 * - DM_mEAD: \f$\mathcal{F}\f$ applied *before* the quadratic sum
 *   (orbital-wise); only rotationally invariant when \f$\mathcal{F}=\rho^2\f$.
 * - DM_Poly: \f$\mathcal{F}\f$ applied *after* the quadratic sum;
 *   invariant for any \f$\mathcal{F}\f$.
 *
 * Available instantiations registered in dm_all.cpp:
 * | Config name    | ChannelFn  | Description                              |
 * |----------------|------------|------------------------------------------|
 * | `DM_PS_SQ`     | F_SQ       | \f$(\varrho_L)^2\f$ — 4th-order term    |
 * | `DM_PS_SQRT`   | F_SQRT     | \f$s\sqrt{\varrho_L}\f$ — Finnis-Sinclair-like |
 * | `DM_PS_RLR`    | F_RLR      | \f$s\,\varrho_L\ln(c\,\varrho_L)\f$    |
 * | `DM_PS_Blip`   | F_Blip     | Blip-basis expansion of \f$\mathcal{F}\f$ |
 *
 * Required Config keys:
 * \ref INITMB \ref CGRIDMB \ref SGRIDMB (plus channel-function keys)
 */
template <typename ChannelFn>
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
