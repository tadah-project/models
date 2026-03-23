#ifndef DM_mEAD_H
#define DM_mEAD_H
#include <tadah/models/descriptors/dm/dm_base.h>

namespace tadah {
namespace models {

/** Modified Embedded Atom Descriptor
 *
 * \warning **Rotational invariance** holds only when \f$\mathcal{F}=\rho^2\f$
 * (i.e.\ F_SQ).  For any other channel function and \f$L>0\f$ the descriptor
 * is **not** rotationally invariant because \f$\mathcal{F}\f$ is applied to
 * individual orbital densities before the invariance-restoring quadratic sum.
 * A runtime warning is emitted in that case.  Use DM_mBlipL for a descriptor
 * that is invariant for arbitrary \f$\mathcal{F}_L\f$.
 *
 * REQUIRED KEYS: SGRIDMB, CGRIDMB, and KEYS OF THE EMBEDDING FUNCTION
 *
 * This descriptor has a mathematical form very similar to DM_EAD
 * but allows the usage of a custom-defined embedding function, \f$ \mathcal{F} \f$, in place of the default quadratic one.
 *
 * Available implementations: 
 *   1. DM_mRLR uses F_RLR
 *   2. DM_mSQRT uses F_SQRT
 *
 * \f[
 * V_i^{L,\eta,r_s} = \sum_{l_x,l_y,l_z}^{l_x+l_y+l_z=L} \frac{L!}{l_x!l_y!l_z!}
 * \mathcal{F}\Big( \rho_i^{\eta,r_s,l_x,l_y,l_z} \Big)
 * \f]
 *
 * where the density \f$ \rho \f$ is calculated using Gaussian Type Orbitals:
 * \f[
 * \rho_i^{\eta,r_s,l_x,l_y,l_z} = \sum_{j \neq i} x_{ij}^{l_x} y_{ij}^{l_y} z_{ij}^{l_z}
 * \exp{\Big(-\eta(r_{ij}-r_s)^2\Big)} f_c(r_{ij})
 * \f]
 *
 * \ref CGRIDMB parameters control the position \f$ r_s \f$ of the Gaussian basis function.
 *
 * \ref SGRIDMB parameters control the width \f$ \eta \f$ of the Gaussian basis function.
 *
 * e.g., \f$L_{max}=2\f$ will calculate descriptors with \f$ L=0,1,2 \f$ (s, p, d orbitals).
 *
 * ```
 * # TYPEMB params: L, size(cgrid), size(sgrid),
 * # size(cembfunc), size(sembfunc), list of element pairs
 * TYPEMB    DM_mRLR 0 7 7 1 1 Ta Ta    
 * RCTYPEMB  Cut_Tanh
 * RCUTMB    7.5
 * SGRIDMB   -2 7 0.1 10
 * CGRIDMB   0 0 0 0 0 0 0
 * SEMBFUNC  1.2
 * CEMBFUNC  0.5
 * ```
 *
 * ```
 * # TYPEMB params: L, size(cgrid), size(sgrid),
 * # size(cembfunc), size(sembfunc), list of element pairs
 * TYPEMB    DM_mSQRT 1 5 5 0 1 * * 
 * RCTYPEMB  Cut_Cos
 * RCUTMB    3.0
 * CGRIDMB   -1 5 0 3.0
 * SGRIDMB   -2 5 1.0 10.0
 * SEMBFUNC   1.5
 * ```
 *
 * Required Config keys:
 * \ref INITMB \ref CGRIDMB \ref SGRIDMB
 */

template <typename F>
class DM_mEAD: public DM_Base {
private:
  int Lmax=0;
  std::string lab="DM_mEAD";
  tadah::core::v_type sgrid;
  tadah::core::v_type cgrid;
  std::vector<std::vector<std::vector<int>>> orbitals;
  std::vector<std::vector<double>> fac;
  size_t gen_atomic_orbitals(int L);
  double rc;
  F emb;

public:
  DM_mEAD();
  DM_mEAD(const tadah::core::Context&);
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
};
} // namespace models
} // namespace tadah

#include "dm_mEAD.hpp"

#endif
