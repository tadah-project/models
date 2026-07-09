#ifndef DM_REAM_H
#define DM_REAM_H
#include <tadah/models/descriptors/dm/dm_base.h>

namespace tadah {
namespace models {

/** \brief Refittable EAM many-body descriptor (flexible density + embedding).
 *
 * The native-Tadah counterpart of \ref DM_EAM, used by `tadah refit`. Where
 * DM_EAM reproduces the *tabulated* density and embedding via fixed cubic
 * splines, DM_REAM represents BOTH as flexible basis-function expansions so an
 * EAM model can be reshaped and re-optimised with the standard linear model and
 * HPO machinery.
 *
 * The EAM many-body energy is
 * \f[
 *   V_i = F\!\Big(\rho_i\Big), \qquad \rho_i = \sum_{j\neq i} \rho(r_{ij}).
 * \f]
 *
 * Because \f$F(\rho)\f$ with \f$\rho=\sum_j\rho(r)\f$ is *bilinear*, only one of
 * the two functions can carry the linear weights. DM_REAM therefore:
 *
 *  - represents the **density** as a fixed (baked) expansion
 *    \f$\rho(r)=\sum_{n} a_n\,\psi_n(r)\,f_c(r)\f$ with coefficients
 *    \f$a_n\f$ supplied in \ref AMPGRIDMB (these are NOT linear weights — they are
 *    found by a 1-D least-squares pre-fit and stored in the descriptor), and
 *  - represents the **embedding** as
 *    \f$F(\rho)=\sum_{m} c_m\,g_m(\rho)\f$ where the \f$g_m\f$ are the descriptor
 *    output components (one per \ref CEMBFUNC entry) and the \f$c_m\f$ are the
 *    M_BLR **linear weights**.
 *
 * Keys (single element, \f$L=0\f$ scalar density):
 *  - \ref CGRIDMB / \ref SGRIDMB : density basis centres / widths (size Nρ)
 *  - \ref AMPGRIDMB                       : baked density coefficients \f$a_n\f$ (size Nρ)
 *  - \ref CEMBFUNC / \ref SEMBFUNC      : embedding basis centres / widths (size M)
 *
 * `TYPEMB` layout (the d_base.h keys contract):
 * `DM_REAM <basis> <nC> <nS> <nAmp> <nEc> <nEs> <nKd> <nKe> <El> <El>`
 * — one family param (the basis code) followed by one SIZE token per keys
 * entry ({CGRIDMB, SGRIDMB, AMPGRIDMB, CEMBFUNC, SEMBFUNC}) and per
 * data_keys entry ({KDEGMB, KDEGEMB}), which makes DM_REAM splittable by
 * DM_mJoin and HPO-addressable via get_arg_pos. `<basis>` is `0` (Blip,
 * default), `1` (Gaussian), `2` (Knot5: one-sided quintic/cubic knot
 * functions) or `3` (Blip5: symmetric quintic B-spline bumps, C^4) and
 * selects the family used for BOTH the density and the embedding
 * expansions. Knot pots carry SIZE 0 for the width channels (SGRIDMB /
 * SEMBFUNC ABSENT — eta is a pure gauge) and use the KDEG* data lists for
 * per-atom degrees (3|5); bump pots carry SIZE 0 for the KDEG* channels.
 * The amplitudes AMPGRIDMB are excluded from the LINEAR solve by
 * bilinearity but are first-class NONLINEAR HPO parameters
 * (`OPTIM AMPGRIDMB`). For Knot5 the per-channel orientation is part of the family
 * definition: the radial density basis is LEFT-active, \f$\psi_n(r) =
 * (\mu_n - r)_+^5\f$ (vanishes C^4-smoothly at its knot, so with all knots
 * \f$\le r_{cut}\f$ the density vanishes smoothly at the cutoff), while the
 * embedding basis is RIGHT-active, \f$g_m(\rho) = (\rho - \mu_m)_+^5\f$
 * (F free at high density, F(0) anchored by the fit). SGRID* entries are a
 * pure gauge for this family and are conventionally 1.0.
 *
 * \note Like \ref D2_Blip / \ref DM_Blip, the per-pair density carries the
 *       symmetric WATOM weight `weights[Zi][Zj]` (1.0 for a single-element
 *       potential), so multi-element generalisation and HPO WATOM search work.
 */
class DM_REAM: public DM_Base {
private:
  std::string lab="DM_REAM";
  // Runtime-selected radial/embedding basis family (from the TYPEMB <basis>
  // wire code). BasisKind is the shared enum; the per-basis math is reached
  // through the single dispatcher in basis_dispatch.h — no scattered
  // if(basis==...) ladders live in this class any more.
  BasisKind basis = BasisKind::Blip;
  // density expansion (baked)
  tadah::core::v_type crho;   // CGRIDMB
  tadah::core::v_type srho;   // SGRIDMB
  tadah::core::v_type arho;   // AMPGRIDMB  (fixed coefficients)
  // embedding expansion (linear weights live in the model)
  tadah::core::v_type cemb;   // CEMBFUNC
  tadah::core::v_type semb;   // SEMBFUNC
  // optional per-atom degrees for the knot family (bcode 2): 3|5, empty = 5
  tadah::core::v_type drho_deg;  // KDEGMB
  tadah::core::v_type demb_deg;  // KDEGEMB
  int rho_deg(std::size_t n) const {
    return n < drho_deg.size() ? static_cast<int>(drho_deg[n]) : 5;
  }
  int emb_deg(std::size_t m) const {
    return m < demb_deg.size() ? static_cast<int>(demb_deg[m]) : 5;
  }

  // ---- inline basis helpers (family selected by `basis`) ----
  /** Density radial value (cutoff-weighted) for centre n at distance r. */
  double psi(std::size_t n, double r, double fc) const;
  /** Density radial r-derivative (cutoff product rule) for centre n. */
  double dpsi(std::size_t n, double r, double fc, double fcp) const;
  /** Sum_n a_n psi_n(r) f_c(r). */
  double rho_raw(double r, double fc) const;
  /** d/dr [ Sum_n a_n psi_n(r) f_c(r) ]. */
  double rho_raw_deriv(double r, double fc, double fcp) const;
  /** Embedding value g_m(rho) (no cutoff). */
  double gemb(std::size_t m, double rho) const;
  /** Embedding derivative g'_m(rho) w.r.t. rho. */
  double dgemb(std::size_t m, double rho) const;
  /** Family-aware grid validation (ctor + reinit): size consistency, the
   *  knot family's width-channel absence, KDEG* in {3,5}. */
  void validate_grids() const;

  bool warned_rcut_ = false;  ///< one-shot guard for the truncation warning
  /** Last-resort guard (ctor + reinit), knot family only: a density knot
   *  beyond rcut is silently hard-truncated by Cut_Dummy — warn once. The
   *  honest fixes live at the command level (refit/edit/hpo derive RCUTMB
   *  from the largest knot); no auto-raise here because widening fcut
   *  without resizing the caller's neighbour list would not recover the
   *  missing pairs. */
  void warn_if_knots_beyond_rcut();

public:
  DM_REAM();
  DM_REAM(const tadah::core::Context&);
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
  std::optional<BasisKind> basis_kind() const override { return basis; }
};
} // namespace models
} // namespace tadah
#endif
