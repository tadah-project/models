#ifndef DM_BASE_H
#define DM_BASE_H

#include <tadah/models/descriptors/d_base.h>
#include <string>

namespace tadah {
namespace models {

/** \brief Base class for all many-body type descriptors.
 *
 * All many-body descriptors **must** inherit this class.
 */
class DM_Base: public D_Base {

protected:
  size_t rfidx=0; // index of the first density component in the rhoi array
  size_t drfidx=0; // index of the first component of the derivative of the density component in the rhoi array
  size_t rhoisize=0;   // size of the density (rho) part
  size_t rhoipsize=0;  // size of the derivative (dF/drho) part
public:
  DM_Base() {};
  DM_Base(const tadah::core::Context &c): D_Base(c) {
    //if (c.get<bool>("INITMB")) {
    if (c.contains("RCUTMBMAX") && c.contains("RCTYPEMB")) {
      double rcutmb = c.get<double>("RCUTMBMAX");
      set_fcut(tadah::core::factory<Cut_Base,double&>( c.get<std::string>("RCTYPEMB"), rcutmb ), true);
    }
    else {
      set_fcut(new Cut_Dummy(0), true);
    }
  }
  virtual ~DM_Base() {};
  virtual std::vector<std::string> get_init_atoms(const tadah::core::Context &c) override;

  /** \brief Calculate \ref AED
   *
   * Calculate Atomic Energy Descriptor for the atom local environment.
   */
  virtual void calc_aed(
    tadah::core::rho_type& rho,
    tadah::core::aed_type &aed) const=0;

  /** \brief Calculate \ref FD
   *
   * Calculate Force Descriptor between to atoms.
   *
   * This method works for half NN lists
   * and linear models.
   *
   * \f[
   * \mathbf{fd}_{ij} \mathrel{+}= 
   * \frac{\partial \mathbf{X}_{ij}}{\partial \mathbf{r}_i} +
   * \frac{\partial \mathbf{X}_{ji}}{\partial \mathbf{r}_i}
   * \f]
   */
  virtual void calc_dXijdri_dXjidri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d &vec_ij,
    tadah::core::rho_type& rhoi,
    tadah::core::rho_type& rhoj,
    tadah::core::fd_type &fd_ij,
    const double scale=1) const=0;

  /** \brief Calculate \ref FD
   *
   * Calculate Force Descriptor between to atoms.
   *
   * This method works for full NN lists
   * and all models.
   *
   * mode=0 inidcates only x-direction is computed and 
   * the resulting descriptor must be multiplied by delx/r.
   *
   * \f[
   * \mathbf{fd}_{ij} \mathrel{+}=
   * \frac{\partial \mathbf{X}_{ij}}{\partial \mathbf{r}_i}
   * \f]
   */
  virtual void calc_dXijdri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d &vec_ij,
    tadah::core::rho_type& rhoi,
    tadah::core::fd_type &fd_ij,
    const double scale=1) const=0;

  /** \brief Resize arrays for a density and F' */
  virtual void init_rhoi(tadah::core::rho_type& rhoi)=0;

  /** \brief Calculate density */
  virtual void calc_rho(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d &vec_ij,
    tadah::core::rho_type& rho,
    const double scale=1) const =0;

  /** \brief Closed-form per-pair contribution to drho_i/dtheta.
   *
   *  Many-body analog of D2_Base::calc_dDdtheta. The MB design-matrix
   *  builder accumulates `rho_i` in a per-pair pre-pass (calc_rho),
   *  then computes `aed = F(rho)` in a per-atom finalisation step
   *  (calc_aed). The dPhi/dtheta walker mirrors that structure with
   *  TWO methods:
   *
   *   1. calc_drhoi_dtheta  — per-pair drho_i/dtheta accumulator
   *      (analog of calc_rho).
   *   2. calc_daed_dtheta   — per-atom (rho_i, drho_i/dtheta) -> daed/dtheta
   *      reducer (analog of calc_aed).
   *
   *  Calling convention mirrors D2_Base::calc_dDdtheta:
   *   - `theta_name` is the OPTIM key, upper-case (e.g. "CGRIDMB").
   *   - `arg_pos` is 1-based for vector parameters, -1 for scalars.
   *   - `drhoi` is the per-atom accumulator, pre-zeroed by the caller.
   *   - Returns `true` iff this call wrote anything into `drhoi`.
   *   - Cross-cutting parameters (RCUTMB) flow through
   *     `fcut->calc_drcut`; each DM descriptor that uses a non-Cut_Dummy
   *     cutoff contributes additively.
   *
   *  Default returns false (no-op). DM_Dummy and DM descriptors that
   *  are numerical-only (DM_EAM, DM_mEAD) keep the default.
   */
  virtual bool calc_drhoi_dtheta(
      const std::string& /*theta_name*/,
      int                /*arg_pos*/,
      const int          /*Zi*/,
      const int          /*Zj*/,
      const double       /*rij*/,
      const double       /*rij_sq*/,
      const tadah::core::Vec3d& /*vec_ij*/,
      tadah::core::rho_type& /*drhoi*/,
      const double       /*scale=1*/) const {
    return false;
  }

  /** \brief Closed-form per-atom contribution to daed/dtheta.
   *
   *  Companion to calc_drhoi_dtheta. Called by the dPhi/dtheta walker
   *  once per atom AFTER:
   *   - calc_rho        has filled `rhoi[rfidx..]   = rho_i`,
   *   - calc_aed        has filled `rhoi[drfidx..]  = F'(rho_i)`,
   *   - calc_drhoi_dtheta has filled `drhoi[rfidx..] = drho_i/dtheta`.
   *
   *  Chain rule for F(rho_i, theta):
   *
   *    daed/dtheta = (dF/drho) (drho/dtheta) + (dF/dtheta)|_{rho fixed}
   *
   *  The first term reads F'(rho) from rhoi[drfidx..] and dotproducts
   *  with drhoi[rfidx..]. The second term is non-zero only for
   *  cross-cutting parameters that enter the per-L multinomial
   *  prefactor `f[L,o] = L!/(lx! ly! lz!) / (4 rc)^L` — namely RCUTMB.
   *
   *  Returns `true` iff this call wrote anything into `daed_out`.
   *
   *  Default returns false.
   */
  virtual bool calc_daed_dtheta(
      const std::string& /*theta_name*/,
      int                /*arg_pos*/,
      const tadah::core::rho_type& /*rhoi*/,
      const tadah::core::rho_type& /*drhoi*/,
      tadah::core::aed_type& /*daed_out*/,
      const double       /*scale=1*/) const {
    return false;
  }

  /** \brief Closed-form per-pair contribution to d(fd_ij)/dtheta — the
   *  many-body sibling of D2_Base::calc_dXijdri_dtheta. Walker calls
   *  this once per (atom i, neighbour j) pair AFTER:
   *   - calc_rho            has filled `rhoi[rfidx..]   = rho_i`,
   *   - calc_aed            has filled `rhoi[drfidx..]  = f[L,o] F'(rho_i)`,
   *   - calc_drhoi_dtheta   has filled `drhoi_dtheta[rfidx..]
   *                                       = drho_i/dtheta`.
   *
   *  Unlike the D2 sibling (which returns a single radial column to be
   *  axis-expanded by the walker), the MB descriptor's calc_dXijdri
   *  writes all three axis components directly because the per-orbital
   *  polynomial poly = x^lx * y^ly * z^lz couples to direction. This
   *  override mirrors that signature — it writes the full axis-expanded
   *  d(fd_ij)/dtheta into `fd_ij_out(row, 0..2)`.
   *
   *  Chain rule sketch (for DM_mEAD<F>, expanding path):
   *
   *     fd_ij(p, k) = scale * w_j * f[L,o] * F'(rho_i, p) * term_k,
   *       term_k = g * (dpoly/dr_k) + dg * poly * (v_k / r)
   *
   *  with g, dg the radial Gaussian and its r-derivative (both
   *  cutoff-weighted). Differentiating w.r.t. theta sees three groups:
   *
   *    1. Radial theta (CGRIDMB[k], SGRIDMB[k], RCUTMB) — touch (g, dg)
   *       directly AND rho via calc_drhoi_dtheta. The F'(rho) chain
   *       needs F''(rho) via F_Base::calc_d2F.
   *    2. Embedding theta (CEMBFUNC[i], SEMBFUNC[i]) — touch F'
   *       directly via F_Base::calc_dFp_dtheta; rho is fixed.
   *    3. Cross-cutting RCUTENV (F_mEnv only) — same direct F' branch.
   *
   *  Calling convention mirrors calc_drhoi_dtheta / calc_daed_dtheta:
   *   - `theta_name` is the OPTIM key, upper-case.
   *   - `arg_pos` is 1-based for vector parameters, -1 for scalars.
   *   - `fd_ij_out` is the per-pair, axis-expanded buffer, pre-zeroed by
   *     the caller.
   *   - Returns `true` iff this call wrote anything into `fd_ij_out`.
   *
   *  Default returns false (descriptor doesn't own this theta on the
   *  force side; orchestrator routes to needs_fd). */
  virtual bool calc_dXijdri_dtheta(
      const std::string& /*theta_name*/,
      int                /*arg_pos*/,
      const int          /*Zi*/,
      const int          /*Zj*/,
      const double       /*rij*/,
      const double       /*rij_sq*/,
      const tadah::core::Vec3d& /*vec_ij*/,
      const tadah::core::rho_type& /*rhoi*/,
      const tadah::core::rho_type& /*drhoi_dtheta*/,
      tadah::core::fd_type& /*fd_ij_out*/,
      const double       /*scale=1*/) const {
    return false;
  }


  /** \brief This is just a convenient wrapper.
   */
  void calc_dXijdri(const int Zi, const int Zj, tadah::core::Vec3d &del, tadah::core::rho_type &rhoi, tadah::core::fd_type &fdij, double scale=1) const{
    double rsq = del*del;
    double r = sqrt(rsq);
    //(p->*func)(r,rsq,del,1.0,0.0,rhoi,fdij);
    (*this).calc_dXijdri(Zi, Zj,r,rsq,del,rhoi,fdij,scale);
  }

  /** \brief This is just a convenient wrapper.
   */
  void calc_rho_mb(const int Zi, const int Zj, tadah::core::Vec3d &del, tadah::core::rho_type &rhoi, double scale=1) const {
    double rsq = del*del;
    double r = sqrt(rsq);
    ((*this).calc_rho)(Zi, Zj,r,rsq,del,rhoi,scale);
  }
  /** \brief Return size of the density array */
  virtual size_t rhoi_size() const{
    return rhoisize;
  }
  /** \brief Return size of the derivative of the embedding energy array */
  virtual size_t rhoip_size() const{
    return rhoipsize > 0 ? rhoipsize : rhoisize;
    // return rhoipsize;
  }
  // virtual size_t size() {
  //   return s;
  // }
  virtual std::string label() const=0;

  virtual size_t get_rfidx() const;
  virtual void set_rfidx(size_t rfidx_);
  virtual size_t get_drfidx() const;
  virtual void set_drfidx(size_t drfidx_);
};
} // namespace models
} // namespace tadah
#endif
