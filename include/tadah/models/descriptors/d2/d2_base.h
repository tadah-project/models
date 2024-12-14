#ifndef D2_BASE_H
#define D2_BASE_H

#include <tadah/models/descriptors/d_base.h>

/** \brief Base class for all two-body type descriptors.
 *
 * All two-body descriptors **must** inherit this class.
 */
class D2_Base: public D_Base {
protected:
  size_t fidx=0;    // first index in aed and fd arrays for this descriptor
public:
  virtual ~D2_Base() {};

  /** \brief Calculate \ref AED
   *
   * Calculate Atomic Energy Descriptor for the atom local environment.
   */
  virtual void calc_aed(
    const double rij,
    const double rij_sq,
    const double fc_ij,
    aed_type2& aed)=0;

  /** \brief Calculate \ref FD
   *
   * Calculate Force Descriptor for the atom local environment.
   *
   * Computes x-direction only.
   * The resulting vector must be scaled by the unit directional vector delij/rij
   */
  virtual void calc_dXijdri(
    const double rij,
    const double rij_sq,
    const double fc_ij,
    const double fcp_ij,
    fd_type &fd_ij)=0;

  /** \brief Calculate \ref AED + \ref FD */
  virtual void calc_all(
    const double rij,
    const double rij_sq,
    const double fc_ij,
    const double fcp_ij,
    aed_type2& aed,
    fd_type &fd_ij)=0;
  virtual std::string label()=0;

  /** Central difference approximation to \ref FD 
   *
   * Calculate Force Descriptor for the atom local environment
   * using central difference approximation.
   * 
   * fc: a pointer to initialised cutoff object.
   * r: a separation between two particles
   * fd: appropriate size container to store computation
   * h: central difference parameter
   */
  template <typename CUT>
  void calc_fd_approx(CUT *fc, double r, fd_type &fd, double h=1e-8) {
    aed_type2 aed_n((*this).size()); // r+h
    aed_type2 aed_p((*this).size()); // r-h
    aed_n.set_zero();
    aed_p.set_zero();
    double rn = r-h;
    double rp = r+h;
    double fcp = fc->calc(rp);
    double fcn = fc->calc(rn);
    (*this).calc_aed(rn,rn*rn,fcn,aed_n);
    (*this).calc_aed(rp,rp*rp,fcp,aed_p);
    fd(0) = (0.5/h)*(aed_p - aed_n);
  }
  virtual void set_fidx(size_t fidx_) { fidx=fidx_; }
  virtual size_t get_fidx() { return fidx; }
};
#endif
