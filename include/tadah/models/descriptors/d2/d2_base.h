#ifndef D2_BASE_H
#define D2_BASE_H

#include <tadah/models/cut_all.h>
#include <tadah/models/descriptors/d_base.h>

namespace tadah {
namespace models {

/** \brief Base class for all two-body type descriptors.
 *
 * All two-body descriptors **must** inherit this class.
 */
class D2_Base: public D_Base {
public:
  D2_Base() {};
  D2_Base(const tadah::core::Context &c): D_Base(c) {
    //if (c.get<bool>("INIT2B")) {
    if (c.contains("RCUT2BMAX") && c.contains("RCTYPE2B")) {
      double rcut2b = c.get<double>("RCUT2BMAX");
      set_fcut(tadah::core::factory<Cut_Base,double>( c.get<std::string>("RCTYPE2B"), rcut2b ), true);
    }
    else {
      set_fcut(new Cut_Dummy(0), true);
    }
  }
  virtual ~D2_Base() {};
  virtual std::vector<std::string> get_init_atoms(const tadah::core::Context &c) override;

  /** \brief Calculate \ref AED
   *
   * Calculate Atomic Energy Descriptor for the atom local environment.
   */
  virtual void calc_aed(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    tadah::core::aed_type& aed,
    const double scale=1)=0;

  /** \brief Calculate \ref FD
   *
   * Calculate Force Descriptor for the atom local environment.
   *
   * For two body this is essentiall df/dr where f is a descriptor function.
   * Do not include negative sign (as you mifht think E=-df/dr). The negative sign 
   * is included elsewhere in the code.
   *
   * Computes x-direction only.
   *
   * The resulting vector must be scaled by the unit directional vector delij/rij
   * where delij = r_i - r_j, and rij is |r_i - r_j|
   */
  virtual void calc_dXijdri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    tadah::core::fd_type &fd_ij,
    const double scale=1)=0;

  /** \brief Calculate \ref AED + \ref FD */
  virtual void calc_all(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    tadah::core::aed_type& aed,
    tadah::core::fd_type &fd_ij,
    const double scale=1)=0;

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
  void calc_fd_approx(const int Zi, const int Zj, double r, tadah::core::fd_type &fd, double h=1e-8) {
    tadah::core::aed_type aed_n((*this).size()); // r+h
    tadah::core::aed_type aed_p((*this).size()); // r-h
    aed_n.set_zero();
    aed_p.set_zero();
    double rn = r-h;
    double rp = r+h;
    (*this).calc_aed(Zi,Zj,rn,rn*rn,aed_n);
    (*this).calc_aed(Zi,Zj,rp,rp*rp,aed_p);

    fd(0) = (0.5/h)*(aed_p - aed_n);
  }
};
} // namespace models
} // namespace tadah
#endif
