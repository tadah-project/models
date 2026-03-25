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
  size_t rhoisize=0;
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
    return rhoisize;
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
