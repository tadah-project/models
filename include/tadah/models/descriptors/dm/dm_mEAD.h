#ifndef DM_mEAD_H
#define DM_mEAD_H
#include <tadah/models/descriptors/dm/dm_base.h>


/** Modified Embedded Atom Descriptor
 *
 * REQUIRED KEYS: SGRIDMB, CGRIDMB + KEYS OF THE THE EMBEDDING FUNCTION
 
 * This descriptor has mathematical form very similar to DM_EAD.
 * but allow usage of custom defined embedding function, \f$ \mathcal{F} \f$, in place of default quadratic one.
 *
 * \f[
 * V_i^{L,\eta,r_s} = \sum_{l_x,l_y,l_z}^{l_x+l_y+l_z=L} \frac{L!}{l_x!l_y!l_z!}
 * \mathcal{F}\Big( \rho_i^{\eta,r_s,l_x,l_y,l_z} \Big)
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
 *
 * e.g. \f$L_{max}=2\f$ will calculate descriptors with \f$ L=0,1,2 \f$ (s,p,d orbitals).
 *
 *
 * Required Config keys:
 * \ref INITMB \ref CGRIDMB \ref SGRIDMB
 */
template <typename F>
class DM_mEAD: public DM_Base {
private:
  int Lmax=0;
  std::string lab="DM_mEAD";
  v_type sgrid;
  v_type cgrid;
  std::vector<std::vector<std::vector<int>>> orbitals;
  std::vector<std::vector<double>> fac;
  size_t gen_atomic_orbitals(int L);
  double fact(int n);
  double my_pow(double,int);
  double rc;
  F emb;

public:
  DM_mEAD();
  DM_mEAD(Config&);
  void calc_aed(
    rho_type& rho,
    aed_type &aed) override;
  void calc_dXijdri_dXjidri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const Vec3d &vec_ij,
    rho_type& rhoi,
    rho_type& rhoj,
    fd_type &fd_ij,
    const double scale=1) override;
  void calc_dXijdri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const Vec3d &vec_ij,
    rho_type& rhoi,
    fd_type &fd_ij, 
    const double scale=1) override;
  std::string label() override;
  void init_rhoi(rho_type& rhoi) override;
  void calc_rho(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const Vec3d &vec_ij,
    rho_type& rho,
    const double scale=1) override;
  void init() override;
};

#include "dm_mEAD.hpp"
#endif
