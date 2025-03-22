#ifndef DM_EAM_H
#define DM_EAM_H
#include <tadah/models/descriptors/dm/dm_base.h>

namespace tadah {
namespace models {

/** \brief many-body part for the Embedded Atom Method descriptor.
 *
 * \f[
 * V_i = F\Bigg(\sum_{j \neq i} \rho(r_{ij}) \Bigg)
 * \f]
 *
 * This descriptor will load tabulated values for the density \f$ \rho \f$
 * and embedded energy \f$ F \f$ from the provided \ref SETFL file.
 *
 * This descriptor is usually used together with the two-body descriptor \ref D2_EAM
 * although this is not required and user can mix it with any other descriptors
 * or use it on its own.
 *
 * This descriptor will enforce cutoff distance as specified in a \ref SETFL file.
 * Set \ref RCUTMB to the same value to suppress the warning message.
 *
 * Required tadah::core::Context keys:
 * \ref INITMB \ref SETFL
 */
class DM_EAM: public DM_Base {
private:
  struct eam_file {
    std::string file_path;
    std::vector<double> frho;
    std::vector<double> rhor;
    std::vector<double> z2r;
    int nrho=0;
    double drho=0;
    int nr;
    double dr;
    double rdr;
    double rdrho;
    double rcut;
    int atomic_number;
    double atomic_mass;
    double lattice_param;
    std::string lattice;
    double rhomax;

  };
  eam_file ef;

  std::vector<std::vector<double>> frho_spline;
  std::vector<std::vector<double>> rhor_spline;
  std::vector<std::vector<double>> z2r_spline;

  std::string lab="DM_EAM";
  void read_setfl();
  void gen_splines(int &n, double &delta, std::vector<double> &f, std::vector<std::vector<double>> &spline);

public:
  DM_EAM();
  DM_EAM(tadah::core::Context&);
  void calc_aed(
    tadah::core::rho_type& rho,
    tadah::core::aed_type &aed) override;
  void calc_dXijdri_dXjidri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d &vec_ij,
    tadah::core::rho_type& rhoi,
    tadah::core::rho_type& rhoj,
    tadah::core::fd_type &fd_ij,
    const double scale=1) override;
  void calc_dXijdri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d &vec_ij,
    tadah::core::rho_type& rhoi,
    tadah::core::fd_type &fd_ij, 
    const double scale=1) override;
  std::string label() override;
  void init_rhoi(tadah::core::rho_type& rhoi) override;
  void calc_rho(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d &vec_ij,
    tadah::core::rho_type& rho,

    const double scale=1) override;
  void init() override;
};
} // namespace models
} // namespace tadah
#endif
