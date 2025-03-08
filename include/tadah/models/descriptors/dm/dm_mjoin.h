#ifndef DM_MJOIN_H
#define DM_MJOIN_H

#include "tadah/core/registry.h"
#include <cstddef>
#include <tadah/models/descriptors/dm/dm_base.h>
#include <tadah/models/descriptors/d_mjoin.h>

namespace tadah {
namespace models {

/** \brief Meta many-body descriptor for combining multiple DM descriptors.
 *
 * This descriptor provides an interface for concatenating various many-body descriptors.
 * The resulting descriptor can then be used by Tadah! like any standard many-body descriptor.
 *
 * Each descriptor must have a specified type in a configuration file, along with a cutoff function, 
 * cutoff distance, and other optional keys that are typically expected for this descriptor,
 * such as SGRIDMB and CGRIDMB.
 *
 * When listing descriptors under the TYPEMB key, include parameters relevant to this descriptor.
 *
 * Here is an example of configuring these descriptors:
 * 
 * ```
 * TYPEMB    DM_mJoin         # Meta descriptor for concatenating many-body descriptors
 * TYPEMB    DM_EAD 1 5 5 * *    # L number, cgrid, sgrid, list of element pairs
 * RCTYPEMB  Cut_Cos
 * RCUTMB    3.0
 * CGRIDMB   -1 5 0 3.0       # Grid for DM_EAD, blips centers, auto-generated
 * SGRIDMB   -2 5 1.0 10.0    # Grid for DM_EAD, blips widths, auto-generated
 * 
 * TYPEMB    DM_Blip 0 7 7 Ta Ta   # L number, cgrid, sgrid, list of element pairs
 * RCTYPEMB  Cut_Tanh
 * RCUTMB    7.5
 * SGRIDMB   -2 7 0.1 10      # Grid for DM_Blip, blips widths, auto-generated
 * CGRIDMB   0 0 0 0 0 0 0    # Grid for DM_Blip, blips centers
 * ```
 *
 * Note: Grids can be specified on a single line, and the order of the grids should match the order of descriptors.
 *
 * There is no limit to the number of descriptors that can be concatenated.
 *
 * \details
 * - Ensure the types and grids are correctly specified in the configuration file.
 * - The cutoff functions (RCTYPEMB) and distances (RCUTMB) must be defined for each descriptor.
 * - Both SGRIDMB and CGRIDMB should be included if relevant, with their sizes matching the given descriptors.
 */
class DM_mJoin : public DM_Base, public D_mJoin {
private:
  std::vector<DM_Base*> ds;
  std::vector<Config> configs;
  std::string lab = "DM_mJoin";
public:
  DM_mJoin();
  ~DM_mJoin();
  DM_mJoin(Config &c);
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
  void set_fidx(size_t fidx_) override;
  void set_rfidx(size_t rfidx_) override;
} // namespace models
} // namespace tadah

  void init() override;
};

#endif
