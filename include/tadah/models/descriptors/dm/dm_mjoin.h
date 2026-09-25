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
  std::vector<tadah::core::Context> contexts;
  std::string lab = "DM_mJoin";
public:
  DM_mJoin();
  ~DM_mJoin();
  DM_mJoin(tadah::core::Context &c);
  DM_mJoin(const tadah::core::Context &c);
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
  /** Every joined descriptor applies to every element pair. */
  bool applies_to_every_element_pair() const override;
  void init_rhoi(tadah::core::rho_type& rhoi) override;
  void calc_rho(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    const tadah::core::Vec3d &vec_ij,
    tadah::core::rho_type& rho,
    const double scale=1) const override;
  void set_fidx(size_t fidx_) override;
  void set_rfidx(size_t rfidx_) override;
  void set_drfidx(size_t drfidx_) override;

  void init() override;
  bool reinit(const tadah::core::Context& ctx) override;

  /** EMBF0 of the joined DM_REAM children (one global key: the children read
   *  copies of the same context); nullopt when no child has an embedding
   *  reference. */
  std::optional<EmbF0Mode> embf0_mode() const override;
  /** Sum of the children's raw F(0) (each indexes its own block of w_full);
   *  nullopt when no child reports one. */
  std::optional<double>
  raw_embedding_f0(const tadah::core::t_type &w_full) const override;
};

} // namespace models
} // namespace tadah
#endif
