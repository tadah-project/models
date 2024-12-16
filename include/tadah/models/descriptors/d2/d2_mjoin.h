#ifndef D2_MJOIN_H
#define D2_MJOIN_H
#include "tadah/core/registry.h"
#include <cstddef>
#include <tadah/models/descriptors/d2/d2_base.h>

/** \brief Meta two-body descriptor for combining multiple D2 descriptors.
 *
 * This descriptor provides a convenient interface for concatenating multiple two-body descriptors. 
 * The resulting descriptor can then be used by Tadah! like any standard two-body descriptor.
 *
 * Each descriptor must have a specified type in a configuration file, along with a cutoff function, 
 * cutoff distance, and optionally SGRID2B and CGRID2B values if applicable.
 *
 * When listing descriptors under the TYPE2B key, you must include two numbers for each descriptor 
 * that indicate the sizes of the CGRID2B and SGRID2B respectively. In some cases, these numbers can be zero. 
 * The numbers must correspond to the lengths of the provided grids/parameters in the configuration file.
 *
 * Here is an example of how to configure these descriptors:
 * 
 * ```
 * TYPE2B    D2_mJoin D2_MIE 0 2 D2_Blip 6 6 
 * RCTYPE2B  Cut_Cos Cut_Tanh
 * RCUT2B    3.0 7.5
 * SGRID2B   11 6   # D2_MIE parameters
 * SGRID2B   -2 6 0.1 10   # Grid for D2_Blip, blips widths
 * CGRID2B   0 0 0 0 0 0   # Grid for D2_Blip, blip centers
 * ```
 *
 * Note: Grids can be specified on a single line, and the order of the grids is important. 
 * In the example above, the first grid corresponds to D2_MIE, followed by the grid for D2_Blip.
 *
 * There is no limit to the number of descriptors that can be concatenated.
 *
 * \details
 * - Ensure the types and grids are correctly specified in the configuration file.
 * - The cutoff functions (RCTYPE2B) and distances (RCUT2B) must be defined for each descriptor.
 * - Both SGRID2B and CGRID2B should be included if relevant, with their sizes matching the given descriptors.
 */
class D2_mJoin : public D2_Base {
private:
  std::vector<D2_Base*> ds;
  std::vector<Config> configs;
  std::string lab = "D2_mJoin";
public:
  D2_mJoin();
  ~D2_mJoin();
  D2_mJoin(Config &c);
  
  void calc_aed(const int Zj, const double rij, const double rij_sq, aed_type2 &aed, const double scale=1) override;
  void calc_dXijdri(const int Zj, const double rij, const double rij_sq, fd_type &fd_ij, const double scale=1) override;
  void calc_all(const int Zj, const double rij, const double rij_sq, aed_type2 &aed, fd_type &fd_ij, const double scale=1) override;
  std::string label() override;
  void set_fidx(size_t fidx_) override;
};
#endif
