#ifndef D2_BP_H
#define D2_BP_H
#include <tadah/models/descriptors/d2/d2_base.h>
#include <vector>

namespace tadah {
namespace models {

/** \brief Behler-Parrinello two-body descriptor.
 *
 * \f[
 * V_i^{\eta,r_s} = \sum_{j \neq i} \exp{\Big(-\eta(r_{ij}-r_s)^2\Big)}f_c(r_{ij})
 * \f]
 *
 * \ref CGRID2B parameters control position \f$ r_s \f$ of the gaussian basis function.
 *
 * \ref SGRID2B parameters control width \f$ \eta \f$ of the gaussian basis function.
 *
 * This is essentially a \f$ G^1_i \f$ descriptor from the below paper with
 * an exception that it can use any cutoff function defined in Ta-dah!:
 *
 * <div class="csl-entry">Behler, J., Parrinello, M. (2007).
 * Generalized neural-network representation of high-dimensional
 * potential-energy surfaces. <i>Physical Review Letters</i>,
 * <i>98</i>(14), 146401. https://doi.org/10.1103/PhysRevLett.98.146401</div>
 *
 * Required tadah::core::Context keys:
 * \ref INIT2B \ref CGRID2B \ref SGRID2B
 */
class D2_BP : public D2_Base {
private:
  std::string lab="D2_BP";
  tadah::core::v_type cgrid;
  tadah::core::v_type sgrid;

public:
  D2_BP();
  D2_BP(const tadah::core::Context &context);
  void calc_aed(
    const int Zi,
    const int Zj,
    const double rij,
    const double ,
    tadah::core::aed_type &aed,
    const double scale=1) const override;
  void calc_dXijdri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    tadah::core::fd_type &fd_ij,
    const double scale=1) const override;
  void calc_all(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    tadah::core::aed_type &aed,
    tadah::core::fd_type &fd_ij,
    const double scale=1) const override;

  std::string label() const override;
  void init() override;
  bool reinit(const tadah::core::Context& ctx) override;

  /** Hook for the HPO PRIOR initialisation strategy.
   *
   * D2_BP ships no descriptor-specific defaults: physically-meaningful
   * values for CGRID2B / SGRID2B / RCUT2B are tightly coupled to the
   * RCUT and the user's coordination shells, and the user's bounds
   * encode that range better than any global constant could. PRIOR
   * therefore falls back to the bounds-midpoint for D2_BP parameters,
   * which is always inside the search space.
   *
   * Overriding this method in a project-specific subclass remains the
   * intended extension point for descriptor authors with strong
   * literature priors for their parameters. */
  std::optional<double>
  suggested_value(const std::string& /*key*/, int /*arg_pos*/) const override {
    return std::nullopt;
  }

  std::optional<DescriptorMeta> meta() const override;
};
} // namespace models
} // namespace tadah
#endif
