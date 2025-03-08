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
 * Required Config keys:
 * \ref INIT2B \ref CGRID2B \ref SGRID2B
 */
class D2_BP : public D2_Base {
private:
  std::string lab="D2_BP";
  v_type etas;
  v_type mius;

public:
  D2_BP();
  D2_BP(Config &config);
  void calc_aed(
    const int Zi,
    const int Zj,
    const double rij,
    const double ,
    aed_type &aed,
    const double scale=1) override;
  void calc_dXijdri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    fd_type &fd_ij,
    const double scale=1) override;
  void calc_all(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    aed_type &aed,
    fd_type &fd_ij,
    const double scale=1) override;
} // namespace models
} // namespace tadah

  std::string label() override;
  void init() override;
};
#endif
