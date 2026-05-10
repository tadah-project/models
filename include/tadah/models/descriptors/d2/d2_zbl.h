#ifndef D2_ZBL_H
#define D2_ZBL_H

#include <tadah/models/descriptors/d2/d2_base.h>

namespace tadah {
namespace models {

/**
 * \brief ZBL Descriptor
 *
 * The ZBL (Ziegler-Biersack-Littmark) potential is an empirical potential used to model short-range interactions between atoms.
 *
 * The constant term \f$ \frac{e^2}{4 \pi \varepsilon_0 } \f$ is set to 1 and will be fitted as needed.
 *
 * The simplified expression for the ZBL potential is given by:
 *
 * \f[
 * V(r) = \frac{Z_1 Z_2}{r} \phi\left(\frac{r}{a}\right)
 * \f]
 *
 * where \f$ a \f$ is the screening length, expressed as:
 * \f[
 * a = \frac{s_0 a_0}{Z_1^{p_0} + Z_2^{p_1}}
 * \f]
 *
 * Here, \f$ a_0 \f$, \f$ s_0 \f$, \f$ p_0 \f$, and \f$ p_1 \f$ are adjustable hyperparameters.
 * Setting any of these to -1 uses the default values:
 *
 * - \f$ a_0 = 0.52917721067 \, \text{Å} \f$
 * - \f$ s_0 = 0.88534 \f$
 * - \f$ p_0 = 0.23 \f$
 * - \f$ p_1 = 0.23 \f$
 *
 * The screening function \f$ \phi \f$ is defined as:
 * \f[
 * \phi(x) = 0.1818 e^{-3.2x} + 0.5099 e^{-0.9423x} + 0.2802 e^{-0.4029x} + 0.02817 e^{-0.2016x}
 * \f]
 *
 * Required tadah::core::Context Key: \ref INIT2B \ref TYPE2B
 *
 * - TYPE2B D2_ZBL \f$ a_0 \f$ \f$ s_0 \f$ \f$ p_0 \f$ \f$ p_1 \f$ ELEMENT1 ELEMENT2
 *
 * Examples:
 *
 * - TYPE2B D2_ZBL 0.53 0.90 0.23 0.23 Kr Kr  # Custom parameters
 * - TYPE2B D2_ZBL -1 -1 -1 -1 Kr Kr  # Default values
 * - TYPE2B D2_ZBL 0.53 -1 -1 -1 Kr Kr  # Mix of default and custom
 *
 */
class D2_ZBL : public D2_Base {
private:
  std::string lab="D2_ZBL";
  double phi(const double x) const; // x = r/a
  // dphi/dx
  double dphi(const double x) const; // x = r/a
  double screening_length(int, int) const; // a = 
  double a0=0.52917721067; // bohr radius
  double s0=0.88534; // bohr radius
  double p0=0.23; // bohr radius
  double p1=0.23; // bohr radius
  double **a = nullptr; // precomputed screening length
  int types[119]; // local index of Zi
  int ntypes=0; // number of elements for this calc
public:
  D2_ZBL();
  D2_ZBL(const tadah::core::Context &context);
  ~D2_ZBL();
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
  std::optional<DescriptorMeta> meta() const override;
};
} // namespace models
} // namespace tadah
#endif
