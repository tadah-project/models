#ifndef D2_MIE_H
#define D2_MIE_H

#include <tadah/models/descriptors/d2/d2_base.h>

/** \brief Mie descriptor
 *
 * \f[
 * V_i = \sum_{j \neq i} C \epsilon \Bigg(\Big(\frac{\sigma}{r_{ij}}\Big)^{n} - \Big(\frac{\sigma}{r_{ij}}\Big)^m\Bigg)
 * \f]
 *
 * where
 * \f[
 * C=\frac{n}{n-m}\Big( \frac{n}{m} \Big)^{\frac{m}{n-m}}
 * \f]
 *
 *
 * Any cutoff can be used
 *
 * Required Config Key: \ref INIT2B \ref SGRID2B
 *
 * e.g. SGID2B 12 6
 *
 * will result in Lennard-Jones type descriptor
 */
class D2_MIE : public D2_Base {
private:
  std::string lab="D2_MIE";
  double m,n;
  v_type etas;
public:
  D2_MIE();
  D2_MIE(Config &config);
  void calc_aed(
    const int Zj,
    const double rij,
    const double ,
    aed_type &aed,
    const double scale=1) override;
  void calc_dXijdri(
    const int Zj,
    const double rij,
    const double rij_sq,
    fd_type &fd_ij,
    const double scale=1) override;
  void calc_all(
    const int Zj,
    const double rij,
    const double rij_sq,
    aed_type &aed,
    fd_type &fd_ij,
    const double scale=1) override;
  std::string label() override;
  void init() override;
};
#endif
