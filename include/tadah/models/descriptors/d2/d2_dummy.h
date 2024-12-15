#ifndef D2_DUMMY_H
#define D2_DUMMY_H
#include <tadah/models/descriptors/d2/d2_base.h>

/**
 * Dummy two-body descriptor.
 *
 * Use it to satisfy DescriptorsCalc requirements in case
 * when two-body descriptor is not required.
 *
 */
class D2_Dummy : public D2_Base {
private:
  std::string lab="D2_Dummy";
public:
  D2_Dummy();
  D2_Dummy(Config &);
  void calc_aed(
    const int Zj,
    const double rij,
    const double ,
    aed_type2 &aed,
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
    aed_type2 &aed,
    fd_type &fd_ij,
    const double scale=1) override;
  std::string label() override;
};
#endif
