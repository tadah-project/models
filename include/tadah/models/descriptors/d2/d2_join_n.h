#ifndef D2_JOIN_N_H
#define D2_JOIN_N_H

#include <tadah/models/descriptors/d2/d2_base.h>
#include <tuple>
#include <string>

/**
 * Meta multi-body descriptor for joining multiple D2 descriptors.
 *
 * <DEVELOPMENT CODE>
 * Only works with D2 without grids so essentially it is useless
 */
template<typename... Descriptors>
class D2_Join_N : public D2_Base {
private:
  std::tuple<Descriptors...> descriptors;
  size_t s = 0;
  std::string lab;

  template<size_t Index = 0>
  typename std::enable_if<Index == sizeof...(Descriptors), void>::type
  calculate_sizes() {}

  template<size_t Index = 0>
  typename std::enable_if<Index < sizeof...(Descriptors), void>::type
  calculate_sizes() {
    s += std::get<Index>(descriptors).size();
    lab += std::get<Index>(descriptors).label();
    calculate_sizes<Index + 1>();
  }

  template<size_t Index = 0>
  typename std::enable_if<Index == sizeof...(Descriptors), void>::type
  calc_aed_impl(const double, const double, const double, aed_type2 &) {}

  template<size_t Index = 0>
  typename std::enable_if<Index < sizeof...(Descriptors), void>::type
  calc_aed_impl(const double rij, const double rij_sq, const double fc_ij, aed_type2 &aed) {
    std::get<Index>(descriptors).calc_aed(rij, rij_sq, fc_ij, aed);
    calc_aed_impl<Index + 1>(rij, rij_sq, fc_ij, aed);
  }

  template<size_t Index = 0>
  typename std::enable_if<Index == sizeof...(Descriptors), void>::type
  calc_dXijdri_impl(const double, const double, const double, const double, fd_type &) {}

  template<size_t Index = 0>
  typename std::enable_if<Index < sizeof...(Descriptors), void>::type
  calc_dXijdri_impl(const double rij, const double rij_sq, const double fc_ij, const double fcp_ij, fd_type &fd_ij) {
    std::get<Index>(descriptors).calc_dXijdri(rij, rij_sq, fc_ij, fcp_ij, fd_ij);
    calc_dXijdri_impl<Index + 1>(rij, rij_sq, fc_ij, fcp_ij, fd_ij);
  }

  template<size_t Index = 0>
  typename std::enable_if<Index == sizeof...(Descriptors), void>::type
  calc_all_impl(const double, const double, const double, const double, aed_type2 &, fd_type &) {}

  template<size_t Index = 0>
  typename std::enable_if<Index < sizeof...(Descriptors), void>::type
  calc_all_impl(const double rij, const double rij_sq, const double fc_ij, const double fcp_ij, aed_type2 &aed, fd_type &fd_ij) {
    std::get<Index>(descriptors).calc_all(rij, rij_sq, fc_ij, fcp_ij, aed, fd_ij);
    calc_all_impl<Index + 1>(rij, rij_sq, fc_ij, fcp_ij, aed, fd_ij);
  }

  template<size_t Index = 0>
  typename std::enable_if<Index == sizeof...(Descriptors), void>::type
  set_fidx_impl(size_t) {}

  template<size_t Index = 0>
  typename std::enable_if<Index < sizeof...(Descriptors), void>::type
  set_fidx_impl(size_t fidx) {
    std::get<Index>(descriptors).set_fidx(fidx);
    set_fidx_impl<Index + 1>(fidx + std::get<Index>(descriptors).size());
  }

  template<typename T>
  T create_descriptor(Config &c) {
    return T(c);
  }

public:
  D2_Join_N(Config &c) : descriptors(create_descriptor<Descriptors>(c)...) {
    calculate_sizes<0>();
  }

  void calc_aed(const double rij, const double rij_sq, const double fc_ij, aed_type2 &aed) override {
    calc_aed_impl<0>(rij, rij_sq, fc_ij, aed);
  }

  void calc_dXijdri(const double rij, const double rij_sq, const double fc_ij, const double fcp_ij, fd_type &fd_ij) override {
    calc_dXijdri_impl<0>(rij, rij_sq, fc_ij, fcp_ij, fd_ij);
  }

  void calc_all(const double rij, const double rij_sq, const double fc_ij, const double fcp_ij, aed_type2 &aed, fd_type &fd_ij) override {
    calc_all_impl<0>(rij, rij_sq, fc_ij, fcp_ij, aed, fd_ij);
  }
  
  size_t size() override {
    return s;
  }
  
  std::string label() override {
    return lab;
  }

  void set_fidx(size_t fidx_) override {
    set_fidx_impl<0>(fidx_);
  }
};

#endif

