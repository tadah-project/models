#ifndef D2_JOIN_H
#define D2_JOIN_H
#include <cstddef>
#include <tadah/models/descriptors/d2/d2_base.h>

/**
 * Meta two-body descriptor for joining two different D2 descriptors.
 *
 */
template<typename D1, typename D2>
class D2_Join : public D2_Base {
private:
  D1 d1;
  D2 d2;
  size_t s=0;
  //std::string lab="D2_Join";
  std::string lab;
  int verbose;
public:
  D2_Join();
  D2_Join(Config &c)
  {
    // prepare configs for d1 and d2
    // first we get joined grid and we have to split it between d1 and d2
    // we have to deal with -1 -2 algorithms and manual cases
    // TYPE2B 1st element is type followed be grid sizes (after expansion)
    v_type cgrid;
    cgrid.resize(c.size("CGRID2B"));
    c.get<v_type>("CGRID2B",cgrid);
    v_type sgrid;
    sgrid.resize(c.size("SGRID2B"));
    c.get<v_type>("SGRID2B",sgrid);
    Config c1=c;
    Config c2=c;
    c1.remove("CGRID2B");
    c1.remove("SGRID2B");
    c2.remove("CGRID2B");
    c2.remove("SGRID2B");
    try {
      size_t c1_cgrid_size = c.get<double>("CGRID2B",0) < 0 ? 4 : c.get<size_t>("TYPE2B",1);
      size_t c1_sgrid_size = c.get<double>("SGRID2B",0) < 0 ? 4 : c.get<size_t>("TYPE2B",1);
      size_t c2_cgrid_size = c.get<double>("CGRID2B",c1_cgrid_size) < 0 ? 4 : c.get<size_t>("TYPE2B",2);
      size_t c2_sgrid_size = c.get<double>("SGRID2B",c1_sgrid_size) < 0 ? 4 : c.get<size_t>("TYPE2B",2);
      for (size_t i=0; i<c1_cgrid_size; ++i) {
        c1.add("CGRID2B",cgrid[i]);
      }
      for (size_t i=0; i<c1_sgrid_size; ++i) {
        c1.add("SGRID2B",sgrid[i]);
      }
      for (size_t i=0; i<c2_cgrid_size; ++i) {
        c2.add("CGRID2B",cgrid[i+c1_cgrid_size]);
      }
      for (size_t i=0; i<c2_sgrid_size; ++i) {
        c2.add("SGRID2B",sgrid[i+c1_sgrid_size]);
      }
    }
    catch (const std::runtime_error& e) {
      std::cerr << "Grid Configuration Error: Ensure two grid generators are "
        << "used for both D2 descriptors or provide manual grids. Avoid "
        << "relying on a single auto grid generation for both descriptors. "
        << "Details: " << e.what() << std::endl;
    }
    d1 = D1(c1);
    d2 = D2(c2);

    // update main config file with generated grids
    cgrid.resize(c1.size("CGRID2B")+c2.size("CGRID2B"));
    sgrid.resize(c1.size("SGRID2B")+c2.size("SGRID2B"));
    c.remove("CGRID2B"); c.remove("SGRID2B");
    for (size_t i=0; i<c1.size("CGRID2B"); ++i) {
      c.add("CGRID2B",c1.get<double>("CGRID2B",i));
      c.add("SGRID2B",c1.get<double>("SGRID2B",i));
    }
    for (size_t i=0; i<c2.size("CGRID2B"); ++i) {
      c.add("CGRID2B",c2.get<double>("CGRID2B",i));
      c.add("SGRID2B",c2.get<double>("SGRID2B",i));
    }

    s=d1.size()+d2.size();
    lab=d1.label()+d1.label();
  }
  void calc_aed(
    const double rij,
    const double rij_sq,
    const double fc_ij,
    aed_type2 &aed) {
    d1.calc_aed(rij,rij_sq,fc_ij,aed);
    d2.calc_aed(rij,rij_sq,fc_ij,aed);
  }
  void calc_dXijdri(
    const double rij,
    const double rij_sq,
    const double fc_ij,
    const double fcp_ij,
    fd_type &fd_ij)
  {
    d1.calc_dXijdri(rij,rij_sq,fc_ij,fcp_ij,fd_ij);
    d2.calc_dXijdri(rij,rij_sq,fc_ij,fcp_ij,fd_ij);
  }

  void calc_all(
    const double rij,
    const double rij_sq,
    const double fc_ij,
    const double fcp_ij,
    aed_type2 &aed,
    fd_type &fd_ij)
  {
    d1.calc_all(rij,rij_sq,fc_ij,fcp_ij,aed,fd_ij);
    d2.calc_all(rij,rij_sq,fc_ij,fcp_ij,aed,fd_ij);
  }
  size_t size() {
    return s;
  }
  std::string label() {
    return lab;
  }
  void set_fidx(size_t fidx_) {
    d1.set_fidx(fidx_);
    d2.set_fidx(d1.size()+d1.get_fidx());
  }
};
#endif
