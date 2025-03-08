#ifndef D2_JOIN_H
#define D2_JOIN_H
#include <cstddef>
#include <tadah/models/descriptors/d2/d2_base.h>

namespace tadah {
namespace models {

/**
 * Meta two-body descriptor for joining two different D2 descriptors.
 *
 */
template<typename D1, typename D2>
class D2_Join : public D2_Base {
private:
  D1 *d1 = nullptr;
  D2 *d2 = nullptr;
  //std::string lab="D2_Join";
  std::string lab;
public:
  D2_Join();
  ~D2_Join() {
    if (d1) delete d1;
    if (d2) delete d2;
  }
  D2_Join(Config &c): D2_Base(c)
  {
    // prepare configs for d1 and d2
    // first we get joined grid and we have to split it between d1 and d2
    // we have to deal with -1 -2 algorithms and manual cases
    // TYPE2B 1st element is type followed be grid sizes (after expansion)
    //
    // TODO
    // We want to be able to use this class in a recurence realtion
    // where d1 <- any D2 but not D2_Join
    // and d2 can be any including D2_Join
    // This means that we trim config for d1 and pass
    // whatever is left to d2
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
      size_t c2_cgrid_size = c.size("CGRID2B") - c1_cgrid_size;
      size_t c2_sgrid_size = c.size("SGRID2B") - c1_sgrid_size;
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
    // deal with cutoffs
    c1.remove("RCTYPE2B"); c2.remove("RCTYPE2B");
    c1.remove("RCUT2B"); c2.remove("RCUT2B");
    std::vector<std::string> rctype2b(c.size("RCTYPE2B"));
    std::vector<std::string> rcut2b(c.size("RCUT2B"));
    c.get("RCTYPE2B",rctype2b);
    c.get("RCUT2B",rcut2b);
    c1.add("RCTYPE2B", rctype2b[0]);
    c1.add("RCUT2B", rcut2b[0]);
    rctype2b.erase(rctype2b.begin());
    rcut2b.erase(rcut2b.begin());
    c2.add("RCTYPE2B", rctype2b);
    c2.add("RCUT2B", rcut2b);
    // std::cout << c1 << std::endl;
    // std::cout << c2 << std::endl;
    d1 = new D1(c1);
    d2 = new D2(c2);

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

    s=d1->size()+d2->size();
    lab="KUTAbezS"+d1->label()+d1->label();
  }
  void calc_aed(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    aed_type &aed,
    const double scale=1) override {
    d1->calc_aed(Zi,Zj,rij,rij_sq,aed,scale);
    d2->calc_aed(Zi,Zj,rij,rij_sq,aed,scale);
  }
  void calc_dXijdri(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    fd_type &fd_ij,
    const double scale=1) override {
    d1->calc_dXijdri(Zi,Zj,rij,rij_sq,fd_ij,scale);
    d2->calc_dXijdri(Zi,Zj,rij,rij_sq,fd_ij,scale);
  }
  void calc_all(
    const int Zi,
    const int Zj,
    const double rij,
    const double rij_sq,
    aed_type &aed,
    fd_type &fd_ij,
    const double scale=1) override {
    d1->calc_all(Zi,Zj,rij,rij_sq,aed,fd_ij,scale);
    d2->calc_all(Zi,Zj,rij,rij_sq,aed,fd_ij,scale);
  }
  std::string label() override {
    return lab;
  }
  void set_fidx(size_t fidx_) override {
    d1->set_fidx(fidx_);
    d2->set_fidx(d1->size()+d1->get_fidx());
} // namespace models
} // namespace tadah

  }
  void init() override {};
};
#endif
