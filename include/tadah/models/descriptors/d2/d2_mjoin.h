#ifndef D2_MJOIN_H
#define D2_MJOIN_H
#include "tadah/core/registry.h"
#include <cstddef>
#include <tadah/models/descriptors/d2/d2_base.h>

/** @brief Meta two-body descriptor for combining multiple D2 descriptors.
 *
 * This descriptor provides a convenient interface for concatenating multiple two-body descriptors. 
 * The resulting descriptor can then be used by Tadah! like any standard two-body descriptor.
 *
 * Each descriptor must have a specified type in a configuration file, along with a cutoff function, 
 * cutoff distance, and optionally SGRID2B and CGRID2B values if applicable.
 *
 * When listing descriptors under the TYPE2B key, you must include two numbers for each descriptor 
 * that indicate the sizes of the CGRID2B and SGRID2B respectively. In some cases, these numbers can be zero. 
 * The numbers must correspond to the lengths of the provided grids in the configuration file.
 *
 * @example 
 * Here is an example of how to configure these descriptors:
 * 
 * ```
 * TYPE2B    D2_mJoin D2_MIE 0 2 D2_Blip 6 6 
 * RCTYPE2B  Cut_Cos Cut_Tanh
 * RCUT2B    3.0 7.5
 * SGRID2B   11 6   # Grid for D2_MIE
 * SGRID2B   -2 6 0.1 10   # Grid for D2_Blip
 * CGRID2B   0 0 0 0 0 0   # Grid for D2_Blip
 * ```
 *
 * Note: Grids can be specified on a single line, and the order of the grids is important. 
 * In the example above, the first grid corresponds to D2_MIE, followed by the grid for D2_Blip.
 *
 * There is no limit to the number of descriptors that can be concatenated.
 *
 * @details
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
  ~D2_mJoin() {
    for (auto d:ds) if (d) delete d;
  }
  D2_mJoin(Config &c): D2_Base(c)
  {
    v_type cgrid;
    cgrid.resize(c.size("CGRID2B"));
    c.get<v_type>("CGRID2B",cgrid);
    v_type sgrid;
    sgrid.resize(c.size("SGRID2B"));
    c.get<v_type>("SGRID2B",sgrid);
    size_t N = (c.size("TYPE2B")-1)/3;

    size_t cgrid_size=0;
    size_t sgrid_size=0;
    std::vector<std::string> type2b(c.size("TYPE2B"));
    c.get("TYPE2B",type2b);
    std::vector<std::string> rctype2b(c.size("RCTYPE2B"));
    std::vector<std::string> rcut2b(c.size("RCUT2B"));
    c.get("RCTYPE2B",rctype2b);
    c.get("RCUT2B",rcut2b);

    for (size_t i=0; i<N; ++i) {
      Config c1=c;
      c1.remove("CGRID2B"); c1.remove("SGRID2B"); c1.remove("TYPE2B");
      c1.remove("RCTYPE2B");  c1.remove("RCUT2B");
      try {
        // we have to deal with len(CGRID2B)!=len(SGRID2B)
        if (cgrid_size<c.size("CGRID2B")) {
          size_t c1_cgrid_size = c.get<double>("CGRID2B",cgrid_size) < 0 ? 4 : c.get<size_t>("TYPE2B",3*i+2);
          for (size_t j=0; j<c1_cgrid_size; ++j) {
            c1.add("CGRID2B",cgrid[j+cgrid_size]);
          }
          cgrid_size+=c1_cgrid_size;
        }
        if (sgrid_size<c.size("SGRID2B")) {
          size_t c1_sgrid_size = c.get<double>("SGRID2B",sgrid_size) < 0 ? 4 : c.get<size_t>("TYPE2B",3*i+3);
          for (size_t j=0; j<c1_sgrid_size; ++j) {
            c1.add("SGRID2B",sgrid[j+sgrid_size]);
          }
          sgrid_size+=c1_sgrid_size;
        }
      }
      catch (const std::runtime_error& e) {
        std::cerr << "Grid Configuration Error: Ensure two grid generators are "
          << "used for both D2 descriptors or provide manual grids. Avoid "
          << "relying on a single auto grid generation for multiple descriptors. "
          << "Details: " << e.what() << std::endl;
      }

      // deal with types
      c1.add("TYPE2B", type2b[3*i+1]);
      // deal with cutoffs
      c1.add("RCTYPE2B", rctype2b[i]);
      c1.add("RCUT2B", rcut2b[i]);
      ds.push_back(CONFIG::factory<D2_Base,Config&>( c1.get<std::string>("TYPE2B"), c1 ));
      configs.push_back(c1);
      s+=ds.back()->size();
    }

    // update main config file with generated grids
    cgrid.resize(cgrid_size);
    sgrid.resize(sgrid_size);
    c.remove("CGRID2B"); c.remove("SGRID2B");
    for (const auto &c1 : configs) {
      if (c1.exist("CGRID2B"))
        for (size_t i=0; i<c1.size("CGRID2B"); ++i)
          c.add("CGRID2B",c1.get<double>("CGRID2B",i));
      if (c1.exist("SGRID2B"))
        for (size_t i=0; i<c1.size("SGRID2B"); ++i)
          c.add("SGRID2B",c1.get<double>("SGRID2B",i));
    }
  }
  void calc_aed(
    const int Zj,
    const double rij,
    const double rij_sq,
    aed_type2 &aed,
    const double scale=1) override {
    for (auto d : ds) d->calc_aed(Zj,rij,rij_sq,aed,scale);
  }
  void calc_dXijdri(
    const int Zj,
    const double rij,
    const double rij_sq,
    fd_type &fd_ij,
    const double scale=1) override {
    for (auto d : ds) d->calc_dXijdri(Zj,rij,rij_sq,fd_ij,scale);
  }
  void calc_all(
    const int Zj,
    const double rij,
    const double rij_sq,
    aed_type2 &aed,
    fd_type &fd_ij,
    const double scale=1) override {
    for (auto d : ds) d->calc_all(Zj,rij,rij_sq,aed,fd_ij,scale);
  }
  std::string label() override {
    return lab;
  }
  void set_fidx(size_t fidx_) override {
    ds[0]->set_fidx(fidx_);
    size_t s=0;
    for (size_t i=1; i<ds.size(); ++i) {
      s+=ds[i-1]->size();
      ds[i]->set_fidx(fidx_+s);
    }
  }
};
#endif

