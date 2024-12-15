#include <tadah/models/descriptors/d2/d2_mjoin.h>
#include <iostream>
#include <stdexcept>

D2_mJoin::D2_mJoin() {}

D2_mJoin::~D2_mJoin() {
  for (auto d : ds) if (d) delete d;
}

D2_mJoin::D2_mJoin(Config &c) : D2_Base(c) {
  v_type cgrid;
  cgrid.resize(c.size("CGRID2B"));
  c.get<v_type>("CGRID2B", cgrid);
  v_type sgrid;
  sgrid.resize(c.size("SGRID2B"));
  c.get<v_type>("SGRID2B", sgrid);
  size_t N = (c.size("TYPE2B") - 1) / 3;

  size_t cgrid_size = 0;
  size_t sgrid_size = 0;
  std::vector<std::string> type2b(c.size("TYPE2B"));
  c.get("TYPE2B", type2b);
  std::vector<std::string> rctype2b(c.size("RCTYPE2B"));
  std::vector<std::string> rcut2b(c.size("RCUT2B"));
  c.get("RCTYPE2B", rctype2b);
  c.get("RCUT2B", rcut2b);

  for (size_t i = 0; i < N; ++i) {
    Config c1 = c;
    c1.remove("CGRID2B"); c1.remove("SGRID2B"); c1.remove("TYPE2B");
    c1.remove("RCTYPE2B"); c1.remove("RCUT2B");

    try {
      if (cgrid_size < c.size("CGRID2B")) {
        size_t c1_cgrid_size = c.get<double>("CGRID2B", cgrid_size) < 0 ? 4 : c.get<size_t>("TYPE2B", 3 * i + 2);
        for (size_t j = 0; j < c1_cgrid_size; ++j) {
          c1.add("CGRID2B", cgrid[j + cgrid_size]);
        }
        cgrid_size += c1_cgrid_size;
      }
      if (sgrid_size < c.size("SGRID2B")) {
        size_t c1_sgrid_size = c.get<double>("SGRID2B", sgrid_size) < 0 ? 4 : c.get<size_t>("TYPE2B", 3 * i + 3);
        for (size_t j = 0; j < c1_sgrid_size; ++j) {
          c1.add("SGRID2B", sgrid[j + sgrid_size]);
        }
        sgrid_size += c1_sgrid_size;
      }
    } catch (const std::runtime_error& e) {
      std::cerr << "Grid Configuration Error: Ensure two grid generators are "
        << "used for both D2 descriptors or provide manual grids. Avoid "
        << "relying on a single auto grid generation for multiple descriptors. "
        << "Details: " << e.what() << std::endl;
    }

    c1.add("TYPE2B", type2b[3 * i + 1]);
    c1.add("RCTYPE2B", rctype2b[i]);
    c1.add("RCUT2B", rcut2b[i]);
    ds.push_back(CONFIG::factory<D2_Base, Config&>(c1.get<std::string>("TYPE2B"), c1));
    configs.push_back(c1);
    s += ds.back()->size();
  }

  cgrid.resize(cgrid_size);
  sgrid.resize(sgrid_size);
  c.remove("CGRID2B"); c.remove("SGRID2B");
  for (const auto &c1 : configs) {
    if (c1.exist("CGRID2B"))
      for (size_t i = 0; i < c1.size("CGRID2B"); ++i)
        c.add("CGRID2B", c1.get<double>("CGRID2B", i));
    if (c1.exist("SGRID2B"))
      for (size_t i = 0; i < c1.size("SGRID2B"); ++i)
        c.add("SGRID2B", c1.get<double>("SGRID2B", i));
  }
}

void D2_mJoin::calc_aed(const int Zj, const double rij, const double rij_sq, aed_type2 &aed, const double scale) {
  for (auto d : ds) d->calc_aed(Zj, rij, rij_sq, aed, scale);
}

void D2_mJoin::calc_dXijdri(const int Zj, const double rij, const double rij_sq, fd_type &fd_ij, const double scale) {
  for (auto d : ds) d->calc_dXijdri(Zj, rij, rij_sq, fd_ij, scale);
}

void D2_mJoin::calc_all(const int Zj, const double rij, const double rij_sq, aed_type2 &aed, fd_type &fd_ij, const double scale) {
  for (auto d : ds) d->calc_all(Zj, rij, rij_sq, aed, fd_ij, scale);
}

std::string D2_mJoin::label() {
  return lab;
}

void D2_mJoin::set_fidx(size_t fidx_) {
  ds[0]->set_fidx(fidx_);
  size_t s = 0;
  for (size_t i = 1; i < ds.size(); ++i) {
    s += ds[i - 1]->size();
    ds[i]->set_fidx(fidx_ + s);
  }
}

