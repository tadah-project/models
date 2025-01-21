#include <cstddef>
#include <tadah/models/descriptors/d2/d2_mjoin.h>

D2_mJoin::D2_mJoin() {
  init();
}
D2_mJoin::~D2_mJoin() {
  for (auto d : ds) if (d) delete d;
}

D2_mJoin::D2_mJoin(Config &c) : D2_Base(c) {
  init();
  configs = parse_config(c, "TYPE2B");
  expand_grids(c, configs, "TYPE2B");

  for (auto &c1 : configs) {
    ds.push_back(CONFIG::factory<D2_Base, Config&>(c1.get<std::string>("TYPE2B"), c1));
    s += ds.back()->size();
  }
}

void D2_mJoin::calc_aed(const int Zi, const int Zj, const double rij,
                        const double rij_sq, aed_type &aed, const double scale) {
  for (auto d : ds) {
    d->calc_aed(Zi, Zj, rij, rij_sq, aed, scale);
  }
}

void D2_mJoin::calc_dXijdri(const int Zi, const int Zj, const double rij,
                            const double rij_sq, fd_type &fd_ij, const double scale) {
  for (auto d : ds) {
    d->calc_dXijdri(Zi, Zj, rij, rij_sq, fd_ij, scale);
  }
}

void D2_mJoin::calc_all(const int Zi, const int Zj, const double rij,
                        const double rij_sq, aed_type &aed, fd_type &fd_ij, const double scale) {
  for (auto d : ds) {
    d->calc_all(Zi, Zj, rij, rij_sq, aed, fd_ij, scale);
  }
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

void D2_mJoin::init() {
  keys.push_back("TYPE2B");
  keys.push_back("RCTYPE2B");
  keys.push_back("RCUT2B");
  nparams=std::numeric_limits<size_t>::max();
}
