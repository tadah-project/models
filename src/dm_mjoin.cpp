#include <cstddef>
#include <tadah/models/descriptors/dm/dm_mjoin.h>

DM_mJoin::DM_mJoin() {
  init();
}
DM_mJoin::~DM_mJoin() {
  for (auto d : ds) if (d) delete d;
}

DM_mJoin::DM_mJoin(Config &c) : DM_Base(c) {
  init();
  configs = parse_config(c, "TYPEMB");
  c.remove("CGRIDMB"); c.remove("SGRIDMB");
  c.remove("CEMBFUNC"); c.remove("SEMBFUNC");

  for (auto &c1 : configs) {
    ds.push_back(CONFIG::factory<DM_Base, Config&>(c1.get<std::string>("TYPEMB"), c1));
    s += ds.back()->size();
    rhoisize += ds.back()->rhoi_size();

    // update main grid with expanded grid
    if (c1.exist("CGRIDMB"))
      for (size_t i = 0; i < c1.size("CGRIDMB"); ++i)
        c.add("CGRIDMB", c1.get<double>("CGRIDMB", i));
    if (c1.exist("SGRIDMB"))
      for (size_t i = 0; i < c1.size("SGRIDMB"); ++i)
        c.add("SGRIDMB", c1.get<double>("SGRIDMB", i));
    if (c1.exist("CEMBFUNC"))
      for (size_t i = 0; i < c1.size("CEMBFUNC"); ++i)
        c.add("CEMBFUNC", c1.get<double>("CEMBFUNC", i));
    if (c1.exist("SEMBFUNC"))
      for (size_t i = 0; i < c1.size("SEMBFUNC"); ++i)
        c.add("SEMBFUNC", c1.get<double>("SEMBFUNC", i));
  }
}


void DM_mJoin::calc_aed(rho_type& rho, aed_type &aed) {
  size_t rho_fidx = 0;
  for (auto d : ds) {
    size_t rho_size = 2 * d->rhoi_size();
    rho_type rho_ptr(&rho[rho_fidx], rho_size, 1);
    d->calc_aed(rho_ptr, aed);
    rho_fidx += rho_size;
  }
}
void DM_mJoin::calc_dXijdri_dXjidri(
  const int Zi,
  const int Zj,
  const double rij,
  const double rij_sq,
  const Vec3d &vec_ij,
  rho_type& rhoi,
  rho_type& rhoj,
  fd_type &fd_ij,
  const double scale) {
  size_t rho_fidx = 0;
  for (auto d : ds) {
    if (rij > d->get_rcut()) return;
    size_t rho_size = 2 * d->rhoi_size();
    rho_type rhoi_ptr(&rhoi[rho_fidx],rho_size);
    rho_type rhoj_ptr(&rhoj[rho_fidx],rho_size);
    d->calc_dXijdri_dXjidri(Zi,Zj,rij,rij_sq,vec_ij,rhoi_ptr,rhoj_ptr,fd_ij,scale); 
    rho_fidx += rho_size;
  }
}
void DM_mJoin::calc_dXijdri(
  const int Zj,
  const double rij,
  const double rij_sq,
  const Vec3d &vec_ij,
  rho_type& rhoi,
  fd_type &fd_ij, 
  const double scale) {
  size_t rho_fidx = 0;
  for (auto d : ds) { 
    if (rij > d->get_rcut()) return;
    size_t rho_size = 2 * d->rhoi_size();
    rho_type rhoi_ptr(&rhoi[rho_fidx],rho_size);
    d->calc_dXijdri(Zj,rij,rij_sq,vec_ij,rhoi_ptr,fd_ij,scale);
    rho_fidx += rho_size;
  }
}
void DM_mJoin::init_rhoi(rho_type& rhoi) {
  rhoi.resize(2*rhoisize);
  rhoi.set_zero();
}
void DM_mJoin::calc_rho(
  const int Zj,
  const double rij,
  const double rij_sq,
  const Vec3d &vec_ij,
  rho_type& rho,
  const double scale) {
  size_t rho_fidx = 0;
  for (auto d : ds) {
    if (rij > d->get_rcut()) return;
    size_t rho_size = 2 * d->rhoi_size();
    rho_type rho_ptr(&rho[rho_fidx],rho_size);
    d->calc_rho(Zj,rij,rij_sq,vec_ij,rho_ptr,scale);
    rho_fidx += rho_size;
  }
}
std::string DM_mJoin::label() {
  return lab;
}

void DM_mJoin::set_fidx(size_t fidx_) {
  ds[0]->set_fidx(fidx_);
  size_t s = 0;
  for (size_t i = 1; i < ds.size(); ++i) {
    s += ds[i - 1]->size();
    ds[i]->set_fidx(fidx_ + s);
  }
}

void DM_mJoin::init() {
  // not really sure about this
  keys.push_back("TYPEMB");
  keys.push_back("RCTYPEMB");
  keys.push_back("RCUTMB");
  nparams=std::numeric_limits<size_t>::max();
}

