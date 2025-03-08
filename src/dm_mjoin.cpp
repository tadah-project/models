#include <cstddef>
#include <tadah/models/descriptors/dm/dm_mjoin.h>
namespace tadah {
namespace models {


DM_mJoin::DM_mJoin() {
  init();
}
DM_mJoin::~DM_mJoin() {
  for (auto d : ds) if (d) delete d;
}

DM_mJoin::DM_mJoin(tadah::core::Config &c) : DM_Base(c) {
  init();
  configs = parse_config(c, "TYPEMB");
  expand_grids(c, configs, "TYPEMB");

  for (auto &c1 : configs) {
    ds.push_back(tadah::core::factory<DM_Base, tadah::core::Config&>(c1.get<std::string>("TYPEMB"), c1));
    s += ds.back()->size();
    rhoisize += ds.back()->rhoi_size();
  }
  set_rfidx(0);
}

void DM_mJoin::calc_aed(tadah::core::rho_type& rho, tadah::core::aed_type &aed) {
  for (auto d : ds) {
    d->calc_aed(rho, aed);
  }
}
void DM_mJoin::calc_dXijdri_dXjidri(
  const int Zi,
  const int Zj,
  const double rij,
  const double rij_sq,
  const tadah::core::Vec3d &vec_ij,
  tadah::core::rho_type& rhoi,
  tadah::core::rho_type& rhoj,
  tadah::core::fd_type &fd_ij,
  const double scale) {
  for (auto d : ds) {
    d->calc_dXijdri_dXjidri(Zi,Zj,rij,rij_sq,vec_ij,rhoi,rhoj,fd_ij,scale); 
  }
}
void DM_mJoin::calc_dXijdri(
  const int Zi,
  const int Zj,
  const double rij,
  const double rij_sq,
  const tadah::core::Vec3d &vec_ij,
  tadah::core::rho_type& rhoi,
  tadah::core::fd_type &fd_ij, 
  const double scale) {
  for (auto d : ds) { 
    d->calc_dXijdri(Zi,Zj,rij,rij_sq,vec_ij,rhoi,fd_ij,scale);
  }
}
void DM_mJoin::init_rhoi(tadah::core::rho_type& rhoi) {
  rhoi.resize(2*rhoisize);
  rhoi.set_zero();
}
void DM_mJoin::calc_rho(
  const int Zi,
  const int Zj,
  const double rij,
  const double rij_sq,
  const tadah::core::Vec3d &vec_ij,
  tadah::core::rho_type& rho,
  const double scale) {
  for (auto d : ds) {
      d->calc_rho(Zi,Zj,rij,rij_sq,vec_ij,rho,scale);
  }
}
std::string DM_mJoin::label() {
  return lab;
}

void DM_mJoin::set_fidx(size_t fidx_) {
  fidx = fidx_;
  ds[0]->set_fidx(fidx_);
  size_t s = 0;
  for (size_t i = 1; i < ds.size(); ++i) {
    s += ds[i - 1]->size();
    ds[i]->set_fidx(fidx_ + s);
  }
}
void DM_mJoin::set_rfidx(size_t rfidx_) {
  rfidx = rfidx_;
  ds[0]->set_rfidx(rfidx_);
  size_t s = 0;
  for (size_t i = 1; i < ds.size(); ++i) {
    s += ds[i - 1]->rhoi_size();
    s += ds[i - 1]->rhoip_size();
    ds[i]->set_rfidx(rfidx_ + s);
  }
}

void DM_mJoin::init() {
  // not really sure about this
  keys.push_back("TYPEMB");
  keys.push_back("RCTYPEMB");
  keys.push_back("RCUTMB");
  nparams=std::numeric_limits<size_t>::max();
}


} // namespace models
} // namespace tadah

