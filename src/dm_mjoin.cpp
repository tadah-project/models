#include <cstddef>
#include <tadah/models/descriptors/dm/dm_mjoin.h>
namespace tadah {
namespace models {

DM_mJoin::DM_mJoin(const tadah::core::Context& c)             // new, const version
  : DM_mJoin::DM_mJoin(const_cast<tadah::core::Context&>(c)) {}      // delegates to old one

DM_mJoin::DM_mJoin() {
  init();
}
DM_mJoin::~DM_mJoin() {
  for (auto d : ds) if (d) delete d;
}

DM_mJoin::DM_mJoin(tadah::core::Context &c) : DM_Base(c) {
  init();
  contexts = parse_context(c, "TYPEMB");
  expand_grids(c, contexts, "TYPEMB");

  for (auto &c1 : contexts) {
    ds.push_back(tadah::core::factory<DM_Base, tadah::core::Context&>(c1.get<std::string>("TYPEMB"), c1));
    s += ds.back()->size();
    rhoisize += ds.back()->rhoi_size();
    rhoipsize += ds.back()->rhoip_size();
  }
  set_rfidx(0);
  set_drfidx(rhoisize);
}

void DM_mJoin::calc_aed(tadah::core::rho_type& rho, tadah::core::aed_type &aed) const {
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
  const double scale) const {
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
  const double scale) const {
  for (auto d : ds) { 
    d->calc_dXijdri(Zi,Zj,rij,rij_sq,vec_ij,rhoi,fd_ij,scale);
  }
}
void DM_mJoin::init_rhoi(tadah::core::rho_type& rhoi) {
  // this should not call child descriptors
  rhoi.resize(rhoisize + rhoipsize);
  rhoi.set_zero();
}
void DM_mJoin::calc_rho(
  const int Zi,
  const int Zj,
  const double rij,
  const double rij_sq,
  const tadah::core::Vec3d &vec_ij,
  tadah::core::rho_type& rho,
  const double scale) const {
  for (auto d : ds) {
      d->calc_rho(Zi,Zj,rij,rij_sq,vec_ij,rho,scale);
  }
}
bool DM_mJoin::applies_to_every_element_pair() const {
  for (auto d : ds)
    if (!d->applies_to_every_element_pair()) return false;
  return true;
}

std::string DM_mJoin::label() const {
  return lab;
}

std::optional<EmbF0Mode> DM_mJoin::embf0_mode() const {
  for (auto d : ds)
    if (auto m = d->embf0_mode()) return m;
  return std::nullopt;
}

std::optional<double>
DM_mJoin::raw_embedding_f0(const tadah::core::t_type &w_full) const {
  std::optional<double> acc;
  for (auto d : ds)
    if (auto v = d->raw_embedding_f0(w_full)) acc = acc.value_or(0.0) + *v;
  return acc;
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
    // s += ds[i - 1]->rhoip_size();
    ds[i]->set_rfidx(rfidx_ + s);
  }
}
void DM_mJoin::set_drfidx(size_t drfidx_) {
  drfidx = drfidx_;
  ds[0]->set_drfidx(drfidx_);
  size_t s = 0;
  for (size_t i = 1; i < ds.size(); ++i) {
    // s += ds[i - 1]->rhoi_size();
    s += ds[i - 1]->rhoip_size();
    ds[i]->set_drfidx(drfidx_ + s);
  }
}

void DM_mJoin::init() {
  // not really sure about this
  keys.push_back("TYPEMB");
  keys.push_back("RCTYPEMB");
  keys.push_back("RCUTMB");
  nparams=std::numeric_limits<size_t>::max();
}
bool DM_mJoin::reinit(const tadah::core::Context& ctx) {
  // Re-derive per-child sub-contexts from the updated master context.
  // Copy ctx because expand_grids may add default keys to the context.
  // Note: contexts[0] is the first child (DM_mJoin itself is filtered by parse_context).
  tadah::core::Context ctx_copy = ctx;
  contexts = D_mJoin::parse_context(ctx_copy, "TYPEMB");
  D_mJoin::expand_grids(ctx_copy, contexts, "TYPEMB");
  // Update the outer (container-level) fcut with the global RCUTMBMAX.
  // This is used by DescriptorsCalc::calc_rho() as the MB prefilter.
  if (fcut) fcut->set_rcut(ctx.get<double>("RCUTMBMAX"));
  bool any_changed = false;
  for (size_t i = 0; i < ds.size(); ++i)
    any_changed |= ds[i]->reinit(contexts[i]);
  // TODO : why would this be required? I don't see the point.
  // Recompute aggregated sizes from children (always, not just when changed).
  // s = 0; rhoisize = 0;
  // for (auto d : ds) { s += d->size(); rhoisize += d->rhoi_size(); }
  // Re-propagate rfidx offsets so each child reads from the correct rho slice.
  // if (any_changed) set_rfidx(rfidx); // TODO this looks wrong its either zero or caller should pass in the correct rfidx initial value.
  return any_changed;
}


} // namespace models
} // namespace tadah
