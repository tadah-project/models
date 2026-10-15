#include <cstddef>
#include <memory>
#include <tadah/models/descriptors/dm/dm_mjoin.h>
#include <tadah/models/descriptors/dm/embed_spec.h>
namespace tadah {
namespace models {

DM_mJoin::DM_mJoin(const tadah::core::Context& c):DM_Base(c) {
  init();auto copy=c;build_children(copy);
}

DM_mJoin::DM_mJoin() {
  init();
}
DM_mJoin::~DM_mJoin() {
  for (auto d : ds) if (d) delete d;
}

DM_mJoin::DM_mJoin(tadah::core::Context &c) : DM_Base(c) {
  init();
  build_children(c);
}
void DM_mJoin::build_children(tadah::core::Context &c) {
  contexts = parse_context(c, "TYPEMB");
  expand_grids(c, contexts, "TYPEMB");
  if(!c.blocks("DESCRIPTOR").empty())
    set_fcut(new Cut_Dummy(c.get<double>("RCUTMBMAX")),true);

  std::vector<std::unique_ptr<DM_Base>> pending;
  for (auto &c1 : contexts) {
    pending.emplace_back(tadah::core::factory<DM_Base, tadah::core::Context&>(c1.get<std::string>("TYPEMB"), c1));
    s += pending.back()->size();
    rhoisize += pending.back()->rhoi_size();
    rhoipsize += pending.back()->rhoip_size();
  }
  ds.reserve(pending.size());for(auto &d:pending)ds.push_back(d.release());
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
  std::optional<EmbF0Mode> mode;
  for(auto d:ds) if(auto current=d->embf0_mode()) {
    if(mode && *mode!=*current)return std::nullopt;
    mode=current;
  }
  return mode;
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
std::string DM_mJoin::cutoff_fixed_by(std::size_t i) const {
  return i < ds.size() ? ds[i]->cutoff_fixed_by(0) : std::string();
}

bool DM_mJoin::reinit(const tadah::core::Context& ctx) {
  if(!ctx.blocks("DESCRIPTOR").empty()) {
    auto copy=ctx;resolve_embed_blocks(copy);
    DM_mJoin candidate(copy);
    candidate.set_fidx(fidx);candidate.set_rfidx(rfidx);candidate.set_drfidx(rfidx+candidate.rhoisize);
    const bool changed=s!=candidate.s || rhoisize!=candidate.rhoisize || rhoipsize!=candidate.rhoipsize;
    std::unique_ptr<Cut_Base> cutoff(new Cut_Dummy(copy.get<double>("RCUTMBMAX")));
    ds.swap(candidate.ds);contexts.swap(candidate.contexts);
    std::swap(s,candidate.s);std::swap(rhoisize,candidate.rhoisize);std::swap(rhoipsize,candidate.rhoipsize);
    drfidx=rfidx+rhoisize;set_fcut(cutoff.release(),true);
    return changed;
  }
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
  // The sizes are the children's sums, as the constructor adds them: a
  // child's descriptor or density may have changed size.
  const size_t old_s = s, old_rhoisize = rhoisize, old_rhoipsize = rhoipsize;
  s = 0;
  rhoisize = 0;
  rhoipsize = 0;
  for (auto d : ds) {
    s += d->size();
    rhoisize += d->rhoi_size();
    rhoipsize += d->rhoip_size();
  }
  // Place the children again from this descriptor's first indices, with
  // F'(rho) after the whole density, as the constructor places them (a
  // child's reinit may also have moved its own F'(rho)).
  set_fidx(fidx);
  set_rfidx(rfidx);
  set_drfidx(rfidx + rhoisize);
  return any_changed || s != old_s || rhoisize != old_rhoisize ||
         rhoipsize != old_rhoipsize;
}


} // namespace models
} // namespace tadah
