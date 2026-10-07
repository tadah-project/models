#include <cstddef>
#include <tadah/models/descriptors/d2/d2_mjoin.h>
namespace tadah {
namespace models {

D2_mJoin::D2_mJoin(const tadah::core::Context& c)             // new, const version
  : D2_mJoin::D2_mJoin(const_cast<tadah::core::Context&>(c)) {}      // delegates to old one

D2_mJoin::D2_mJoin() {
  init();
}
D2_mJoin::~D2_mJoin() {
  for (auto d : ds) if (d) delete d;
}

D2_mJoin::D2_mJoin(tadah::core::Context &c) : D2_Base(c) {
  init();
  contexts = parse_context(c, "TYPE2B");
  expand_grids(c, contexts, "TYPE2B");

  for (auto &c1 : contexts) {
    ds.push_back(tadah::core::factory<D2_Base, tadah::core::Context&>(c1.get<std::string>("TYPE2B"), c1));
    s += ds.back()->size();
  }
}

void D2_mJoin::calc_aed(const int Zi, const int Zj, const double rij,
                        const double rij_sq, tadah::core::aed_type &aed, const double scale) const {
  for (auto d : ds) {
    d->calc_aed(Zi, Zj, rij, rij_sq, aed, scale);
  }
}

void D2_mJoin::calc_dXijdri(const int Zi, const int Zj, const double rij,
                            const double rij_sq, tadah::core::fd_type &fd_ij, const double scale) const {
  for (auto d : ds) {
    d->calc_dXijdri(Zi, Zj, rij, rij_sq, fd_ij, scale);
  }
}

void D2_mJoin::calc_all(const int Zi, const int Zj, const double rij,
                        const double rij_sq, tadah::core::aed_type &aed, tadah::core::fd_type &fd_ij, const double scale) const {
  for (auto d : ds) {
    d->calc_all(Zi, Zj, rij, rij_sq, aed, fd_ij, scale);
  }
}

bool D2_mJoin::applies_to_every_element_pair() const {
  for (auto d : ds)
    if (!d->applies_to_every_element_pair()) return false;
  return true;
}

std::string D2_mJoin::label() const {
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
std::string D2_mJoin::cutoff_fixed_by(std::size_t i) const {
  return i < ds.size() ? ds[i]->cutoff_fixed_by(0) : std::string();
}

bool D2_mJoin::reinit(const tadah::core::Context& ctx) {
  // Re-derive per-child sub-contexts from the updated master context.
  // Copy ctx because expand_grids may add default keys to the context.
  tadah::core::Context ctx_copy = ctx;
  contexts = D_mJoin::parse_context(ctx_copy, "TYPE2B");
  D_mJoin::expand_grids(ctx_copy, contexts, "TYPE2B");
  // Update the outer (container-level) fcut with the global RCUT2BMAX.
  // This is used by DescriptorsCalc::calc() as the 2B prefilter.
  if (fcut) fcut->set_rcut(ctx.get<double>("RCUT2BMAX"));
  bool any_changed = false;
  for (size_t i = 0; i < ds.size(); ++i)
    any_changed |= ds[i]->reinit(contexts[i]);
  // Recompute aggregated size from children (always, not just when changed).
  s = 0;
  for (auto d : ds) s += d->size();
  // Re-propagate fidx offsets so each child writes to the correct aed slot.
  if (any_changed) set_fidx(fidx);
  return any_changed;
}

} // namespace models
} // namespace tadah

