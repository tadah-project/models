#include <tadah/models/eam/pot_curve_edit.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#include <tadah/core/logger.h>
#include <tadah/models/descriptors/basis_dispatch.h>
#include <tadah/models/descriptors/d_base.h>
#include <tadah/models/descriptors/d2/d2_blip.h>
#include <tadah/models/descriptors/d2/d2_blip5.h>
#include <tadah/models/descriptors/d2/d2_knot5.h>
#include <tadah/models/descriptors/dm/dm_ream.h>

namespace tadah {
namespace models {
namespace potedit {

namespace {

// Support half-width of the bump families in the normalised coordinate
// (Blip: (-2,2), Blip5: (-3,3)); knots are one-sided.
double bump_halfwidth(BasisKind k) {
  return k == BasisKind::Blip5 ? 3.0 : 2.0;
}

std::string fmt(double v) {
  std::ostringstream os;
  os.precision(12);
  os << v;
  return os.str();
}

std::vector<double> read_grid(const core::Context &c, const std::string &key) {
  core::v_type v;
  D_Base::get_grid(c, key, v);
  return std::vector<double>(v.begin(), v.end());
}

std::vector<int> read_degs(const core::Context &c, const std::string &key,
                           std::size_t n_expected) {
  std::vector<int> d;
  if (!c.contains(key) || c.size(key) == 0) return d;
  d.resize(c.size(key));
  for (std::size_t i = 0; i < d.size(); ++i) {
    d[i] = static_cast<int>(c.get<double>(key, i));
    if (d[i] != 3 && d[i] != 5)
      throw std::runtime_error("pot edit: " + key + " entry " +
                               std::to_string(i + 1) + " is " +
                               std::to_string(d[i]) + " (allowed: 3, 5).");
  }
  if (d.size() != n_expected)
    throw std::runtime_error("pot edit: " + key + " has " +
                             std::to_string(d.size()) + " entries but the " +
                             "grid has " + std::to_string(n_expected) + ".");
  return d;
}

const char *d2_name(BasisKind k) {
  switch (k) {
    case BasisKind::Knot5: return "D2_Knot5";
    case BasisKind::Blip5: return "D2_Blip5";
    case BasisKind::Blip:  return "D2_Blip";
    default: break;
  }
  throw std::runtime_error("pot edit: unsupported two-body basis family.");
}

void add_all(core::Context &c, const std::string &key,
             const std::vector<double> &vals) {
  c.remove(key);
  for (double v : vals) c.add(key, v);
}

void add_degs(core::Context &c, const std::string &key,
              const std::vector<int> &deg) {
  c.remove(key);
  for (int d : deg) c.add(key, static_cast<double>(d));
}

}  // namespace

// ---------------------------------------------------------------------------
// CurveView evaluation
// ---------------------------------------------------------------------------

double CurveView::support_infinity() {
  return std::numeric_limits<double>::max();
}

double CurveView::value(double x) const {
  if (!right_active && rcut > 0.0 && x > rcut) return 0.0;  // Cut_Dummy
  double acc = 0.0;
  for (std::size_t k = 0; k < n(); ++k) {
    acc += coeff[k] * (right_active
        ? embedding_value(family, x, width(k), mu[k], degree(k))
        : radial_value_bare(family, x, width(k), mu[k], degree(k)));
  }
  return acc;
}

double CurveView::deriv(double x) const {
  if (!right_active && rcut > 0.0 && x > rcut) return 0.0;
  double acc = 0.0;
  for (std::size_t k = 0; k < n(); ++k) {
    acc += coeff[k] * (right_active
        ? embedding_deriv(family, x, width(k), mu[k], degree(k))
        : radial_deriv_bare(family, x, width(k), mu[k], degree(k)));
  }
  return acc;
}

double CurveView::deriv2(double x) const {
  if (!right_active && rcut > 0.0 && x > rcut) return 0.0;
  double acc = 0.0;
  for (std::size_t k = 0; k < n(); ++k) {
    acc += coeff[k] * (right_active
        ? embedding_deriv2(family, x, width(k), mu[k], degree(k))
        : radial_deriv2_bare(family, x, width(k), mu[k], degree(k)));
  }
  return acc;
}

std::pair<double, double> CurveView::support(std::size_t k) const {
  if (family == BasisKind::Knot5) {
    return right_active
        ? std::make_pair(mu[k], support_infinity())
        : std::make_pair(-support_infinity(), mu[k]);
  }
  const double w = bump_halfwidth(family) / width(k);
  return {mu[k] - w, mu[k] + w};
}

// ---------------------------------------------------------------------------
// Parse / rebuild
// ---------------------------------------------------------------------------

PotEAMView parse_pot_view(const core::Context &ctx) {
  PotEAMView v;
  v.ctx = ctx;

  if (!ctx.contains("TYPE2B") || !ctx.contains("TYPEMB"))
    throw std::runtime_error(
        "pot edit: the potential must carry TYPE2B and TYPEMB keys "
        "(only single-descriptor `tadah refit` potentials are supported).");

  // ---- two-body family ----
  const std::string t2 = ctx.get<std::string>("TYPE2B", 0);
  if (t2 == "D2_Knot5")      v.phi.family = BasisKind::Knot5;
  else if (t2 == "D2_Blip")  v.phi.family = BasisKind::Blip;
  else if (t2 == "D2_Blip5") v.phi.family = BasisKind::Blip5;
  else
    throw std::runtime_error(
        "pot edit: unsupported two-body descriptor '" + t2 +
        "' — only refit families (D2_Knot5, D2_Blip, D2_Blip5) can be "
        "edited (Gaussian/D2_BP and D2_mJoin are not supported).");

  // ---- many-body: DM_REAM only ----
  const std::string tm = ctx.get<std::string>("TYPEMB", 0);
  if (tm != "DM_REAM")
    throw std::runtime_error(
        "pot edit: unsupported many-body descriptor '" + tm +
        "' — only DM_REAM (as written by `tadah refit`) can be edited.");
  const BasisKind mb =
      basis_kind_from_code(ctx.get<int>("TYPEMB", 1));
  if (mb == BasisKind::Gaussian)
    throw std::runtime_error(
        "pot edit: the Gaussian DM_REAM basis is not supported.");

  if (!ctx.contains("ATOM") || ctx.size("ATOM") == 0)
    throw std::runtime_error("pot edit: the potential has no ATOM key.");
  v.element = ctx.get<std::string>("ATOM", 0);

  // ---- phi ----
  v.phi.right_active = false;
  v.phi.mu = read_grid(ctx, "CGRID2B");
  if (v.phi.mu.empty())
    throw std::runtime_error("pot edit: CGRID2B is empty.");
  if (v.phi.family == BasisKind::Knot5) {
    if (ctx.contains("SGRID2B") && ctx.size("SGRID2B") > 0)
      throw std::runtime_error(
          "pot edit: SGRID2B must be absent for the knot family "
          "(eta is a pure gauge).");
    v.phi.deg = read_degs(ctx, "KDEG2B", v.phi.mu.size());
  } else {
    v.phi.eta = read_grid(ctx, "SGRID2B");
    if (v.phi.eta.size() != v.phi.mu.size())
      throw std::runtime_error("pot edit: SGRID2B/CGRID2B size mismatch.");
  }
  v.phi.rcut = ctx.get<double>("RCUT2B");

  // ---- rho ----
  v.rho.family = mb;
  v.rho.right_active = false;
  v.rho.mu = read_grid(ctx, "CGRIDMB");
  v.rho.coeff = read_grid(ctx, "AMPGRIDMB");
  if (v.rho.mu.empty() || v.rho.coeff.size() != v.rho.mu.size())
    throw std::runtime_error(
        "pot edit: CGRIDMB/AMPGRIDMB must be non-empty and equal-sized.");
  if (mb == BasisKind::Knot5) {
    v.rho.deg = read_degs(ctx, "KDEGMB", v.rho.mu.size());
  } else {
    v.rho.eta = read_grid(ctx, "SGRIDMB");
    if (v.rho.eta.size() != v.rho.mu.size())
      throw std::runtime_error("pot edit: SGRIDMB/CGRIDMB size mismatch.");
  }
  v.rho.rcut = ctx.get<double>("RCUTMB");

  // ---- F (embedding) ----
  v.F.family = mb;
  v.F.right_active = (mb == BasisKind::Knot5);
  v.F.mu = read_grid(ctx, "CEMBFUNC");
  if (v.F.mu.empty())
    throw std::runtime_error("pot edit: CEMBFUNC is empty.");
  if (mb == BasisKind::Knot5) {
    v.F.deg = read_degs(ctx, "KDEGEMB", v.F.mu.size());
  } else {
    v.F.eta = read_grid(ctx, "SEMBFUNC");
    if (v.F.eta.size() != v.F.mu.size())
      throw std::runtime_error("pot edit: SEMBFUNC/CEMBFUNC size mismatch.");
  }
  v.F.rcut = 0.0;  // no r cutoff on the embedding

  // ---- weights: [phi | F] ----
  const std::size_t K = ctx.contains("WEIGHTS") ? ctx.size("WEIGHTS") : 0;
  const std::size_t n2b = v.phi.mu.size(), nEc = v.F.mu.size();
  if (K != n2b + nEc)
    throw std::runtime_error(
        "pot edit: WEIGHTS has " + std::to_string(K) + " entries, expected " +
        std::to_string(n2b) + " (CGRID2B) + " + std::to_string(nEc) +
        " (CEMBFUNC) = " + std::to_string(n2b + nEc) + ".");
  std::vector<double> w(K);
  for (std::size_t i = 0; i < K; ++i) w[i] = ctx.get<double>("WEIGHTS", i);
  v.phi.coeff.assign(w.begin(), w.begin() + n2b);
  v.F.coeff.assign(w.begin() + n2b, w.end());

  return v;
}

core::Context rebuild_pot_context(const PotEAMView &v) {
  core::Context c = v.ctx;
  const std::string el = v.element;
  const bool knot2b = (v.phi.family == BasisKind::Knot5);
  const bool knotmb = (v.rho.family == BasisKind::Knot5);

  // ---- two-body: <name> <nC> <nKdeg|nS> El El ----
  c.remove("TYPE2B");
  c.add("TYPE2B", d2_name(v.phi.family));
  c.add("TYPE2B", static_cast<int>(v.phi.mu.size()));
  c.add("TYPE2B", static_cast<int>(knot2b ? v.phi.deg.size()
                                          : v.phi.eta.size()));
  c.add("TYPE2B", el);
  c.add("TYPE2B", el);
  add_all(c, "CGRID2B", v.phi.mu);
  c.remove("SGRID2B");
  if (!knot2b) add_all(c, "SGRID2B", v.phi.eta);
  add_degs(c, "KDEG2B", v.phi.deg);

  // ---- many-body: DM_REAM <bcode> <nC> <nS> <nAmp> <nEc> <nEs> <nKd> <nKe> ----
  c.remove("TYPEMB");
  c.add("TYPEMB", "DM_REAM");
  c.add("TYPEMB", basis_code(v.rho.family));
  c.add("TYPEMB", static_cast<int>(v.rho.mu.size()));
  c.add("TYPEMB", static_cast<int>(knotmb ? 0 : v.rho.eta.size()));
  c.add("TYPEMB", static_cast<int>(v.rho.coeff.size()));
  c.add("TYPEMB", static_cast<int>(v.F.mu.size()));
  c.add("TYPEMB", static_cast<int>(knotmb ? 0 : v.F.eta.size()));
  c.add("TYPEMB", static_cast<int>(v.rho.deg.size()));
  c.add("TYPEMB", static_cast<int>(v.F.deg.size()));
  c.add("TYPEMB", el);
  c.add("TYPEMB", el);
  add_all(c, "CGRIDMB", v.rho.mu);
  c.remove("SGRIDMB");
  if (!knotmb) add_all(c, "SGRIDMB", v.rho.eta);
  add_all(c, "AMPGRIDMB", v.rho.coeff);
  add_degs(c, "KDEGMB", v.rho.deg);
  add_all(c, "CEMBFUNC", v.F.mu);
  c.remove("SEMBFUNC");
  if (!knotmb) add_all(c, "SEMBFUNC", v.F.eta);
  add_degs(c, "KDEGEMB", v.F.deg);

  // ---- weights + HPO coefficient-search seed (mirrors build_pot_context) ----
  std::vector<double> w(v.phi.coeff);
  w.insert(w.end(), v.F.coeff.begin(), v.F.coeff.end());
  add_all(c, "WEIGHTS", w);
  c.remove("FIXINDEX");
  c.remove("FIXWEIGHT");
  if (!w.empty()) {
    c.add("FIXINDEX", std::string("1-") + std::to_string(w.size()));
    for (double x : w) c.add("FIXWEIGHT", x);
  }

  // ---- RCUT: derived for knot channels, preserved for bump channels ----
  if (knot2b) {
    c.remove("RCUT2B");
    c.add("RCUT2B", *std::max_element(v.phi.mu.begin(), v.phi.mu.end()));
  } else {
    c.remove("RCUT2B");
    c.add("RCUT2B", v.phi.rcut);
  }
  if (knotmb) {
    c.remove("RCUTMB");
    c.add("RCUTMB", *std::max_element(v.rho.mu.begin(), v.rho.mu.end()));
  } else {
    c.remove("RCUTMB");
    c.add("RCUTMB", v.rho.rcut);
  }

  // ---- self-check: the rebuilt Context must construct cleanly ----
  try {
    core::Context probe = c;
    if (!probe.contains("VERBOSE")) probe.add("VERBOSE", 0);
    if (knot2b) { D2_Knot5 d(probe); }
    else if (v.phi.family == BasisKind::Blip5) { D2_Blip5 d(probe); }
    else { D2_Blip d(probe); }
    DM_REAM m(probe);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        std::string("pot edit: rebuilt potential failed descriptor "
                    "validation (internal error): ") + e.what());
  }
  return c;
}

// ---------------------------------------------------------------------------
// RCUT derivation (F2)
// ---------------------------------------------------------------------------

std::optional<double> derived_knot_rcut(const core::Context &ctx,
                                        const std::string &suffix) {
  try {
    if (suffix == "2B") {
      if (!ctx.contains("TYPE2B") || !ctx.contains("CGRID2B") ||
          ctx.size("CGRID2B") == 0)
        return std::nullopt;
      if (ctx.get<std::string>("TYPE2B", 0) != "D2_Knot5") return std::nullopt;
      const auto g = read_grid(ctx, "CGRID2B");
      return *std::max_element(g.begin(), g.end());
    }
    if (suffix == "MB") {
      if (!ctx.contains("TYPEMB") || !ctx.contains("CGRIDMB") ||
          ctx.size("CGRIDMB") == 0)
        return std::nullopt;
      if (ctx.get<std::string>("TYPEMB", 0) != "DM_REAM") return std::nullopt;
      if (basis_kind_from_code(ctx.get<int>("TYPEMB", 1)) != BasisKind::Knot5)
        return std::nullopt;
      const auto g = read_grid(ctx, "CGRIDMB");
      return *std::max_element(g.begin(), g.end());
    }
  } catch (const std::exception &) {
    return std::nullopt;  // malformed channel: not ours to police here
  }
  throw std::invalid_argument("derived_knot_rcut: suffix must be 2B or MB");
}

bool normalize_knot_rcut(core::Context &ctx, const std::string &who) {
  bool changed = false;
  for (const char *sfx : {"2B", "MB"}) {
    const auto derived = derived_knot_rcut(ctx, sfx);
    if (!derived) continue;
    const std::string key = std::string("RCUT") + sfx;
    const bool have = ctx.contains(key) && ctx.size(key) > 0;
    const double stated = have ? ctx.get<double>(key) : -1.0;
    if (have &&
        std::fabs(stated - *derived) <=
            1e-9 * std::max(1.0, std::fabs(*derived)))
      continue;
    std::ostringstream msg;
    msg << who << ": " << key << " ";
    if (!have) {
      msg << "is missing";
    } else if (stated < *derived) {
      msg << stated << " is SMALLER than the largest knot (" << fmt(*derived)
          << ") — knots beyond " << key << " were silently truncated";
    } else {
      msg << stated << " differs from the largest knot (" << fmt(*derived)
          << ") — the knot basis vanishes at its largest knot, so the "
          << "effective cutoff is the knot position";
    }
    msg << "; using the derived value " << fmt(*derived)
        << ". To change the effective cutoff, move the CGRID" << sfx
        << " knot positions.";
    core::Logger::instance().log(core::LogLevel::WARNING, msg.str());
    ctx.remove(key);
    ctx.add(key, *derived);
    changed = true;
  }
  return changed;
}

// ---------------------------------------------------------------------------
// ADD
// ---------------------------------------------------------------------------

AddReport add_basis(CurveView &c, const AddSpec &s, bool radial) {
  AddReport rep;
  const bool knot = (c.family == BasisKind::Knot5);
  if (c.n() == 0)
    throw std::runtime_error("pot edit: cannot edit an empty curve.");

  if (knot && s.eta)
    throw std::runtime_error(
        "pot edit: --eta is rejected for the knot family — the width is a "
        "pure gauge (pinned 1); the knot position is the only shape "
        "parameter.");
  if (!knot && s.deg)
    throw std::runtime_error(
        "pot edit: --deg only applies to the knot family (Knot5/Knot3); "
        "bump families have no per-function degree.");
  if (s.eta && *s.eta <= 0.0)
    throw std::runtime_error("pot edit: --eta must be positive.");
  rep.deg_used = s.deg.value_or(5);
  if (rep.deg_used != 3 && rep.deg_used != 5)
    throw std::runtime_error("pot edit: --deg must be 3 or 5.");
  if (radial && s.centre <= 0.0)
    throw std::runtime_error(
        "pot edit: a radial basis-function centre must be positive (got " +
        fmt(s.centre) + ").");

  if (!knot) {
    // width default: mean of the (up to) two nearest neighbours' widths.
    if (s.eta) {
      rep.eta_used = *s.eta;
    } else {
      std::vector<std::size_t> idx(c.n());
      for (std::size_t i = 0; i < idx.size(); ++i) idx[i] = i;
      std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
        return std::fabs(c.mu[a] - s.centre) < std::fabs(c.mu[b] - s.centre);
      });
      double acc = 0.0;
      std::size_t m = std::min<std::size_t>(2, idx.size());
      for (std::size_t i = 0; i < m; ++i) acc += c.eta[idx[i]];
      rep.eta_used = acc / static_cast<double>(m);
    }
    const double w = bump_halfwidth(c.family) / rep.eta_used;
    if (radial && c.rcut > 0.0) {
      if (s.centre - w >= c.rcut)
        throw std::runtime_error(
            "pot edit: the new basis function's support [" +
            fmt(s.centre - w) + ", " + fmt(s.centre + w) +
            "] lies entirely beyond RCUT = " + fmt(c.rcut) +
            " — it could never contribute (RCUT is fixed for blip "
            "potentials: the hard cutoff is part of the fitted model).");
      if (s.centre + w > c.rcut)
        rep.warnings.push_back(
            "the new basis function's support crosses RCUT = " + fmt(c.rcut) +
            " and is truncated there (RCUT is never changed for blip "
            "potentials).");
    }
  } else if (radial) {
    const double mx = *std::max_element(c.mu.begin(), c.mu.end());
    if (s.centre > mx)
      rep.warnings.push_back(
          "the new knot at " + fmt(s.centre) +
          " is beyond the current largest knot (" + fmt(mx) +
          ") — the effective cutoff (RCUT) grows to " + fmt(s.centre) +
          "; downstream HPO/MD neighbour lists must cover it.");
  }

  // Append; grids are positional and index-aligned with the coefficients,
  // sortedness is not required by any consumer.
  rep.index = c.n();
  c.mu.push_back(s.centre);
  c.coeff.push_back(0.0);
  if (!knot) {
    c.eta.push_back(rep.eta_used);
  } else if (rep.deg_used == 3 && c.deg.empty()) {
    c.deg.assign(rep.index, 5);  // materialise: prior entries were quintic
    c.deg.push_back(3);
  } else if (!c.deg.empty()) {
    c.deg.push_back(rep.deg_used);
  }
  // (deg stays empty when everything, including the new knot, is quintic)
  return rep;
}

// ---------------------------------------------------------------------------
// REMOVE
// ---------------------------------------------------------------------------

RemoveReport remove_basis(CurveView &c, double centre_request, bool radial) {
  RemoveReport rep;
  if (c.n() == 0)
    throw std::runtime_error("pot edit: cannot edit an empty curve.");
  if (c.n() == 1)
    throw std::runtime_error(
        "pot edit: cannot remove the last basis function of a curve.");

  // nearest centre; ties -> lower index
  std::size_t r = 0;
  for (std::size_t k = 1; k < c.n(); ++k)
    if (std::fabs(c.mu[k] - centre_request) < std::fabs(c.mu[r] - centre_request))
      r = k;
  rep.removed_idx = r;
  rep.removed_mu = c.mu[r];
  rep.removed_coeff = c.coeff[r];

  const bool knot = (c.family == BasisKind::Knot5);
  if (knot && radial) {
    const double mx = *std::max_element(c.mu.begin(), c.mu.end());
    if (c.mu[r] >= mx) {
      double next = -CurveView::support_infinity();
      for (std::size_t k = 0; k < c.n(); ++k)
        if (k != r && c.mu[k] > next) next = c.mu[k];
      rep.warnings.push_back(
          "removing the OUTERMOST knot: the curve on (" + fmt(next) + ", " +
          fmt(c.mu[r]) + "] cannot be represented by the remaining basis, "
          "and the effective cutoff (RCUT) shrinks to " + fmt(next) + ".");
    }
  }

  if (rep.removed_coeff != 0.0) {
    // ---- affected region ----
    const double r_lo = radial ? kCompensationRLoFrac * c.rcut : 0.0;
    double lo, hi;
    if (!knot) {
      const double w = bump_halfwidth(c.family) / c.width(r);
      lo = c.mu[r] - w;
      hi = c.mu[r] + w;
      if (radial) {
        lo = std::max(lo, r_lo);
        hi = std::min(hi, c.rcut);
      }
    } else if (!c.right_active) {
      lo = r_lo;
      hi = c.mu[r];
    } else {
      double x_lo = *std::min_element(c.mu.begin(), c.mu.end());
      double x_hi = *std::max_element(c.mu.begin(), c.mu.end());
      lo = c.mu[r];
      hi = x_hi + 0.25 * std::max(x_hi - x_lo, 1e-12);
    }

    if (hi > lo) {
      // ---- compensators: up to two remaining BFs nearest in centre ----
      std::vector<std::size_t> J;
      {
        std::vector<std::size_t> idx;
        for (std::size_t k = 0; k < c.n(); ++k)
          if (k != r) idx.push_back(k);
        std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
          return std::fabs(c.mu[a] - c.mu[r]) < std::fabs(c.mu[b] - c.mu[r]);
        });
        for (std::size_t i = 0; i < std::min<std::size_t>(2, idx.size()); ++i)
          J.push_back(idx[i]);
      }

      // ---- stacked value/1st/2nd rows on a dense sample of [lo, hi] ----
      constexpr std::size_t N = 512;
      auto bf = [&](std::size_t k, int order, double x) {
        if (c.right_active) {
          if (order == 0) return embedding_value (c.family, x, c.width(k), c.mu[k], c.degree(k));
          if (order == 1) return embedding_deriv (c.family, x, c.width(k), c.mu[k], c.degree(k));
          return embedding_deriv2(c.family, x, c.width(k), c.mu[k], c.degree(k));
        }
        if (order == 0) return radial_value_bare (c.family, x, c.width(k), c.mu[k], c.degree(k));
        if (order == 1) return radial_deriv_bare (c.family, x, c.width(k), c.mu[k], c.degree(k));
        return radial_deriv2_bare(c.family, x, c.width(k), c.mu[k], c.degree(k));
      };

      const std::size_t nj = J.size();
      // per-order data: t[m][i] (lost contribution), A[m][i][j]
      std::array<std::vector<double>, 3> t;
      std::array<std::vector<std::array<double, 2>>, 3> A;
      std::array<double, 3> wm{0.0, 0.0, 0.0};
      for (int m = 0; m < 3; ++m) {
        t[m].resize(N);
        A[m].resize(N);
        double ss = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
          const double x =
              lo + (hi - lo) * static_cast<double>(i) / (N - 1);
          t[m][i] = c.coeff[r] * bf(r, m, x);
          for (std::size_t j = 0; j < nj; ++j) A[m][i][j] = bf(J[j], m, x);
          ss += t[m][i] * t[m][i];
        }
        wm[m] = std::sqrt(ss / N);
      }
      const double wmax = std::max({wm[0], wm[1], wm[2]});
      for (int m = 0; m < 3; ++m)
        wm[m] = (wm[m] > 1e-14 * wmax && wm[m] > 0.0) ? 1.0 / wm[m] : 0.0;

      // ---- normal equations over the compensators ----
      auto solve = [&](std::size_t nuse, std::array<double, 2> &dc) -> bool {
        double G[2][2] = {{0, 0}, {0, 0}}, h[2] = {0, 0};
        for (int m = 0; m < 3; ++m) {
          if (wm[m] == 0.0) continue;
          const double w2 = wm[m] * wm[m];
          for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t a = 0; a < nuse; ++a) {
              h[a] += w2 * A[m][i][a] * t[m][i];
              for (std::size_t b = 0; b < nuse; ++b)
                G[a][b] += w2 * A[m][i][a] * A[m][i][b];
            }
          }
        }
        if (nuse == 2) {
          const double det = G[0][0] * G[1][1] - G[0][1] * G[1][0];
          if (std::fabs(det) <= 1e-12 * std::fabs(G[0][0] * G[1][1]))
            return false;  // ill-conditioned: caller retries with one
          dc[0] = (h[0] * G[1][1] - h[1] * G[0][1]) / det;
          dc[1] = (h[1] * G[0][0] - h[0] * G[1][0]) / det;
          return true;
        }
        if (G[0][0] <= 0.0) return false;
        dc[0] = h[0] / G[0][0];
        dc[1] = 0.0;
        return true;
      };

      std::array<double, 2> dc{0.0, 0.0};
      std::size_t nuse = nj;
      bool ok = nuse > 0 && solve(nuse, dc);
      if (!ok && nuse == 2) {
        nuse = 1;
        ok = solve(nuse, dc);
        if (ok)
          rep.warnings.push_back(
              "compensation system ill-conditioned with two neighbours; "
              "adjusted only the nearest one.");
      }

      // ---- residuals (stacked weighted + value-only RMS) ----
      auto resid = [&](bool with_dc) {
        double stacked = 0.0, val = 0.0;
        for (int m = 0; m < 3; ++m) {
          for (std::size_t i = 0; i < N; ++i) {
            double d = -t[m][i];
            if (with_dc)
              for (std::size_t j = 0; j < nuse; ++j) d += dc[j] * A[m][i][j];
            stacked += wm[m] * wm[m] * d * d;
            if (m == 0) val += d * d;
          }
        }
        return std::make_pair(std::sqrt(stacked / (3.0 * N)),
                              std::sqrt(val / N));
      };
      std::tie(rep.resid_plain, rep.vresid_plain) = resid(false);
      if (ok) {
        std::tie(rep.resid_comp, rep.vresid_comp) = resid(true);
        for (std::size_t j = 0; j < nuse; ++j) {
          c.coeff[J[j]] += dc[j];
          rep.adjusted.push_back(J[j]);
          rep.delta_c.push_back(dc[j]);
        }
        rep.compensated = true;
      } else {
        rep.resid_comp = rep.resid_plain;
        rep.vresid_comp = rep.vresid_plain;
        rep.warnings.push_back(
            "no usable compensation system — plain removal applied.");
      }
    } else {
      rep.warnings.push_back(
          "the removed basis function only acts below the sampled window "
          "(r < " + fmt(r_lo) + ") — plain removal applied.");
    }
  }

  // ---- erase index r from every parallel array ----
  c.mu.erase(c.mu.begin() + r);
  c.coeff.erase(c.coeff.begin() + r);
  if (!c.eta.empty()) c.eta.erase(c.eta.begin() + r);
  if (!c.deg.empty()) {
    c.deg.erase(c.deg.begin() + r);
    // canonical form: an all-quintic degree array collapses to "absent"
    if (std::all_of(c.deg.begin(), c.deg.end(),
                    [](int d) { return d == 5; }))
      c.deg.clear();
  }
  return rep;
}

}  // namespace potedit
}  // namespace models
}  // namespace tadah
