#include <tadah/models/descriptors/density_reach.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace tadah {
namespace models {

namespace {

using tadah::core::stats::DataStats;

/* @p v with up to 6 significant digits ("2.85788", "0.2", "30"). */
std::string num(double v) {
  std::ostringstream os;
  os << v;
  return os.str();
}

/* "KEY(n)" for value @p i of @p key of @p s, n its 1-based number in the
 * configuration, as an OPTIM line names an entry (SeedView::entry_number:
 * the entries of a joined line's descriptor follow those of the
 * descriptors before it). */
std::string entry(const SeedView &s, const std::string &key, std::size_t i) {
  return key + "(" + std::to_string(s.entry_number(key, i)) + ")";
}

/* The seed's centre of entry @p arg_pos (-1: the first) of CGRIDMB, when
 * the seed has one. */
bool seed_centre(const SeedView &s, int arg_pos, double &c) {
  if (!s.contains("CGRIDMB")) return false;
  const std::size_t i = arg_pos < 0 ? 0 : static_cast<std::size_t>(arg_pos);
  if (i >= s.size("CGRIDMB")) return false;
  c = s.get<double>("CGRIDMB", i);
  return true;
}

/* The neighbour shell nearest to the distance @p c: the peak of the
 * distances of the first element pair of @p d nearest to it, its position
 * @p mu and its width @p w, the FWHM (0.2 times the position for a peak
 * without one); @p what names the width. False when the distances show no
 * peak. */
bool nearest_shell(const DataStats &d, double c, double &mu, double &w,
                   std::string &what) {
  const auto &peaks = d.per_pair[0].r_peaks;
  if (peaks.empty()) return false;
  const auto *best = &peaks.front();
  for (const auto &p : peaks)
    if (std::abs(p.position - c) < std::abs(best->position - c)) best = &p;
  mu = best->position;
  if (best->fwhm > 0.0) {
    w = best->fwhm;
    what = "the FWHM (" + num(w) + " Å) of the neighbour shell nearest to "
           "it (" + num(mu) + " Å)";
  } else {
    w = 0.2 * mu;
    what = num(w) + " Å, 0.2 times the distance of the neighbour shell "
           "nearest to it (" + num(mu) + " Å), which has no FWHM";
  }
  return true;
}

/* The first two neighbour shells of @p d: the first two peaks of the
 * distances of its first element pair, or, without them, the 10th
 * percentile and the median; @p what names them. */
void shells(const DataStats &d, double &r1, double &r2, std::string &what) {
  const auto &pp = d.per_pair[0];
  if (pp.r_peaks.size() >= 2) {
    r1 = pp.r_peaks[0].position;
    r2 = pp.r_peaks[1].position;
    what = "from the first neighbour shell (" + num(r1) +
           " Å) to the second (" + num(r2) + " Å)";
  } else {
    r1 = pp.r_percentiles[tadah::core::stats::IDX_P10];
    r2 = pp.r_percentiles[tadah::core::stats::IDX_P50];
    what = "from the 10th percentile of the neighbour distances (" +
           num(r1) + " Å) to their median (" + num(r2) + " Å)";
  }
}

}  // namespace

NeighbourRange neighbour_range(const DataStats &d, double rcut) {
  NeighbourRange r;
  if (d.per_pair.empty()) return r;
  r.lo = std::numeric_limits<double>::infinity();
  double top = 0.0;
  for (const auto &pp : d.per_pair) {
    r.lo = std::min(r.lo, pp.r_min);
    top = std::max(top, pp.r_max);
  }
  r.hi = rcut > 0.0 ? std::min(top, rcut) : top;
  r.empty = !(r.lo <= r.hi);
  return r;
}

double gaussian_reach_width(double dist) {
  if (dist == 0.0) return std::numeric_limits<double>::infinity();
  return std::log(1.0 / kDensityReachFloor) / (dist * dist);
}

HealthVerdict radial_reach_check(const DataStats &d, const SeedView &s,
                                 BasisKind kind) {
  if (kind != BasisKind::Gaussian && kind != BasisKind::Blip)
    throw std::logic_error(
        "radial_reach_check: a Gaussian or blip radial basis, not " +
        std::string(basis_name(kind)));
  if (!s.contains("CGRIDMB") || d.n_atoms == 0) return HealthOk{};
  const bool has_cut = s.contains("RCUTMB");
  const double rcut = has_cut ? s.get<double>("RCUTMB") : 0.0;
  const NeighbourRange range = neighbour_range(d, rcut);
  const std::string within =
      has_cut ? "within RCUTMB = " + num(rcut) : "of the dataset";
  if (range.empty) {
    if (!d.per_pair.empty() && has_cut)
      return HealthFail{"RCUTMB = " + num(rcut) +
                        " is below the smallest neighbour distance of the "
                        "dataset (" + num(range.lo) +
                        " Å): the density sums over no neighbour, so ρ is 0 "
                        "for every atom."};
    return HealthFail{"no pair of atoms of the dataset is " + within +
                      ": ρ is 0 for every atom."};
  }
  const std::string distances = "every neighbour distance " + within + " (" +
                                num(range.lo) + " to " + num(range.hi) +
                                " Å)";
  const std::size_t n = s.size("CGRIDMB");
  const std::size_t widths = s.contains("SGRIDMB") ? s.size("SGRIDMB") : 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double c = s.get<double>("CGRIDMB", i);
    if (i >= widths)
      return HealthFail{entry(s, "CGRIDMB", i) + " = " + num(c) +
                        " has no width " + entry(s, "SGRIDMB", i) +
                        ": the reach of its channel of ρ cannot be judged."};
    const double eta = s.get<double>("SGRIDMB", i);
    if (!(eta > 0.0)) continue;  // no decay: refused by other rules
    // The distance from the centre to the nearest neighbour distance.
    const double dist =
        c < range.lo ? range.lo - c : (c > range.hi ? c - range.hi : 0.0);
    const std::string channel = entry(s, "CGRIDMB", i) + " = " + num(c) +
                                " with " + entry(s, "SGRIDMB", i) + " = " +
                                num(eta);
    if (kind == BasisKind::Gaussian) {
      if (!(eta > gaussian_reach_width(dist))) continue;
      return HealthFail{channel +
                        ": the Gaussian exp(-SGRIDMB (r - CGRIDMB)^2) is "
                        "below 1e-6 at " + distances +
                        ", so ρ is ~0 for every atom (a smaller SGRIDMB, "
                        "or a centre nearer the neighbours, reaches them)."};
    }
    const double half = 2.0 / eta;  // b(eta (r - c)) is 0 for |r - c| >= 2/eta
    if (dist < half) continue;
    return HealthFail{channel + ": the support of the blip, |r - CGRIDMB| < "
                      "2/SGRIDMB = " + num(half) + " Å, holds no neighbour "
                      "distance " + within + " (" + num(range.lo) + " to " +
                      num(range.hi) + " Å), so ρ is 0 for every atom (a "
                      "smaller SGRIDMB widens the support)."};
  }
  return HealthOk{};
}

BoundSuggestion gaussian_centre_suggestion(const DataStats &d,
                                           const SeedView &s, int arg_pos) {
  if (d.per_pair.empty()) return std::make_pair(0.0, 0.0);
  double c = 0.0;
  const bool seeded = seed_centre(s, arg_pos, c);
  if (seeded && c == 0.0)
    return BoundSuggestion::no_row(
        "Gaussian exp(-SGRIDMB r^2) centred at the origin (the seed's "
        "CGRIDMB 0): its width SGRIDMB, not its centre, sets its reach");
  const auto &pp = d.per_pair[0];
  double mu = 0.0, w = 0.0;
  std::string width;
  if (seeded && nearest_shell(d, c, mu, w, width)) {
    // About the seed's centre, within the neighbour distances the channel
    // sums over, unless the seed's centre lies outside them: the range
    // always holds it.
    const double rcut = s.contains("RCUTMB") ? s.get<double>("RCUTMB") : 0.0;
    const NeighbourRange range = neighbour_range(d, rcut);
    return BoundSuggestion(
        std::max(c - w, std::min(c, range.lo)),
        std::min(c + w, std::max(c, range.hi)),
        "Gaussian exp(-SGRIDMB (r - CGRIDMB)^2) centred on the seed's centre "
        "(" + num(c) + " Å): +- " + width + ", within the neighbour "
        "distances " + num(range.lo) + " to " + num(range.hi) +
        " Å unless the seed's centre lies outside them");
  }
  if (!pp.r_peaks.empty()) {
    mu = pp.r_peaks[0].position;
    w = pp.r_peaks[0].fwhm > 0.0 ? pp.r_peaks[0].fwhm : 0.2 * mu;
    return BoundSuggestion(
        std::max(pp.r_min, mu - w), mu + w,
        "Gaussian exp(-SGRIDMB (r - CGRIDMB)^2) centred on the first "
        "neighbour shell: its peak (" + num(mu) + " Å) +- " + num(w) +
        " Å, from r_min (" + num(pp.r_min) + " Å)");
  }
  using tadah::core::stats::IDX_P10;
  using tadah::core::stats::IDX_P90;
  return BoundSuggestion(
      pp.r_percentiles[IDX_P10], pp.r_percentiles[IDX_P90],
      "Gaussian exp(-SGRIDMB (r - CGRIDMB)^2) centred on the first "
      "neighbour shell: the 10th to 90th percentiles of the neighbour "
      "distances (no peak found)");
}

BoundSuggestion gaussian_width_suggestion(const DataStats &d,
                                          const SeedView &s, int arg_pos) {
  if (d.per_pair.empty()) return std::make_pair(0.0, 0.0);
  const auto &pp = d.per_pair[0];
  double c = 0.0;
  if (seed_centre(s, arg_pos, c) && c == 0.0) {
    double r1 = 0.0, r2 = 0.0;
    std::string between;
    shells(d, r1, r2, between);
    const double gap = r2 * r2 - r1 * r1;
    if (!(gap > 0.0))
      return BoundSuggestion::no_row(
          "Gaussian exp(-SGRIDMB r^2) centred at the origin: the neighbour "
          "distances show no second shell beyond the first");
    const double lo = std::log(1.0 / 0.9) / gap;
    const double drop = std::log(20.0) / gap;
    const double reach = gaussian_reach_width(pp.r_min);
    const double hi = std::min(drop, reach);
    std::string why =
        "Gaussian exp(-SGRIDMB r^2) centred at the origin (the seed's "
        "CGRIDMB 0): its value drops by a factor of 0.9 to 1/20 " + between;
    if (reach < drop)
      why += ", up to " + num(reach) + ", where it is 1e-6 at r_min (" +
             num(pp.r_min) + " Å)";
    if (!(lo < hi))
      return BoundSuggestion::no_row(
          why + ": the widths it allows do not reach r_min");
    return BoundSuggestion(lo, hi, why);
  }
  double mu = 0.0, wn = 0.0;
  std::string width;
  if (seed_centre(s, arg_pos, c) && nearest_shell(d, c, mu, wn, width))
    return BoundSuggestion(
        0.25 / (wn * wn), 9.0 / (wn * wn),
        "Gaussian exp(-SGRIDMB (r - CGRIDMB)^2) centred on the seed's centre "
        "(" + num(c) + " Å): 0.25/w^2 to 9/w^2 with w " + width);
  using tadah::core::stats::IDX_P10;
  using tadah::core::stats::IDX_P90;
  double w = pp.r_percentiles[IDX_P90] - pp.r_percentiles[IDX_P10];
  std::string what = "the spread from the 10th to the 90th percentile of "
                     "the neighbour distances";
  if (!pp.r_peaks.empty() && pp.r_peaks[0].fwhm > 0.0) {
    w = pp.r_peaks[0].fwhm;
    what = "the FWHM of the first neighbour shell";
  }
  if (!(w > 0.0))
    return BoundSuggestion::fixed_range(
        0.1, 1.0,
        "fixed range: the neighbour distances have no spread to take a "
        "width of the Gaussian from");
  return BoundSuggestion(
      0.25 / (w * w), 9.0 / (w * w),
      "Gaussian exp(-SGRIDMB (r - CGRIDMB)^2) centred on the first "
      "neighbour shell: 0.25/w^2 to 9/w^2 with w = " + num(w) + " Å, " +
          what);
}

}  // namespace models
}  // namespace tadah
