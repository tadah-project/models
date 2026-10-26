#ifndef TADAH_MODELS_DESCRIPTORS_DENSITY_REACH_H
#define TADAH_MODELS_DESCRIPTORS_DENSITY_REACH_H

// The reach of a radial density channel, and the OPTIM ranges of its centre
// and width that the seed and the data suggest.
//
// A many-body density sums a radial function g(r) of each channel over the
// neighbours within its cutoff, so a channel is judged by what g gives at the
// neighbour distances of the data, not by where its centre lies: a Gaussian
// exp(-SGRIDMB r^2) centred at the origin reaches the first neighbour shell
// when its width SGRIDMB is small enough, and a blip centred below the
// nearest neighbour reaches it when its support does.
//
// One constant, kDensityReachFloor, is the value below which a radial
// function or a density counts as reaching nothing: the CGRIDMB checks of the
// densities (radial_reach_check), the rho-collapse checks of DM_mBlipEnv and
// DM_Blip, and the widths tadah analysis suggest_bounds proposes for a
// Gaussian at the origin (gaussian_width_suggestion) use it.

#include <string>

#include <tadah/core/stats/data_stats.h>
#include <tadah/models/descriptors/basis_kind.h>
#include <tadah/models/descriptors/param_meta.h>

namespace tadah {
namespace models {

/** \brief The value below which a radial function at a neighbour distance,
 *  or the density of an atom, counts as zero: 1e-6. */
inline constexpr double kDensityReachFloor = 1.0e-6;

/** \brief The neighbour distances a density with cutoff @p rcut sums over:
 *  from the smallest r_min of the element pairs of @p d to the smaller of
 *  their largest r_max and @p rcut (@p rcut <= 0: no cutoff given, r_max).
 *  `empty` when the data have no pair or none closer than @p rcut. */
struct NeighbourRange {
  double lo = 0.0;
  double hi = 0.0;
  bool empty = true;
};
NeighbourRange neighbour_range(const tadah::core::stats::DataStats &d,
                               double rcut);

/** \brief The width eta above which a Gaussian exp(-eta (r - c)^2) is
 *  below kDensityReachFloor at the distance @p dist from its centre:
 *  ln(1 / kDensityReachFloor) / dist^2 (infinite for @p dist 0). */
double gaussian_reach_width(double dist);

/** \brief The reach check of the radial channels of a many-body density:
 *  each centre CGRIDMB(i) with its width SGRIDMB(i) and the cutoff RCUTMB
 *  of @p s, against the neighbour distances of @p d (neighbour_range).
 *
 *  A Gaussian channel (@p kind Gaussian, exp(-SGRIDMB (r - CGRIDMB)^2))
 *  fails when it is below kDensityReachFloor at every distance of the
 *  range; a blip channel (Blip, b(SGRIDMB (r - CGRIDMB)), support |r -
 *  CGRIDMB| < 2/SGRIDMB) when its support holds none. A width of 0 or
 *  less is no reach failure (other rules refuse it). HealthFail names the
 *  first channel that fails by its entry of the configuration
 *  (SeedView::entry_number(): CGRIDMB(3) for the second channel of a
 *  joined line's descriptor after one of another), its centre and width,
 *  the form and the range;
 *  also a range that is empty (RCUTMB below the smallest neighbour
 *  distance, or no pair within it, for data of one atom or more) and a
 *  centre without a width. HealthOk without CGRIDMB, or for data without
 *  atoms. @throws std::logic_error for a @p kind other than Gaussian and
 *  Blip. */
HealthVerdict radial_reach_check(const tadah::core::stats::DataStats &d,
                                 const SeedView &s, BasisKind kind);

/** \brief The OPTIM range of the centre CGRIDMB(@p arg_pos) (-1: the
 *  first) of a Gaussian density channel, by the seed's centre c:
 *  - at the origin (CGRIDMB 0): no row, the width alone sets the reach;
 *  - elsewhere: c +- the FWHM of the neighbour shell nearest to c (the
 *    peak of the neighbour distances nearest to it; 0.2 times its position
 *    for a peak without a FWHM), its low end not below r_min and its high
 *    end not above the neighbour distances within RCUTMB
 *    (neighbour_range()) unless c is (then the range ends at c): the range
 *    holds c;
 *  - without a seed centre: the first peak of the neighbour distances +-
 *    its FWHM, from r_min;
 *  - without peaks: the 10th to 90th percentiles of the distances. */
BoundSuggestion gaussian_centre_suggestion(
    const tadah::core::stats::DataStats &d, const SeedView &s, int arg_pos);

/** \brief The OPTIM range of the width SGRIDMB(@p arg_pos) (-1: the first)
 *  of a Gaussian density channel, by the seed's centre CGRIDMB c of the
 *  same entry:
 *  - at the origin: the widths at which exp(-SGRIDMB r^2) drops by a
 *    factor of 0.9 to 1/20 from the first neighbour shell r1 to the second
 *    r2, ln(1/0.9) / (r2^2 - r1^2) to ln(20) / (r2^2 - r1^2), the upper end
 *    at most gaussian_reach_width(r_min) (the shells: the first two peaks
 *    of the neighbour distances; without them, their 10th percentile and
 *    median);
 *  - elsewhere: 0.25 / w^2 to 9 / w^2 with w the FWHM of the neighbour
 *    shell nearest to c (as gaussian_centre_suggestion() takes it);
 *  - without a seed centre, or without peaks: 0.25 / w^2 to 9 / w^2 with
 *    w the FWHM of the first peak (without peaks, or without a FWHM, the
 *    spread from the 10th to the 90th percentile), a fixed range 0.1 to 1
 *    when the distances have no spread. */
BoundSuggestion gaussian_width_suggestion(
    const tadah::core::stats::DataStats &d, const SeedView &s, int arg_pos);

}  // namespace models
}  // namespace tadah

#endif  // TADAH_MODELS_DESCRIPTORS_DENSITY_REACH_H
