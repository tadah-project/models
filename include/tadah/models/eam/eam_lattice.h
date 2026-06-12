#ifndef TADAH_MODELS_EAM_LATTICE_H
#define TADAH_MODELS_EAM_LATTICE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

namespace tadah {
namespace models {
namespace eam_lattice {

/** \brief Ideal-lattice neighbour distances for monatomic bcc/fcc/hcp crystals.
 *
 * Shared by the refit reachability analysis (which embedding densities / pair
 * distances a compressed crystal can actually visit) and the refit crystal
 * energy self-test. All atoms in these lattices are symmetry-equivalent, so
 * distances from a single basis atom describe every atom. hcp uses the ideal
 * c/a = sqrt(8/3); for energy *differences* between two potentials evaluated
 * on the same geometry this is exact, and for reachability it is a perfectly
 * adequate probe geometry.
 */

/// true if \p lattice is one of the supported monatomic lattices
inline bool supported(const std::string &lattice) {
  return lattice == "bcc" || lattice == "fcc" || lattice == "hcp";
}

/// nearest-neighbour distance for lattice constant \p a
inline double nn_distance(const std::string &lattice, double a) {
  if (lattice == "fcc") return a / std::sqrt(2.0);
  if (lattice == "hcp") return a;
  return a * std::sqrt(3.0) / 2.0;  // bcc (and fallback)
}

/** All neighbour distances (one central atom) within \p rcut for lattice
 *  constant \p a. Unknown lattice falls back to bcc (caller should report). */
inline std::vector<double> neighbour_distances(const std::string &lattice,
                                               double a, double rcut) {
  // lattice vectors (rows) and fractional basis
  double A[3][3];
  std::vector<std::array<double, 3>> basis;
  if (lattice == "fcc") {
    A[0][0] = a; A[0][1] = 0; A[0][2] = 0;
    A[1][0] = 0; A[1][1] = a; A[1][2] = 0;
    A[2][0] = 0; A[2][1] = 0; A[2][2] = a;
    basis = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  } else if (lattice == "hcp") {
    const double c = a * std::sqrt(8.0 / 3.0);
    A[0][0] = a;       A[0][1] = 0;                  A[0][2] = 0;
    A[1][0] = -0.5 * a; A[1][1] = a * std::sqrt(3.0) / 2.0; A[1][2] = 0;
    A[2][0] = 0;       A[2][1] = 0;                  A[2][2] = c;
    basis = {{0, 0, 0}, {1.0 / 3.0, 2.0 / 3.0, 0.5}};
  } else {  // bcc + fallback
    A[0][0] = a; A[0][1] = 0; A[0][2] = 0;
    A[1][0] = 0; A[1][1] = a; A[1][2] = 0;
    A[2][0] = 0; A[2][1] = 0; A[2][2] = a;
    basis = {{0, 0, 0}, {0.5, 0.5, 0.5}};
  }
  const double amin = std::min({std::fabs(A[0][0]), std::fabs(A[1][1]),
                                std::fabs(A[2][2])});
  const int N = static_cast<int>(std::ceil(rcut / amin)) + 1;
  std::vector<double> out;
  for (int i = -N; i <= N; ++i)
    for (int j = -N; j <= N; ++j)
      for (int k = -N; k <= N; ++k)
        for (const auto &b : basis) {
          const double fx = i + b[0], fy = j + b[1], fz = k + b[2];
          const double x = fx * A[0][0] + fy * A[1][0] + fz * A[2][0];
          const double y = fx * A[0][1] + fy * A[1][1] + fz * A[2][1];
          const double z = fx * A[0][2] + fy * A[1][2] + fz * A[2][2];
          const double d = std::sqrt(x * x + y * y + z * z);
          if (d > 1e-9 && d <= rcut) out.push_back(d);
        }
  return out;
}

/** Energy per atom of a monatomic bcc/fcc/hcp crystal at lattice constant
 *  \p a by a direct pair/density sum within \p rcut. NaN for an unsupported
 *  lattice. `phi`,`rho`,`F` may be tabulated originals or fitted functions —
 *  shared by the refit self-test and the compact-mode energy metric. */
inline double crystal_epa(const std::string &lat, double a, double rcut,
                          const std::function<double(double)> &phi,
                          const std::function<double(double)> &rho,
                          const std::function<double(double)> &F) {
  if (!supported(lat)) return std::nan("");
  const auto dists = neighbour_distances(lat, a, rcut);
  double half_pair = 0.0, dens = 0.0;  // sums seen by one (any) atom
  for (double r : dists) {
    half_pair += 0.5 * phi(r);
    dens += rho(r);
  }
  return half_pair + F(dens);  // monatomic: every atom is equivalent
}

}  // namespace eam_lattice
}  // namespace models
}  // namespace tadah
#endif  // TADAH_MODELS_EAM_LATTICE_H
