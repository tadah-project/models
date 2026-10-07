#include <tadah/models/functions/weights_fit.h>

#include <tadah/core/intercept_layout.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace tadah {
namespace models {

namespace {

std::size_t stored_size(const tadah::core::Context &c, const char *key) {
  return c.has_stored(key) ? c.size(key) : 0;
}

std::string atom_list(const tadah::core::Context &c) {
  std::string s;
  const std::size_t n = stored_size(c, "ATOM");
  for (std::size_t i = 0; i < n; ++i)
    s += (i ? " " : "") + c.get<std::string>("ATOM", i);
  return s;
}

}  // namespace

TrainedSizes trained_sizes(const tadah::core::Context &c) {
  TrainedSizes s;
  s.weights = stored_size(c, "WEIGHTS");
  s.nmean = stored_size(c, "NMEAN");
  s.nstdev = stored_size(c, "NSTDEV");
  return s;
}

void require_weights_fit(const tadah::core::Context &c, std::size_t dsize,
                         const Function_Base &f, const TrainedSizes &s) {
  const bool norm = c.contains("NORM") && c.get<bool>("NORM");
  const std::string D = std::to_string(dsize);

  std::vector<std::string> problems;
  const std::size_t n_w = f.weights_count(dsize);
  if (s.weights != n_w)
    problems.push_back("WEIGHTS has " + std::to_string(s.weights) +
                       " values, but " + f.get_label() + " needs " +
                       std::to_string(n_w) + " for a descriptor of DSIZE " +
                       D);
  const std::size_t bvs = f.basis_vector_size();
  if (bvs != 0 && bvs != dsize)
    problems.push_back("the BASIS vectors (SBASIS " + std::to_string(n_w) +
                       ") have " + std::to_string(bvs) +
                       " entries each, but the descriptor has DSIZE " + D);
  if (norm) {
    if (s.nmean != dsize)
      problems.push_back("NORM true: NMEAN has " + std::to_string(s.nmean) +
                         " values, but the descriptor has DSIZE " + D +
                         " entries");
    if (s.nstdev != dsize)
      problems.push_back("NORM true: NSTDEV has " +
                         std::to_string(s.nstdev) +
                         " values, but the descriptor has DSIZE " + D +
                         " entries");
  }
  if (problems.empty()) return;

  const std::size_t n_int = tadah::core::n_intercept_columns(c);
  const std::string atoms = atom_list(c);
  std::string layout;
  if (n_int > 0)
    layout = "DSIZE counts " + std::to_string(n_int) + " intercept column" +
             (n_int > 1 ? "s" : "") +
             ": BIAS true, one per ATOM species, ATOM " + atoms;
  else
    layout = "BIAS false: no intercept column" +
             (atoms.empty() ? std::string() : "; ATOM " + atoms);

  std::string msg = "The potential does not fit its descriptor: ";
  for (std::size_t i = 0; i < problems.size(); ++i)
    msg += (i ? "; " : "") + problems[i];
  msg += " (" + layout + ").";

  /* The sizes of the same descriptor with one intercept column instead of
   * N: a multi-species BIAS potential trained before 1.6.0. */
  if (n_int >= 2 && dsize >= n_int) {
    const std::size_t one = dsize - (n_int - 1);
    bool legacy = bvs == 0 ? s.weights == f.weights_count(one) : bvs == one;
    if (norm && (s.nmean == one || s.nstdev == one)) legacy = true;
    if (legacy)
      msg += " Its sizes are those of one intercept column shared by all "
             "species, the layout in which tadah 1.5.x and earlier (and "
             "1.6.0 development builds before the per-species intercept "
             "fix) trained BIAS true with several species; pair_style "
             "tadah used one intercept column per species with it, so its "
             "LAMMPS energies and forces were wrong. Its weights cannot be "
             "used: retrain it (tadah train -c <its training configuration, "
             "or this potential>).";
  }
  throw std::runtime_error(msg);
}

void require_weights_fit(const tadah::core::Context &c, std::size_t dsize,
                         const Function_Base &f) {
  require_weights_fit(c, dsize, f, trained_sizes(c));
}

} // namespace models
} // namespace tadah
