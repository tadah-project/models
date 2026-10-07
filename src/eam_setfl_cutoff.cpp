#include <tadah/models/descriptors/eam_setfl_cutoff.h>

#include <tadah/core/number_format.h>

#include <cstddef>
#include <stdexcept>

namespace tadah {
namespace models {

void require_setfl_cutoff(const tadah::core::Context &c,
                          const std::string &rcut_key,
                          const std::string &label, const EAMSetflIO &eam) {
  using tadah::core::format_roundtrip;
  const double rc = eam.rcut();
  for (std::size_t i = 0; i < c.size(rcut_key); ++i) {
    const double v = c.get<double>(rcut_key, i);
    if (v == rc) continue;
    throw std::runtime_error(
        label + ": " + rcut_key + " " + format_roundtrip(v) +
        " differs from the cutoff " + format_roundtrip(rc) +
        " of the setfl file '" + eam.path() + "' (SETFL); set " + rcut_key +
        " to the file's cutoff, " + format_roundtrip(rc) + ".");
  }
}

}  // namespace models
}  // namespace tadah
