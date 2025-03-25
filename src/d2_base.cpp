#include <tadah/models/descriptors/d2/d2_base.h>
namespace tadah {
namespace models {


std::vector<std::string> D2_Base::get_init_atoms(const tadah::core::Context &c) {
  if (!c.contains("TYPE2B")) {
    std::cerr << "Warning: TYPE2B key was not found. "
      "Atom pairs must be initialized manually." << std::endl;
    return {};
  }
  std::vector<std::string> init_atoms(c.size("TYPE2B"));
  c.get("TYPE2B", init_atoms);
  init_atoms.erase(init_atoms.begin(), init_atoms.begin() + nparams+1);
  return init_atoms;
}

} // namespace models
} // namespace tadah

