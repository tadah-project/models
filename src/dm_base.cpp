#include <tadah/models/descriptors/dm/dm_base.h>
namespace tadah {
namespace models {

void DM_Base::set_rfidx(size_t rfidx_) { rfidx=rfidx_; }
size_t DM_Base::get_rfidx() { return rfidx; }

std::vector<std::string> DM_Base::get_init_atoms(tadah::core::Config &c) {
  if (!c.exist("TYPEMB")) {
    std::cerr << "Warning: TYPEMB key was not found. "
      "Atom pairs must be initialized manually." << std::endl;
    return {};
  }
  std::vector<std::string> init_atoms(c.size("TYPEMB"));
  c.get("TYPEMB", init_atoms);
  init_atoms.erase(init_atoms.begin(), init_atoms.begin() + nparams + 1);
  return init_atoms;
}

} // namespace models
} // namespace tadah

