#include <tadah/models/descriptors/d3/d3_base.h>
namespace tadah {
namespace models {

//template struct Registry<D3_Base, tadah::core::Config &>;
std::vector<std::string> D3_Base::get_init_atoms(tadah::core::Config &) {
  return {};
}


} // namespace models
} // namespace tadah

