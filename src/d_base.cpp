#include "tadah/core/core_types.h"
#include "tadah/core/utils.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <tadah/models/descriptors/d_base.h>
#include <utility>
#include <vector>
namespace tadah {
namespace models {

core::v_type D_Base::get_grid(const std::vector<std::string>&vold) {
  core::v_type vnew;
  auto it = vold.begin();
  while (it != vold.end()) {
    std::string token = *it++;
    if (core::to_upper(token) == "LIN") {
      size_t cg = std::stoul(*it++);
      double cgmin = std::stod(*it++);
      double cgmax = std::stod(*it++);
      core::v_type v = core::linspace(cgmin, cgmax, cg);
      vnew.insert(vnew.end(), v.begin(), v.end());
    }
    else if (core::to_upper(token) == "GEOM") {
      size_t cg = std::stoul(*it++);
      double cgmin = std::stod(*it++);
      double cgmax = std::stod(*it++);
      core::v_type v = core::geomspace(cgmin, cgmax, cg, exp(1));
      vnew.insert(vnew.end(), v.begin(), v.end());
    }
    else if (!core::is_number(token)) {
      throw std::invalid_argument("Invalid grid specification: " + token);
    }
    else {
      vnew.push_back(std::stod(token));
    }
  }
  return vnew;
}

void D_Base::get_grid(const tadah::core::Context &c, const std::string &key, tadah::core::v_type &v) {
  std::vector<std::string> grid_line(c.size(key));
  c.get<>(key, grid_line);
  v = get_grid(grid_line);
}

D_Base::~D_Base() {
  if (this->manage_memory && fcut) {
    delete fcut;
  }
}

D_Base::D_Base() {
  for (int i = 0; i < 119; ++i) std::fill_n(weights[i], 119, 1.0);
}
D_Base::D_Base(const tadah::core::Context &c):
  verbose(c.get<int>("VERBOSE"))
{
  for (int i = 0; i < 119; ++i) std::fill_n(weights[i], 119, 1.0);
  update_weights(c);
}
void D_Base::update_weights(const tadah::core::Context &c) {
  watom_Z_by_arg_pos.clear();
  watom_pair_by_arg_pos.clear();
  watom_mode = WatomMode::Empty;
  if (!c.contains("WATOM") || c.size("WATOM") == 0) return;

  const size_t n = c.size("WATOM");
  const bool per_species =
    tadah::core::is_number(c.get<std::string>("WATOM", 0));

  if (per_species) {
    if (!c.contains("ATOM")) return;
    watom_mode = WatomMode::PerSpecies;
    std::vector<int> Zs(n);
    for (size_t i = 0; i < n; ++i) {
      Zs[i] = tadah::core::PeriodicTable::find_by_symbol(
                c.get<std::string>("ATOM", i)).Z;
    }
    watom_Z_by_arg_pos = Zs;
    for (size_t i = 0; i < n; ++i) {
      const double wi = c.get<double>("WATOM", i);
      weights[Zs[i]][Zs[i]] = wi;
      for (size_t j = i + 1; j < n; ++j) {
        const double w = std::sqrt(wi * c.get<double>("WATOM", j));
        weights[Zs[i]][Zs[j]] = w;
        weights[Zs[j]][Zs[i]] = w;
      }
    }
  } else {
    if (n % 3 != 0) {
      throw std::runtime_error(
        "WATOM pair format: expected sym sym value triples");
    }
    watom_mode = WatomMode::PairTriple;
    watom_pair_by_arg_pos.reserve(n / 3);
    for (size_t i = 0; i < n; i += 3) {
      const int Zi = tadah::core::PeriodicTable::find_by_symbol(
                       c.get<std::string>("WATOM", i)).Z;
      const int Zj = tadah::core::PeriodicTable::find_by_symbol(
                       c.get<std::string>("WATOM", i + 1)).Z;
      const double w = c.get<double>("WATOM", i + 2);
      weights[Zi][Zj] = w;
      weights[Zj][Zi] = w;
      watom_pair_by_arg_pos.emplace_back(Zi, Zj);
    }
  }
}
int D_Base::get_atom_Z_for_arg_pos(int arg_pos) const {
  if (arg_pos < 1) return 0;
  const size_t idx = static_cast<size_t>(arg_pos - 1);
  if (idx >= watom_Z_by_arg_pos.size()) return 0;
  return watom_Z_by_arg_pos[idx];
}
bool D_Base::has_watom_arg_pos(int arg_pos) const {
  if (arg_pos < 1) return false;
  const size_t idx = static_cast<size_t>(arg_pos - 1);
  switch (watom_mode) {
    case WatomMode::PerSpecies: return idx < watom_Z_by_arg_pos.size();
    case WatomMode::PairTriple: return idx < watom_pair_by_arg_pos.size();
    case WatomMode::Empty:      return false;
  }
  return false;
}
double D_Base::dweight_dw_n(int Zi, int Zj, int arg_pos) const {
  if (arg_pos < 1) return 0.0;
  const size_t idx = static_cast<size_t>(arg_pos - 1);
  if (watom_mode == WatomMode::PerSpecies) {
    if (idx >= watom_Z_by_arg_pos.size()) return 0.0;
    const int Z_target = watom_Z_by_arg_pos[idx];
    const bool zi_eq = (Zi == Z_target);
    const bool zj_eq = (Zj == Z_target);
    if (!zi_eq && !zj_eq) return 0.0;
    if (zi_eq && zj_eq)   return 1.0;
    const double w_n     = weights[Z_target][Z_target];
    const int    other   = zi_eq ? Zj : Zi;
    const double w_other = weights[other][other];
    if (w_n <= 0.0) return 0.0;
    return 0.5 * std::sqrt(w_other / w_n);
  }
  if (watom_mode == WatomMode::PairTriple) {
    if (idx >= watom_pair_by_arg_pos.size()) return 0.0;
    const auto &pair = watom_pair_by_arg_pos[idx];
    const int Za = pair.first;
    const int Zb = pair.second;
    return ((Zi == Za && Zj == Zb) || (Zi == Zb && Zj == Za)) ? 1.0 : 0.0;
  }
  return 0.0;
}
void D_Base::printGrid(core::v_type &v, const std::string &key, core::LogLevel logLevel) {
  core::Logger::instance().log(logLevel, label()," ",key,":\n  ", v);
}

void D_Base::set_fidx(size_t fidx_) { fidx=fidx_; }
size_t D_Base::get_fidx() const { return fidx; }
void D_Base::set_fcut(Cut_Base* cut, bool manage_memory) {
  if (this->manage_memory && fcut != cut) {
    delete fcut;
  }
  fcut = cut;
  this->manage_memory = manage_memory;
}
int D_Base::get_arg_pos(const std::string &key) const {
  auto it = std::find(keys.begin(), keys.end(), key);

  if (it == keys.end()) return -1;  // incl. data_keys: not OPTIM-addressable

  // fparam = leading family params before the per-key size tokens
  const long fparam = static_cast<long>(nparams) -
                      static_cast<long>(keys.size()) -
                      static_cast<long>(data_keys.size());
  return static_cast<int>(std::distance(keys.begin(), it) + fparam);
}
size_t D_Base::size() const { return s; }
double D_Base::get_rcut() const { 
  return fcut->get_rcut(); 
}
bool D_Base::is_init_for_atoms(int Zi, int Zj) const {
  return init_for_atoms_map.is_init(Zi,Zj) ;
}
void D_Base::init_for_atoms(int Zi, int Zj) {
  init_for_atoms_map.init(Zi,Zj);
}
void D_Base::uninit_for_atoms(int Zi, int Zj) {
  init_for_atoms_map.uninit(Zi,Zj);
}
void D_Base::init_for_atoms(const std::vector<std::string> &Zs) {
  if (Zs.size() % 2 != 0) {
    throw std::invalid_argument("The vector size must be even.");
  }
  
  for (size_t i = 0; i < Zs.size(); i += 2) {
    if (Zs[i] == "*" && Zs[i+1] != "*") {
      int Zj = tadah::core::PeriodicTable::find_by_symbol(Zs[i+1]).Z;
      for (int Zi = 1; Zi < 119; ++Zi) {
        init_for_atoms(Zi, Zj);
      }
    } 
    else if (Zs[i] != "*" && Zs[i+1] == "*") {
      int Zi = tadah::core::PeriodicTable::find_by_symbol(Zs[i]).Z;
      for (int Zj = 1; Zj < 119; ++Zj) {
        init_for_atoms(Zi, Zj);
      }
    } 
    else if (Zs[i] == "*" && Zs[i+1] == "*") {
      for (int Zi = 1; Zi < 119; ++Zi) {
        for (int Zj = Zi; Zj < 119; ++Zj) {
          init_for_atoms(Zi, Zj);
        }
      }
    } 
    else {
      int Zi = tadah::core::PeriodicTable::find_by_symbol(Zs[i]).Z;
      int Zj = tadah::core::PeriodicTable::find_by_symbol(Zs[i+1]).Z;
      init_for_atoms(Zi, Zj);
    }
  }
}
std::vector<std::string> D_Base::element_tokens(const std::string &key,
                                                std::vector<std::string> line) const {
  const size_t nsizes = keys.size() + data_keys.size();
  const std::string name =
      line.empty() ? std::string("<no descriptor>") : line[0];
  if (nsizes > nparams) {
    throw std::logic_error(name + " declares " + std::to_string(nsizes) +
                           " list keys but takes " + std::to_string(nparams) +
                           " parameters");
  }
  std::string text = key;
  for (const auto &t : line) text += " " + t;
  const std::string takes =
      key + " " + name + " takes " + std::to_string(nparams) +
      (nparams == 1 ? " parameter" : " parameters") +
      " before its element pairs; ";

  const size_t have = line.empty() ? 0 : line.size() - 1;
  if (have < nparams) {
    throw std::runtime_error(takes + "the line has " + std::to_string(have) +
                             ": " + text);
  }

  // SIZE tokens: line[1 + fparam + k] is the size of the k-th list key.
  const size_t fparam = nparams - nsizes;
  for (size_t k = 0; k < nsizes; ++k) {
    const std::string &token = line[1 + fparam + k];
    if (tadah::core::is_integer(token)) continue;
    const std::string &list_key =
        k < keys.size() ? keys[k] : data_keys[k - keys.size()];
    throw std::runtime_error(takes + "parameter " +
                             std::to_string(fparam + k + 1) +
                             ", the number of " + list_key +
                             " values, must be a non-negative integer, not \"" +
                             token + "\": " + text);
  }

  line.erase(line.begin(),
             line.begin() + static_cast<std::ptrdiff_t>(nparams + 1));
  return line;
}

} // namespace models
} // namespace tadah

