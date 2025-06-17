#include "tadah/core/core_types.h"
#include "tadah/core/utils.h"
#include <string>
#include <tadah/models/descriptors/d_base.h>
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

D_Base::D_Base() {}
D_Base::D_Base(const tadah::core::Context &c):
  verbose(c.get<int>("VERBOSE"))
{
  std::fill_n(weights, 119, 1.0);
  if (c.contains("ATOM")) {
    for (size_t i=0; i<c.size("ATOM"); ++i) {
      std::string symbol = c.get<std::string>("ATOM",i);
      double wi = c.get<double>("WATOM",i);
      int Z = tadah::core::PeriodicTable::find_by_symbol(symbol).Z;
      weights[Z]=wi;
    }
  }
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

  if (it == keys.end()) return -1;

  return std::distance(keys.begin(), it)+(nparams-keys.size());
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
std::vector<std::string> D_Base::get_init_atoms(const tadah::core::Context &c, std::string type) {
  std::vector<std::string> init_atoms(c.size(type));
  c.get(type, init_atoms);
  std::unique_ptr<D_Base> d(tadah::core::factory<D_Base>(init_atoms[0]));
  size_t nparams = d->nparams;
  init_atoms.erase(init_atoms.begin(), init_atoms.begin() + nparams);
  return init_atoms;
}

} // namespace models
} // namespace tadah

