#include "tadah/core/core_types.h"
#include <string>
#include <tadah/models/descriptors/d_base.h>
#include <vector>
void D_Base::get_grid(const Config &c, std::string key, v_type &v) {
  if (c.get<double>(key) < 0) {
    // generate grid automatically if first element
    // is<zero and int. If it is double than assume that
    // grid is provided manually and min value is just smaller
    // than zero.
    double int_part;
    if (std::modf ( c.get<double>(key), &int_part ) == 0.0) {
      // automatic generation
      size_t cg = abs(c.get<int>(key,1));
      double cgmin = c.get<double>(key,2);
      double cgmax = c.get<double>(key,3);
      if (c.get<int>(key) == -1) {
        v = linspace(cgmin, cgmax, cg);
      }
      else if (c.get<int>(key) == -2) {
        v = logspace(cgmin, cgmax, cg, exp(1));
      }
      else {
        throw std::runtime_error(key+" algorithm not supported\n");
      }
    }
    else {
      v.resize(c.size(key));
      c.get<v_type>(key,v);
    }
  }
  else {
    v.resize(c.size(key));
    c.get<v_type>(key,v);
  }
}
v_type D_Base::get_grid(std::vector<std::string>&vold) {
  v_type vnew;
  auto it = vold.begin();
  while (it != vold.end()) {
    std::string token = *it++;
    if (token == "-1") {
      size_t cg = std::stoul(*it++);
      double cgmin = std::stod(*it++);
      double cgmax = std::stod(*it++);
      v_type v = linspace(cgmin, cgmax, cg);
      vnew.insert(vnew.end(), v.begin(), v.end());
    }
    else if (token == "-2") {
      size_t cg = std::stoul(*it++);
      double cgmin = std::stod(*it++);
      double cgmax = std::stod(*it++);
      v_type v = logspace(cgmin, cgmax, cg, exp(1));
      vnew.insert(vnew.end(), v.begin(), v.end());
    }
    else if (std::stod(token) < -2) {
        throw std::runtime_error("Either algorithm is not supported or negative value was encountered in a grid: " + token);
    }
    else {
      vnew.push_back(std::stod(token));
    }
  }
  return vnew;
}
D_Base::~D_Base() {
  if (this->manage_memory && fcut) {
    delete fcut;
  }
}

D_Base::D_Base() {}
D_Base::D_Base(Config &c):
  verbose(c.get<int>("VERBOSE"))
{
  std::fill_n(weights, 119, 1.0);
  if (c.exist("ATOMS")) {
    for (size_t i=0; i<c.size("ATOMS"); ++i) {
      std::string symbol = c.get<std::string>("ATOMS",i);
      double wi = c.get<double>("WATOMS",i);
      int Z = PeriodicTable::find_by_symbol(symbol).Z;
      weights[Z]=wi;
    }
  }
}
void D_Base::set_fidx(size_t fidx_) { fidx=fidx_; }
size_t D_Base::get_fidx() { return fidx; }
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
size_t D_Base::size() { return s; };
double D_Base::get_rcut() { 
  return fcut->get_rcut(); 
}
bool D_Base::is_init_for_atom(int Z) {
  return init_for_atoms_map[Z];
}
bool D_Base::is_init_for_atoms(int Zi, int Zj) {
  return init_for_atoms_map[Zi] && init_for_atoms_map[Zj];
}
void D_Base::init_for_atom(int Z) {
  init_for_atoms_map.set(Z);
}
void D_Base::init_for_atoms(const std::vector<std::string> &Zs) {
  for (const auto &Z_str : Zs){
  int Z = PeriodicTable::find_by_symbol(Z_str).Z;
  init_for_atom(Z);
  }
}
void D_Base::uninit_for_atom(int Z) {
  init_for_atoms_map.reset(Z);
}
std::vector<std::string> D_Base::get_init_atoms(Config &c, std::string type) {
  std::vector<std::string> init_atoms(c.size(type));
  c.get(type, init_atoms);
  std::unique_ptr<D_Base> d(CONFIG::factory<D_Base>(init_atoms[0]));
  size_t nparams = d->nparams;
  init_atoms.erase(init_atoms.begin(), init_atoms.begin() + nparams);
  return init_atoms;
}
