#include <tadah/models/descriptors/d_base.h>
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
D_Base::~D_Base() {
  if (this->manage_memory && fcut) {
    delete fcut;
  }
}

D_Base::D_Base() {}
D_Base::D_Base(Config &c):
  verbose(c.get<int>("VERBOSE"))
{
  for (size_t i=0; i<c.size("ATOMS"); ++i) {
    std::string symbol = c.get<std::string>("ATOMS",i);
    double wi = c.get<double>("WATOMS",i);
    int Z = PeriodicTable::find_by_symbol(symbol).Z;
    weights[Z]=wi;
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
