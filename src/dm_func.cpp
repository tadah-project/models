#include <tadah/models/descriptors/dm/dm_func.h>
#include <cmath>

namespace tadah {
namespace models {

F_RLR::F_RLR(const tadah::core::Context &c, size_t cSize, size_t sSize):
  context(&c)
{
  if (!context->get<bool>("INITMB")) return;
  D_Base::get_grid(*context,"SEMBFUNC",sfuncgrid);
  D_Base::get_grid(*context,"CEMBFUNC",cfuncgrid);
  // check if size is correct, throw otherwise
  if (sfuncgrid.size() != sSize) {
    throw std::runtime_error("F_RLR: size of SEMBFUNC does not match declared size.");
  }
  if (cfuncgrid.size() != cSize) {
    throw std::runtime_error("F_RLR: size of CEMBFUNC does not match declared size.");
  }
}
double F_RLR::calc_F(double rho,size_t p) const {
  return sfuncgrid[p]*std::abs(rho)*std::log(cfuncgrid[p]*std::abs(rho) + eps);
}
double F_RLR::calc_dF(double rho,size_t p) const {
  return sfuncgrid[p]*(std::log(cfuncgrid[p]*std::abs(rho) + eps)+1);
}

//////////////////////////
// F_SQRT
//////////////////////////
F_SQRT::F_SQRT(const tadah::core::Context &c, size_t, size_t sSize):
  context(&c)
{
  if (!context->get<bool>("INITMB")) return;
  if (context->contains("SEMBFUNC") && context->size("SEMBFUNC") > 0) {
    D_Base::get_grid(*context, "SEMBFUNC", sfuncgrid);
    // check if size is correct, throw otherwise
    if (sfuncgrid.size() != sSize) {
      throw std::runtime_error("F_SQRT: size of SEMBFUNC does not match declared size.");
    }
  } else {
    sfuncgrid.resize(sSize);
    for (size_t i=0; i<sSize; ++i) {
      sfuncgrid[i] = 1.0;
    }
  }
}
double F_SQRT::calc_F(double rho,size_t p) const {
    return sfuncgrid[p]*std::sqrt(std::abs(rho) + eps);
}
double F_SQRT::calc_dF(double rho,size_t p) const {
    return 0.5*sfuncgrid[p]/(std::sqrt(std::abs(rho) + eps));
}

} // namespace models
} // namespace tadah
