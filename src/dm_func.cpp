#include <tadah/models/descriptors/dm/dm_func.h>

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
double F_RLR::calc_F(double rho,size_t p) {
  if (rho<=0) return 0;
  return sfuncgrid[p]*rho*log(cfuncgrid[p]*rho);
}
double F_RLR::calc_dF(double rho,size_t p) {
  if (rho<=0) return 0;
  return sfuncgrid[p]*(log(cfuncgrid[p]*rho)+1);
}

//////////////////////////
// F_SQRT
//////////////////////////
F_SQRT::F_SQRT(const tadah::core::Context &c, size_t, size_t sSize):
  context(&c)
{
  if (!context->get<bool>("INITMB")) return;
  if (context->contains("SEMBFUNC")) {
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
double F_SQRT::calc_F(double rho,size_t p) {
  if (rho>0)
    return sfuncgrid[p]*std::sqrt(rho);
  return 0;
}
double F_SQRT::calc_dF(double rho,size_t p) {
  if (rho>0)
    return 0.5*sfuncgrid[p]/(std::sqrt(rho));
  return 0;
}

} // namespace models
} // namespace tadah
