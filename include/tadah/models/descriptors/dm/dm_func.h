#ifndef DM_func_H
#define DM_func_H

#include <tadah/core/context.h>
#include <tadah/models/descriptors/d_base.h>

#include <cstdlib>

namespace tadah {
namespace models {

class F_Base {
public:
  /**
   * @brief Calculates the embedding function value.
   * @param rho Electron density.
   * @param p Index for model parameters.
   * @return Function value F(rho).
   */
  virtual double calc_F(double rho, size_t p) const = 0;

  /**
   * @brief Calculates the first derivative of the embedding function.
   * @param rho Electron density.
   * @param p Index for model parameters.
   * @return Derivative dF/d(rho).
   */
  virtual double calc_dF(double rho, size_t p) const = 0;
};

/**
 * @class F_RLR
 * @brief Implements an embedding function of the form: \f$ s \rho \log(c \rho) \f$
 *
 * This class supports embedding functions characterized by two main parameters:
 * - **SEMBFUNC**: Controls the depth, \f$ s \f$, of the embedding function. 
 * - **CEMBFUNC**: Determines the x-intercept, with the x-intercept at \f$ 1/c \f$.
 * 
 *  Require: size(**SEMBFUNC**)=size(**CEMBFUNC**)=size(**[C/S]GRIDMB**)
 * 
 * The number of keys for these parameters must match the entries in the mEAD descriptor.
 *
 */
class F_RLR:  public F_Base {
private:
  const tadah::core::Context *context=nullptr;
  tadah::core::v_type sfuncgrid;   // controls depth
  tadah::core::v_type cfuncgrid;   // control x-intercep at 1/c
  double eps=1e-10; // small value to avoid log(0)
public:
  F_RLR(const tadah::core::Context &context, size_t cSize, size_t sSize);
  F_RLR() = default;
  double calc_F(double rho,size_t p) const override;
  double calc_dF(double rho,size_t p) const override;
  //~F_RLR() {

  // Do not delete context!
  //}
};

/**
 * @class F_SQRT
 * @brief Implements \f$ s \sqrt{\rho} \f$.
 *
 * Optional parameter:
 * - **SEMBFUNC**: Controls the strength, \f$ s \f$, of the embedding function.
 *   If no; value is provided, the default is 1 for every s/cgrid point.
 *
 *  Require: size(**SEMBFUNC**)=(0 or size(**[C/S]GRIDMB**)) and size(**CEMBFUNC**)=0
 *
 */
class F_SQRT : public F_Base {
private:
  const tadah::core::Context *context=nullptr;
  tadah::core::v_type sfuncgrid;   // controls depth
  double eps=1e-10; // small value to avoid sqrt(0)

public:
  /**
   * @brief Constructs F_SQRT from a Context object.
   * @param context Reference to a Context instance.
   *
   * The constructor extracts or initializes relevant parameter for the model from
   * the provided configuration.
   *
   * SEMBFUNC valuas are set to 1 if not specified in the Context.
   */
  F_SQRT(const tadah::core::Context &context, size_t, size_t sSize);

  F_SQRT() = default;

  /**
   * @brief Calculates the embedding function value using a square-root dependence.
   * @param rho Electron density.
   * @param p Index for model parameters. Here p=0
   * @return Embedding function value \f$ F(\rho) = s_{p} \sqrt{ c_{p} \rho } \f$.
   */
  double calc_F(double rho, size_t p) const override;

  /**
   * @brief Calculates the derivative of the square-root embedding function.
   * @param rho Electron density.
   * @param p Index for model parameters. Here p=0
   * @return First derivative \f$ dF/d\rho = \frac{s_{p} \, {2 \sqrt{\rho}} \f$.
   */
  double calc_dF(double rho, size_t p) const override;
};

} // namespace models
} // namespace tadah

#endif
