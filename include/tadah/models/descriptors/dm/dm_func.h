#ifndef DM_func_H
#define DM_func_H

#include <tadah/core/config.h>
#include <tadah/models/descriptors/d_base.h>

#include <cstdlib>

class F_Base {
public:
  /**
   * @brief Calculates the embedding function value.
   * @param rho Electron density.
   * @param p Index for model parameters.
   * @return Function value F(rho).
   */
  virtual double calc_F(double rho, size_t p) = 0;

  /**
   * @brief Calculates the first derivative of the embedding function.
   * @param rho Electron density.
   * @param p Index for model parameters.
   * @return Derivative dF/d(rho).
   */
  virtual double calc_dF(double rho, size_t p) = 0;
};

/**
 * @class F_RLR
 * @brief Implements an embedding function of the form: \f$ s \rho \log(c \rho) \f$
 *
 * This class supports embedding functions characterized by two main parameters:
 * - **SEMBFUNC**: Controls the depth, \f$ s \f$, of the embedding function. 
 * - **CEMBFUNC**: Determines the x-intercept, with the x-intercept at \f$ 1/c \f$.
 * 
 *  Require: size(**SEMBFUNC**)=size(**CEMBFUNC**)=1
 * 
 * The number of keys for these parameters must match the entries in the mEAD descriptor.
 *
 */
class F_RLR:  public F_Base {
private:
  Config *config=nullptr;
  v_type sgrid;   // controls depth
  v_type cgrid;   // control x-intercep at 1/c
public:
  F_RLR(Config &conf);
  F_RLR();
  double calc_F(double rho,size_t p) override;
  double calc_dF(double rho,size_t p) override;
  //~F_RLR() {
  // Do not delete config!
  //}
};

/**
 * @class F_SQRT
 * @brief Implements \f$ s \sqrt{\rho} \f$.
 *
 * Optional parameter:
 * - **SEMBFUNC**: Controls the strength, \f$ s \f$, of the embedding function.
 *   If no; value is provided, the default is 1.
 *
 *  Require: size(**SEMBFUNC**)=(0 or 1) and size(**CEMBFUNC**)=0
 *
 */
class F_SQRT : public F_Base {
private:
  Config *config = nullptr; ///< Configuration object (not owned by this class).
  v_type sgrid;             ///< Controls depth or scaling factor.

public:
  /**
   * @brief Constructs F_SQRT from a Config object.
   * @param conf Reference to a Config instance.
   *
   * The constructor extracts or initializes relevant parameter for the model from
   * the provided configuration.  SEMBFUNC is set to 1 if not specified in the Config.
   */
  F_SQRT(Config &conf);

  /**
   * @brief Default constructor for F_SQRT.
   *
   * This version does not initialize from a Config object. Users must
   * manually set sgrid if this constructor is used. SEMBFUNC is set to 1.
   */
  F_SQRT();

  /**
   * @brief Calculates the embedding function value using a square-root dependence.
   * @param rho Electron density.
   * @param p Index for model parameters. Here p=0
   * @return Embedding function value \f$ F(\rho) = s_{p} \sqrt{ c_{p} \rho } \f$.
   */
  double calc_F(double rho, size_t p) override;

  /**
   * @brief Calculates the derivative of the square-root embedding function.
   * @param rho Electron density.
   * @param p Index for model parameters. Here p=0
   * @return First derivative \f$ dF/d\rho = \frac{s_{p} \, {2 \sqrt{\rho}} \f$.
   */
  double calc_dF(double rho, size_t p) override;
};


#endif
