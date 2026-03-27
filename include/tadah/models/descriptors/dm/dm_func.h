#ifndef DM_func_H
#define DM_func_H

#include <tadah/core/context.h>
#include <tadah/models/descriptors/d_base.h>

#include <cstdlib>

namespace tadah {
namespace models {

class F_Base {
public:
  /** \brief True only for \f$\mathcal{F}(\rho)=\rho^2\f$ (F_SQ).
   *
   * DM_mEAD applies \f$\mathcal{F}\f$ to individual orbital densities
   * \f$\rho^{\eta,r_s,l_x,l_y,l_z}\f$ before the rotational-invariance
   * restoring quadratic sum.  This is invariant only when
   * \f$\mathcal{F}=\rho^2\f$; any other choice breaks SO(3) invariance for
   * \f$L>0\f$.
   */
  static constexpr bool rotationally_invariant = false;

  /**
   * @brief Whether this embedding function expands dimensions.
   *
   * When true, the embedding function takes a single rho value (from
   * radial grid of size 1) and produces multiple output channels
   * (one per embedding grid point).  DM_mEAD uses this to adjust
   * its inner loops accordingly.
   */
  static constexpr bool expands_dimensions = false;

  /**
   * @brief Returns the number of output channels per orbital per L.
   *
   * For most embedding functions this equals the radial grid size (ng).
   * For expanding functions (e.g. F_Blip) this equals the embedding
   * function grid size.
   */
  virtual size_t output_dim() const = 0;

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

  virtual ~F_Base() = default;
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
  size_t output_dim() const override final { return sfuncgrid.size(); }
  double calc_F(double rho,size_t p) const override final;
  double calc_dF(double rho,size_t p) const override final;
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

  size_t output_dim() const override final { return sfuncgrid.size(); }

  /**
   * @brief Calculates the embedding function value using a square-root dependence.
   * @param rho Electron density.
   * @param p Index for model parameters. Here p=0
   * @return Embedding function value \f$ F(\rho) = s_{p} \sqrt{ c_{p} \rho } \f$.
   */
  double calc_F(double rho, size_t p) const override final;

  /**
   * @brief Calculates the derivative of the square-root embedding function.
   * @param rho Electron density.
   * @param p Index for model parameters. Here p=0
   * @return First derivative \f$ dF/d\rho = \frac{s_{p} \, {2 \sqrt{\rho}} \f$.
   */
  double calc_dF(double rho, size_t p) const override final;
};

class F_SQ : public F_Base {
private:
  const tadah::core::Context *context=nullptr;
  size_t out_dim = 0;

public:
  /** F_SQ = ρ² is the only channel function that preserves rotational
   *  invariance in DM_mEAD for all L. */
  static constexpr bool rotationally_invariant = true;

  F_SQ(const tadah::core::Context &context, size_t, size_t);

  F_SQ() = default;

  size_t output_dim() const override final { return out_dim; }

  double calc_F(double rho, size_t p) const override final;

  double calc_dF(double rho, size_t p) const override final;
};

  /**
   * @class F_Blip
   * @brief Implements a "blip" embedding function.
   *
   * F_Blip is a dimension-expanding embedding function:
   * it requires CGRID/SGRID of size 1 (a single radial basis function
   * producing one ρ per orbital) and CEMBFUNC/SEMBFUNC of size ≥ 1.
   * Each (orbital, embedding-grid-point) pair produces one descriptor
   * component, so the descriptor dimension per L is CEMBFUNC.size().
   *
   * @param context Reference to a Context instance.
   */
class F_Blip:  public F_Base {
private:
  const tadah::core::Context *context=nullptr;
  tadah::core::v_type sfuncgrid;   // controls width of the "blip"
  tadah::core::v_type cfuncgrid;   // control blip offset
public:
  static constexpr bool expands_dimensions = true;

  F_Blip(const tadah::core::Context &context, size_t cSize, size_t sSize);
  F_Blip() = default;
  size_t output_dim() const override final { return sfuncgrid.size(); }
  double calc_F(double rho,size_t p) const override final;
  double calc_dF(double rho,size_t p) const override final;

};

} // namespace models
} // namespace tadah

#endif
