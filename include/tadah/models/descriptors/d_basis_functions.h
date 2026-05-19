#ifndef D_BASIS_F_H
#define D_BASIS_F_H

#include <math.h>

namespace tadah {
namespace models {

/** \defgroup BasisFunctions Radial Basis Functions
 * \brief Inline Gaussian and Blip radial basis functions and their derivatives.
 *
 * This header provides two families of radial basis functions used to expand
 * atomic densities and two-body descriptors:
 *
 * - **Gaussian** functions \f$ G(r;\eta,\mu) \f$
 * - **Blip** (B-spline) functions \f$ \mathcal{B}(x_c) \f$
 *
 * Both families are available in plain and cutoff-weighted overloads.
 * @{
 */

// =========================================================================
// GAUSSIANS
// =========================================================================

/** \brief Gaussian basis function.
 *
 * \f[
 *   G(r;\,\eta,\mu) = e^{-\eta(r-\mu)^2}
 * \f]
 *
 * \param r   interatomic distance \f$ r_{ij} \f$
 * \param eta width parameter \f$ \eta > 0 \f$
 * \param miu centre \f$ \mu \f$ (often written \f$ r_s \f$)
 */
inline double G(double r, double eta, double miu) {
  return exp(-eta*(r-miu)*(r-miu));
}

/** \brief Gaussian basis function weighted by a cutoff.
 *
 * \f[
 *   G(r;\,\eta,\mu)\,f_c(r) = e^{-\eta(r-\mu)^2}\,f_c(r)
 * \f]
 *
 * \param r   interatomic distance
 * \param eta width parameter \f$ \eta \f$
 * \param miu centre \f$ \mu \f$
 * \param f   cutoff value \f$ f_c(r) \f$
 */
inline double G(double r, double eta, double miu, double f) {
  return exp(-eta*(r-miu)*(r-miu))*f;
}

/** \brief Derivative of the Gaussian basis function w.r.t.\ \f$ r \f$.
 *
 * \f[
 *   \frac{d G}{d r} = -2\eta(r-\mu)\,e^{-\eta(r-\mu)^2}
 * \f]
 *
 * \param r   interatomic distance
 * \param eta width parameter \f$ \eta \f$
 * \param miu centre \f$ \mu \f$
 */
inline double dG(double r, double eta, double miu) {
  return -2.0*eta*(r-miu)*exp(-eta*(r-miu)*(r-miu));
}

/** \brief Derivative of the cutoff-weighted Gaussian w.r.t.\ \f$ r \f$.
 *
 * Product rule applied to \f$ G(r)\,f_c(r) \f$:
 * \f[
 *   \frac{d}{dr}\bigl[G(r;\eta,\mu)\,f_c(r)\bigr]
 *   = e^{-\eta(r-\mu)^2}\!\left[-2\eta(r-\mu)\,f_c(r)+f_c'(r)\right]
 * \f]
 *
 * \param r   interatomic distance
 * \param eta width parameter \f$ \eta \f$
 * \param miu centre \f$ \mu \f$
 * \param f   cutoff value \f$ f_c(r) \f$
 * \param fp  cutoff derivative \f$ f_c'(r) \f$
 */
inline double dG(double r, double eta, double miu, double f, double fp) {
  return exp(-eta*(r-miu)*(r-miu) )*(-2.0*f*eta*(r-miu)+fp);
}

// =========================================================================
// BLIPS
// =========================================================================

/** \brief Blip (cubic B-spline) basis function.
 *
 * The blip function is a piecewise cubic polynomial with compact support on
 * \f$(-2,\,2)\f$ in the normalised coordinate \f$ x_c = \eta(r_{ij}-r_s) \f$:
 *
 * \f[
 *   \mathcal{B}(x_c) =
 *   \begin{cases}
 *     \dfrac{1}{4}(2+x_c)^3       & -2 < x_c \le -1 \\[6pt]
 *     1 - \dfrac{3}{2}x_c^2 + \dfrac{3}{4}|x_c|^3 & -1 < x_c \le 1 \\[6pt]
 *     \dfrac{1}{4}(2-x_c)^3       & 1 < x_c < 2 \\[4pt]
 *     0                           & \text{otherwise}
 *   \end{cases}
 * \f]
 *
 * The function is \f$ C^1 \f$ everywhere and \f$ C^2 \f$ everywhere except at
 * \f$ x_c = \pm 1 \f$.
 *
 * \param x_c normalised coordinate \f$ x_c = \eta(r_{ij}-r_s) \f$
 *
 * \sa Hernández, E., Gillan, M., Goringe, C. (1997).
 *     Basis functions for linear-scaling first-principles calculations.
 *     *Phys. Rev. B*, **55**(20), 13485–13493.
 *     https://doi.org/10.1103/PhysRevB.55.13485
 */
inline double B(double x_c)
{
  if (-2.0 < x_c && x_c <= -1.0) {
    double b = 2.0 + x_c;
    return 0.25 * b*b*b;
  }
  else if (-1.0 < x_c && x_c <= 1.0) {
    double x2 = x_c*x_c;
    double x3 = x2*x_c;
    return 1.0 - 1.5*x2 + 0.75*std::fabs(x3);
  }
  else if (1.0 < x_c && x_c < 2.0) {
    double b = 2.0 - x_c;
    return 0.25 * b*b*b;
  }
  else {
    return 0.0;
  }
}

/** \brief Blip basis function evaluated at \f$ x_c = \eta(r-r_s) \f$.
 *
 * Convenience wrapper: computes \f$ x_c \f$ from physical coordinates and
 * delegates to B(double).
 *
 * \param r   interatomic distance
 * \param eta inverse width \f$ \eta \f$; full support in \f$ r \f$ is \f$ 4/\eta \f$
 * \param miu centre \f$ r_s \f$
 */
inline double B(double r, double eta, double miu) {
  return B(eta*(r-miu));
}

/** \brief Cutoff-weighted blip basis function.
 *
 * \f[
 *   \mathcal{B}(x_c)\,f_c(r)
 * \f]
 *
 * \param x_c normalised coordinate \f$ x_c = \eta(r_{ij}-r_s) \f$
 * \param f   cutoff value \f$ f_c(r) \f$
 */
inline double B(double x_c, double f)
{
  if (-2.0 < x_c && x_c <= -1.0) {
    double b = 2.0 + x_c;
    return 0.25 * b*b*b *f;
  }
  else if (-1.0 < x_c && x_c <= 1.0) {
    double x2 = x_c*x_c;
    double x3 = x2*x_c;
    return f*(1.0 - 1.5*x2 + 0.75*std::fabs(x3));
  }
  else if (1.0 < x_c && x_c < 2.0) {
    double b = 2.0 - x_c;
    return f*(0.25 * b*b*b);
  }
  else {
    return 0.0;
  }
}

/** \brief Derivative of the blip basis function w.r.t.\ \f$ r \f$.
 *
 * Given \f$ x_c = \eta(r-r_s) \f$, by the chain rule:
 * \f[
 *   \frac{d\mathcal{B}}{dr} = \eta\,\frac{d\mathcal{B}}{dx_c}
 * \f]
 * where
 * \f[
 *   \frac{d\mathcal{B}}{dx_c} =
 *   \begin{cases}
 *     \dfrac{3}{4}(2+x_c)^2             & -2 < x_c \le -1 \\[6pt]
 *     -3x_c + \dfrac{9}{4}\,x_c|x_c|   & -1 < x_c \le 1  \\[6pt]
 *     -\dfrac{3}{4}(2-x_c)^2            & 1 < x_c < 2     \\[4pt]
 *     0                                  & \text{otherwise}
 *   \end{cases}
 * \f]
 *
 * \param x_c normalised coordinate \f$ x_c = \eta(r-r_s) \f$
 * \param eta inverse width \f$ \eta \f$
 */
inline double dB(double x_c, double eta)
{
  if (-2.0 < x_c && x_c <= -1.0) {
    double d = 2.0+x_c;
    return 0.75*eta*d*d;
  }
  else if (-1.0 < x_c && x_c <= 1.0) {
    return -3.0*eta*x_c + 2.25*eta*x_c*std::fabs(x_c);
  }
  else if (1.0 < x_c && x_c < 2.0) {
    double d = 2.0-x_c;
    return -0.75*eta*d*d;
  }
  else {
    return 0.0;
  }
}

/** \brief Derivative of the cutoff-weighted blip function w.r.t.\ \f$ r \f$.
 *
 * Product rule applied to \f$ f_c(r)\,\mathcal{B}(x_c) \f$:
 * \f[
 *   \frac{d}{dr}\bigl[f_c(r)\,\mathcal{B}(x_c)\bigr]
 *   = f_c'(r)\,\mathcal{B}(x_c) + f_c(r)\,\frac{d\mathcal{B}}{dr}
 * \f]
 *
 * \param x_c normalised coordinate \f$ x_c = \eta(r-r_s) \f$
 * \param eta inverse width \f$ \eta \f$
 * \param f   cutoff value \f$ f_c(r) \f$
 * \param fp  cutoff derivative \f$ f_c'(r) \f$
 */
inline double dB(double x_c, double eta, double f, double fp)
{
  if (-2.0 < x_c && x_c <= -1.0) {
    double d = 2.0+x_c;
    return f*(0.75*eta*d*d)+fp*(0.25 * d*d*d);
  }
  else if (-1.0 < x_c && x_c <= 1.0) {
    return f*(-3.0*eta*x_c + 2.25*eta*x_c*std::fabs(x_c))
      +fp*(1.0 - 1.5*x_c*x_c + 0.75*std::fabs(x_c*x_c*x_c));
  }
  else if (1.0 < x_c && x_c < 2.0) {
    double d = 2.0-x_c;
    return f*(-0.75*eta*d*d)+fp*(0.25 * d*d*d);
  }
  else {
    return 0.0;
  }
}

/** \brief Derivative of the blip function evaluated at physical coordinates.
 *
 * Convenience wrapper: computes \f$ x_c = \eta(r-r_s) \f$ and calls
 * dB(double, double).
 *
 * \param r   interatomic distance
 * \param eta inverse width \f$ \eta \f$
 * \param miu centre \f$ r_s \f$
 */
inline double dB(double r, double eta, double miu) {
  return dB(eta*(r-miu), eta);
}

/** \brief Second derivative of the blip basis function w.r.t.\ \f$ x_c \f$.
 *
 * Piecewise quadratic — discontinuous at \f$ |x_c| = 1 \f$ (B is C^1
 * globally but only C^2 outside |x_c|=1). Used by HPO analytical-gradient
 * code paths that need the F''(rho) chain through Blip embeddings.
 *
 * \param x_c normalised coordinate \f$ x_c = \eta(r-r_s) \f$
 */
inline double d2B(double x_c) {
  if (-2.0 < x_c && x_c <= -1.0) return 1.5 * (2.0 + x_c);
  if (-1.0 < x_c && x_c <= 1.0)  return -3.0 + 4.5 * std::fabs(x_c);
  if (1.0 < x_c && x_c < 2.0)    return 1.5 * (2.0 - x_c);
  return 0.0;
}

/** @} */ // end of BasisFunctions group

} // namespace models
} // namespace tadah
#endif
