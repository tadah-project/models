#ifndef FUNCTION_BASE_h
#define FUNCTION_BASE_h

#include <tadah/core/core_types.h>
#include <tadah/core/context.h>

#include <string>
#include <vector>

namespace tadah {
namespace models {

/** \brief One value of \ref MPARAMS that a function reads
 *  (Function_Base::mparams()). */
struct MParamSlot {
  std::string name;       ///< the parameter, e.g. "sigma" (Kern_RBF) or "d"
  bool integer = false;   ///< read as an integer (an exact integer value)
  bool positive = false;  ///< a searched value (OPTIM MPARAMS) must be > 0
};

/** Base class for Kernels and Basis Functions */
struct Function_Base {

private:
  int verbose;

public:
  Function_Base();
  Function_Base(const tadah::core::Context &c);
  void set_verbose(int v);
  int get_verbose();

  // Derived classes must implement Derived() and Derived(tadah::core::Context)
  virtual ~Function_Base();

  // For Kernels
  virtual double operator() (const tadah::core::aed_type& , const tadah::core::aed_type& ) const;
  virtual tadah::core::aed_type derivative(const tadah::core::aed_type& , const tadah::core::aed_type& ) const;
  virtual double prime(const tadah::core::aed_type& , const tadah::core::aed_type& ,
                       const tadah::core::aed_type& ) const;
  virtual void set_basis(const tadah::core::Matrix );

  /** \brief Return label of this object */
  virtual std::string get_label() const=0;

  /** \brief The \ref MPARAMS values this function reads when it is built,
   *  in order: element i is MPARAMS (i+1).
   *
   *  - Kern_RBF: sigma (positive);
   *  - Kern_Polynomial: d (an integer, positive), gamma, c;
   *  - Kern_Sigmoid: gamma, c;
   *  - every other function (Kern_Linear, Kern_Quadratic, Kern_LQ,
   *    BF_Linear, BF_Polynomial2): none, the default.
   *
   *  tadah hpo refuses at startup an OPTIM MPARAMS (i) of a function that
   *  reads fewer than i values, of an integer value, and a box that is not
   *  positive for a positive one. */
  virtual std::vector<MParamSlot> mparams() const { return {}; }

  /** \brief The number of leading columns of this function's design
   *  matrix that are intercepts, columns 0 .. n-1.
   *
   *  An intercept column is constant over the atoms: every atom adds the
   *  same value to the energy row, and its force and stress rows are 0.
   *  The regression leaves these columns out of the LAMBDA ridge penalty
   *  (LinearRegressor::train), and OALGO 3, which penalises every column
   *  alike, refuses a model with LAMBDA != 0 that has one.
   *
   *  The descriptor has a bias intercept when BIAS is true and ATOM is
   *  set (tadah::core::n_intercept_columns() > 0); its value is 1 for
   *  every atom.
   *  - BF_Linear, Kern_Linear: tadah::core::n_intercept_columns(), the
   *    bias columns of the descriptor;
   *  - BF_Polynomial2: 1 with a bias intercept (column 0, bias x bias),
   *    else 0;
   *  - Kern_Quadratic, Kern_LQ, Kern_Polynomial, Kern_Sigmoid (functions
   *    of b . x): 1 with a bias intercept when the basis vector of
   *    column 0 has no component other than the bias (the vector e0 that
   *    training puts there): K(b0, x) then depends on x only through the
   *    bias, which is 1; else 0;
   *  - Kern_RBF: 0 (K(e0, x) depends on every component of x).
   *
   *  Pure, so that a new function states its own. */
  virtual size_t intercept_columns(const tadah::core::Context &c) const=0;
  virtual double epredict(const tadah::core::t_type &, const tadah::core::aed_type& ) const=0;
  virtual double fpredict(const tadah::core::t_type &, const tadah::core::fd_type &,
                          const tadah::core::aed_type&, const size_t ) const=0;

  virtual tadah::core::force_type fpredict(const tadah::core::t_type &, const tadah::core::fd_type &,
                              const tadah::core::aed_type& ) const=0;
};
} // namespace models
} // namespace tadah
#endif
