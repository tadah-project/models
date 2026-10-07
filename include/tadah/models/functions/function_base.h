#ifndef FUNCTION_BASE_h
#define FUNCTION_BASE_h

#include <tadah/core/core_types.h>
#include <tadah/core/context.h>

#include <cstddef>

namespace tadah {
namespace models {


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
  virtual double epredict(const tadah::core::t_type &, const tadah::core::aed_type& ) const=0;
  virtual double fpredict(const tadah::core::t_type &, const tadah::core::fd_type &,
                          const tadah::core::aed_type&, const size_t ) const=0;

  virtual tadah::core::force_type fpredict(const tadah::core::t_type &, const tadah::core::fd_type &,
                              const tadah::core::aed_type& ) const=0;

  /** \brief Number of model weights for a descriptor of @p dsize entries:
   *  the design-matrix width. BF_Linear and Kern_Linear: dsize;
   *  BF_Polynomial2: dsize (dsize + 1) / 2; a kernel with a basis: the
   *  number of basis vectors, whatever dsize is (0 before it has a basis). */
  virtual std::size_t weights_count(std::size_t dsize) const=0;

  /** \brief Number of descriptor entries of each basis vector the weights
   *  multiply (a kernel's basis.rows()); 0 for a function without a basis
   *  (the basis functions and Kern_Linear act on the descriptor itself). */
  virtual std::size_t basis_vector_size() const;
};
} // namespace models
} // namespace tadah
#endif
