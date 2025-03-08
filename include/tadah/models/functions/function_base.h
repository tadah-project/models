#ifndef FUNCTION_BASE_h
#define FUNCTION_BASE_h

#include <tadah/core/core_types.h>
#include <tadah/core/config.h>

namespace tadah {
namespace models {


/** Base class for Kernels and Basis Functions */
struct Function_Base {

private:
  int verbose;

public:
  Function_Base();
  Function_Base(const tadah::core::Config &c);
  void set_verbose(int v);
  int get_verbose();

  // Derived classes must implement Derived() and Derived(tadah::core::Config)
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
};
} // namespace models
} // namespace tadah
#endif
