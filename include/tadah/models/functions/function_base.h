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
  Function_Base(const Config &c);
  void set_verbose(int v);
  int get_verbose();

  // Derived classes must implement Derived() and Derived(Config)
  virtual ~Function_Base();

  // For Kernels
  virtual double operator() (const aed_type& , const aed_type& ) const;
  virtual aed_type derivative(const aed_type& , const aed_type& ) const;
  virtual double prime(const aed_type& , const aed_type& ,
                       const aed_type& ) const;
  virtual void set_basis(const Matrix );

  /** \brief Return label of this object */
  virtual std::string get_label() const=0;
  virtual double epredict(const t_type &, const aed_type& ) const=0;
  virtual double fpredict(const t_type &, const fd_type &,
                          const aed_type&, const size_t ) const=0;
} // namespace models
} // namespace tadah

  virtual force_type fpredict(const t_type &, const fd_type &,
                              const aed_type& ) const=0;
};
#endif
