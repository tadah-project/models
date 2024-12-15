#ifndef D_BASE_H
#define D_BASE_H

#include <tadah/core/config.h>
#include <tadah/core/core_types.h>
#include <tadah/core/periodic_table.h>
#include <tadah/core/registry.h>
#include <tadah/models/cut_all.h>

/** \brief Base class for all descriptor types.
 *
 */
class D_Base {
protected:
  size_t fidx=0;    // first index in aed and fd arrays for this descriptor
  size_t s=0;
  Cut_Base *fcut=nullptr;

public:
  static void get_grid(const Config &c, std::string key, v_type &v);
  int verbose;
  virtual ~D_Base();

  /** \brief Return dimension of the descriptor.
  */
  virtual size_t size() { return s; };

  /** \brief Return label of this descriptor.
   */
  virtual std::string label()=0;
  double weights[119] = {1};  // ignore zero index; w[1]->H
  D_Base();
  D_Base(Config &c);

  virtual void set_fidx(size_t fidx_);
  virtual size_t get_fidx();

  virtual void set_fcut(Cut_Base* cut, bool manage_memory);

  double get_rcut() { 
    return fcut->get_rcut(); }

private:
  bool manage_memory=false; // who owns fcut
};
#endif
