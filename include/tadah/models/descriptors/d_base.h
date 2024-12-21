#ifndef D_BASE_H
#define D_BASE_H

#include <string>
#include <tadah/core/config.h>
#include <tadah/core/core_types.h>
#include <tadah/core/periodic_table.h>
#include <tadah/core/registry.h>
#include <tadah/models/cut_all.h>
#include <vector>

/** \brief Base class for all descriptor types.
 *
 */
class D_Base {
protected:
  size_t fidx=0;    // first index in aed and fd arrays for this descriptor
  size_t s=0;
  Cut_Base *fcut=nullptr;

public:
  std::vector<std::string> keys;    // keys required by this descriptor
  size_t nparams;  // number of params which follow TYPExB
  static void get_grid(const Config &c, std::string key, v_type &v);
  static v_type get_grid(std::vector<std::string>&);  // expand joined grids
  int verbose;
  virtual ~D_Base();

  // return position of an argument given to TYPExB for a given key
  // if key is not used, returns -1
  int get_arg_pos(const std::string &key) const;

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

  /**
 * @brief Initialize the descriptor.
 *
 * This method sets up the required keys names and number of parameters for the descriptor.
 *
 * Note: The implementation of the descriptor should not rely on these initial values,
 * allowing for greater flexibility in how the Config file is parsed.
 *
 * The keys and parameters are utilized by other classes inheriting from this implementation,
 * such as meta-descriptors.
 *
 */
  virtual void init() = 0;

  virtual void set_fcut(Cut_Base* cut, bool manage_memory);

  double get_rcut() { 
    return fcut->get_rcut(); }

private:
  bool manage_memory=false; // who owns fcut
};
#endif
