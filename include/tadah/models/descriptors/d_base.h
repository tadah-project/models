#ifndef D_BASE_H
#define D_BASE_H

#include <string>
#include <optional>
#include <tadah/core/context.h>
#include <tadah/core/core_types.h>
#include <tadah/core/periodic_table.h>
#include <tadah/core/logger.h>
#include <tadah/core/registry.h>
#include <tadah/models/cut_all.h>
#include <tadah/models/descriptors/param_meta.h>
#include <vector>
#include <bitset>

namespace tadah {
namespace models {

class Bitset2D {
  static constexpr size_t N = 119;
  std::bitset<N * N> data;

public:
  void init(size_t i, size_t j) {
    data[i * N + j] = true;
    data[j * N + i] = true;

  }
  void uninit(size_t i, size_t j) {
    data[i * N + j] = false;
    data[j * N + i] = false;
  }
  bool is_init(size_t i, size_t j) const {
    return data[i * N + j];
  }
};

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
  std::vector<bool> is_optimizable;
  static void get_grid(const tadah::core::Context &c, const std::string &key, tadah::core::v_type &v);
  static tadah::core::v_type get_grid(const std::vector<std::string>&);  // expand joined grids
  int verbose;
  virtual ~D_Base();

  // return position of an argument given to TYPExB for a given key
  // if key is not used, returns -1
  int get_arg_pos(const std::string &key) const;


  /** \brief Return label of this descriptor.
   */
  virtual std::string label() const=0;
  double weights[119] = {1};  // ignore zero index; w[1]->H
  D_Base();
  D_Base(const tadah::core::Context &c);

  /** Refresh per-element weights[] from WATOM in the (possibly updated)
   *  context. Intended for HPO reinit paths where WATOM may change
   *  between iters but the descriptor object is reused. */
  void refresh_atom_weights(const tadah::core::Context &c);

  /** Resolve a 1-based WATOM arg_pos (matching `OPTIM WATOM (n)`) to the
   *  atomic number Z of the n-th ATOM entry in the context that last
   *  populated `weights[]`. Returns 0 when arg_pos is out of range or
   *  the ATOM list was never populated (no context constructor / no
   *  refresh_atom_weights call). Used by calc_dDdtheta /
   *  calc_dXijdri_dtheta / calc_drhoi_dtheta WATOM branches to test
   *  whether a per-pair contribution belongs to the optimised species. */
  int get_atom_Z_for_arg_pos(int arg_pos) const;

  /** \brief Print grid to the log.
   *
   * @param v Vector to print.
   * @param key Key name for the vector.
   */
  void printGrid(core::v_type &v, const std::string &key, core::LogLevel logLevel=core::LogLevel::INFO);

  virtual void set_fidx(size_t fidx_);
  virtual size_t get_fidx() const;

  /**
 * @brief Initialize the descriptor.
 *
 * This method sets up the required keys names and number of parameters for the descriptor.
 *
 * Note: The implementation of the descriptor should not rely on these initial values,
 * allowing for greater flexibility in how the tadah::core::Context file is parsed.
 *
 * The keys and parameters are utilized by other classes inheriting from this implementation,
 * such as meta-descriptors.
 *
 */
  virtual void init() = 0;

  virtual void set_fcut(Cut_Base* cut, bool manage_memory);

  double get_rcut() const;

  /** \brief Return dimension of the descriptor.
  */
  virtual size_t size() const;

  /** \brief Update descriptor parameters in-place from context.
   *
   * Re-reads all parameters that were baked in at construction time.
   * Returns true if the descriptor's output size changed, requiring a
   * full Trainer rebuild.  Default: no-op (for file-based or fixed descriptors).
   */
  virtual bool reinit(const tadah::core::Context& /*ctx*/) { return false; }

  virtual bool is_init_for_atoms(int Zi, int Zj) const;
  virtual void init_for_atoms(int Zi, int Zj);
  virtual void init_for_atoms(const std::vector<std::string> &Zs);
  virtual void uninit_for_atoms(int Zi, int Zj);
  static std::vector<std::string> get_init_atoms(const tadah::core::Context &c, std::string type);
  virtual std::vector<std::string> get_init_atoms(const tadah::core::Context &c)=0;

  /** \brief Suggested physically-motivated default for an OPTIM-able
   * parameter belonging to this descriptor.
   *
   * Used by the HPO PRIOR initialisation strategy. Returns std::nullopt
   * (the base-class default) when no prior is known for the (key,
   * arg_pos) pair; HPO will then fall back to the geometric/arithmetic
   * mid-point of the [low, high] bounds.
   *
   * Concrete descriptors are encouraged to override this whenever the
   * literature suggests a sensible default for a given parameter (e.g.
   * typical RCUT for a given coordination shell, typical SIGMA width
   * for a Gaussian basis).
   *
   * @param key      OPTIM key (e.g. "RCUT2B", "SGRID2B", "CEMBFUNC").
   * @param arg_pos  Numeric ordinal of the parameter within the key's
   *                 value list (1-based, matching the OPTIM (n) syntax).
   *                 Pass -1 when the descriptor is scalar in this key.
   */
  virtual std::optional<double>
  suggested_value(const std::string& /*key*/, int /*arg_pos*/) const {
    return std::nullopt;
  }

  /** \brief Per-descriptor parameter metadata for the HPO bound-validation
   *  pipeline (audit, suggester, per-iter health monitor).
   *
   *  Default `std::nullopt` means "no metadata supplied" — the audit
   *  treats the descriptor as opaque (skips block-level health checks
   *  that depend on metadata) and `suggested_value` (above) remains the
   *  only optional hint mechanism. Concrete descriptors override this
   *  method to declare their parameter physics. See the merry-peacock
   *  plan §"Examples — what `meta()` actually looks like" for canonical
   *  shapes (D2_LJ scalar, D2_Blip vector, DM_mBlipEnv full).
   */
  virtual std::optional<DescriptorMeta> meta() const {
    return std::nullopt;
  }

private:

  bool manage_memory=false; // who owns fcut
  Bitset2D init_for_atoms_map;  // 0th is unused

  /* Z values in ATOM-list order, populated by refresh_atom_weights.
   * Empty when the descriptor was default-constructed (no Context).
   * `OPTIM WATOM (n)` resolves to `watom_Z_by_arg_pos[n-1]`. */
  std::vector<int> watom_Z_by_arg_pos;
};
} // namespace models
} // namespace tadah
#endif
