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
  // ---- the descriptor key contract (mJoin + HPO addressing) ----
  // `keys`: the descriptor's OPTIMIZABLE per-instance list parameters.
  //   Every entry is (a) splittable by the *_mJoin meta-descriptors and
  //   (b) addressable by HPO via get_arg_pos/OPTIM.
  // `data_keys`: per-instance list DATA that mJoin must split but that is
  //   NOT a continuous optimisation dimension (e.g. the knot-degree lists
  //   KDEG*: a discrete 3|5 structure choice). Never exposed to OPTIM.
  // TYPE2B/TYPEMB token layout (consumed by D_mJoin::expand_grids and
  // get_init_atoms):
  //   <Label> <fparam family params> <one SIZE token per `keys` entry>
  //           <one SIZE token per `data_keys` entry> <El El ...>
  // with nparams = fparam + keys.size() + data_keys.size() and
  // get_arg_pos(key) = index-in-keys + fparam. A size token of 0 means the
  // pot carries no entries (and no key) for that channel — e.g. the knot
  // families' width channels, where eta is a pure gauge.
  std::vector<std::string> keys;       // optimizable list parameters
  std::vector<std::string> data_keys;  // split-only list data (no OPTIM)
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

  /** Symmetric per-pair WATOM weights. `weights[Zi][Zj] == weights[Zj][Zi]`.
   *  Filled to 1.0 in every constructor; populated by update_weights() from
   *  either per-species (geometric mean) or pair-triple WATOM input.
   *  Index 0 is unused (no element has Z=0). */
  double weights[119][119];

  /** Convenience read-only accessor for the symmetric pair weight. */
  inline double get_weight(int Zi, int Zj) const { return weights[Zi][Zj]; }

  /** Tracks which WATOM input form populated `weights`. Drives arg_pos→Z
   *  vs arg_pos→pair resolution in dweight_dw_n. */
  enum class WatomMode { Empty, PerSpecies, PairTriple };
  WatomMode watom_mode = WatomMode::Empty;

  D_Base();
  D_Base(const tadah::core::Context &c);

  /** Populate `weights[][]` from the (possibly updated) WATOM key in `c`.
   *  Accepts two input formats:
   *    - Per-species (parallel to ATOM): pair weights computed as the
   *      geometric mean `weights[Zi][Zj] = sqrt(w_i * w_j)`, diagonal
   *      `weights[Zi][Zi] = w_i`. Sets watom_mode = PerSpecies and
   *      populates watom_Z_by_arg_pos.
   *    - Pair triples (sym sym value)+ : symmetric direct writes. Sets
   *      watom_mode = PairTriple and populates watom_pair_by_arg_pos.
   *  Always resets `weights` to 1.0 first when called via the
   *  constructor; subsequent calls (HPO reinit) write only the listed
   *  pairs so callers can rely on prior state being preserved on
   *  unlisted slots. */
  void update_weights(const tadah::core::Context &c);

  /** Resolve a 1-based per-species WATOM arg_pos (matching
   *  `OPTIM WATOM (n)`) to the atomic number Z of the n-th ATOM entry
   *  populated by the last update_weights call. Returns 0 when arg_pos
   *  is out of range, the ATOM list was never populated, or the
   *  descriptor is in PairTriple mode. */
  int get_atom_Z_for_arg_pos(int arg_pos) const;

  /** Chain factor `d weights[Zi][Zj] / d θ_{arg_pos}` for the current
   *  WATOM mode.
   *    PerSpecies : geometric-mean derivative w.r.t. per-element w_n;
   *      returns 1.0 when Zi==Zj==Z_target (diagonal direct write),
   *      0.5 * sqrt(w_other / w_n) when exactly one of Zi,Zj equals
   *      Z_target, 0 otherwise.
   *    PairTriple : 1.0 iff (Zi,Zj) matches the arg_pos-th triple's
   *      pair (either order), 0 otherwise.
   *    Empty / out-of-range arg_pos : 0. */
  double dweight_dw_n(int Zi, int Zj, int arg_pos) const;

  /** Returns true when `arg_pos` is in range for the current watom_mode
   *  (i.e. the descriptor owns the WATOM derivative for that θ). Lets
   *  descriptor calc_dDdtheta branches distinguish "I own this θ but
   *  contribute zero for this (Zi,Zj) pair" (return true) from "this θ
   *  is unconfigured, fall back to FD" (return false). */
  bool has_watom_arg_pos(int arg_pos) const;

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

  /** \brief Per-descriptor parameter metadata for the HPO bound-validation
   *  pipeline (audit, suggester, per-iter health monitor).
   *
   *  Default `std::nullopt` means "no metadata supplied" — the audit
   *  treats the descriptor as opaque (skips block-level health checks
   *  that depend on metadata). Concrete descriptors override this
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

  /* Z values in ATOM-list order, populated by update_weights when the
   * per-species WATOM input form is detected.
   * Empty when watom_mode != PerSpecies (default-constructed, pair-triple
   * input, or no WATOM key). `OPTIM WATOM (n)` resolves to
   * `watom_Z_by_arg_pos[n-1]` only in PerSpecies mode. */
  std::vector<int> watom_Z_by_arg_pos;

  /* (Zi, Zj) pairs in the order they appear in WATOM, populated by
   * update_weights when the pair-triple WATOM input form is detected.
   * Empty when watom_mode != PairTriple. `OPTIM WATOM (k)` resolves to
   * watom_pair_by_arg_pos[k-1] only in PairTriple mode. */
  std::vector<std::pair<int,int>> watom_pair_by_arg_pos;
};
} // namespace models
} // namespace tadah
#endif
