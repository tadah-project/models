#ifndef D_MJOIN_H
#define D_MJOIN_H
#include <tadah/core/context.h>
#include <vector>

namespace tadah {
namespace models {

/** \brief Base class to be used by meta descriptors for joining
*/
class D_mJoin {
protected:
  struct D_Type {
    std::string type;
    std::string rc_type;
    double rcut;
    std::vector<std::string> params;
    tadah::core::Context c;
  };
  std::vector<D_Type> dt;
  std::string lab = "D_mJoin";
public:
  virtual ~D_mJoin();
  D_mJoin();
  /** \brief One context per descriptor of the @p type_str line (TYPE2B or
   *  TYPEMB) of @p c, each with that descriptor's TYPE line, its RCTYPE and
   *  RCUT value and the list values (CGRID*, SGRID*, AMPGRIDMB, CEMBFUNC,
   *  SEMBFUNC, KDEG*) its SIZE parameters ask for, taken in order.
   *
   *  \throws std::runtime_error naming the descriptor's line when the line
   *  lacks the descriptor's parameters or a SIZE parameter is not a
   *  non-negative integer (D_Base::element_tokens), when a list key that a
   *  SIZE parameter asks values of is not set ('TYPE2B D2_BP takes 2 CGRID2B
   *  values and 2 SGRID2B values, but CGRID2B and SGRID2B are not set:
   *  TYPE2B D2_BP 2 2 Ta Ta'; a SIZE of 0 needs no key) or has fewer
   *  values left than the SIZE parameter asks for (on a joined line the
   *  descriptors before it take theirs first: 'TYPE2B D2_BP (descriptor 2
   *  of 2) takes 2 CGRID2B values, but CGRID2B has 1 value left of 3 (the
   *  descriptors before it take 2): TYPE2B D2_BP 2 2 Ta Ta'; the values are
   *  counted after LIN/GEOM expansion), when RCTYPE or RCUT has fewer or
   *  more values than the line has descriptors, and when a list key of the
   *  line's type has more values than the descriptors take (a key that no
   *  descriptor of the line reads included: 'TYPE2B D2_LJ takes 0 CGRID2B
   *  values, but CGRID2B has 3 values: TYPE2B D2_LJ Ta Ta').
   *
   *  A potential (@p c with WEIGHTS) is read as it was written: values its
   *  model does not read are reported once per key and process
   *  (core::report_ignored_key_once: 'Part of RCUT2B (6) is ignored: TYPE2B
   *  names 2 descriptors (D2_Blip, D2_Knot5), which take its first 2
   *  values; the potential works unchanged without it') instead of
   *  refused. Those are the RCTYPE, RCUT and list values past the ones a
   *  joined line's descriptors take, the RCTYPE values past the first and
   *  the RCUT values but the largest (its cutoff) of a single descriptor,
   *  and the values of a key no descriptor of the line reads. A single
   *  descriptor reads every value of the list keys it declares: their
   *  count is checked as in a configuration.
   *
   *  The list keys a descriptor reads whole (D_Base::context_keys, the
   *  CEMBFUNC/SEMBFUNC of a DM_Poly channel function) go whole, expanded,
   *  to a single line's descriptor. A joined line gives such a descriptor
   *  none; one that cannot run without them is refused ('TYPEMB DM_PS_RLR
   *  (descriptor 2 of 2) reads CEMBFUNC and SEMBFUNC whole and cannot be
   *  joined: ...').
   */
  static std::vector<tadah::core::Context> parse_context(const tadah::core::Context &c, std::string type_str);

  /** \brief The checks of parse_context() on the TYPE2B and TYPEMB lines
   *  of @p c that hold one descriptor (D2_mJoin and DM_mJoin check their
   *  own lines), for a model built from @p c (DC_Selector): a potential
   *  read by tadah predict or pair_style tadah reports the values its
   *  descriptors do not read, as a joined line's does.
   */
  static void check_single_lines(const tadah::core::Context &c);

  static void expand_grids(tadah::core::Context &c, std::vector<tadah::core::Context> &contexts, std::string type_str);

  static void expand_grids(tadah::core::Context &c, std::string type_str);
  /** \brief The TYPE2B and TYPEMB lines of @p c split and checked
   *  (parse_context()), and their grids expanded in @p c, for tadah train
   *  and tadah hpo. For a potential (@p c with WEIGHTS: training from a
   *  potential) the values its descriptors do not read are reported and
   *  left out of @p c, so the training context builds the potential's
   *  model: a single descriptor keeps its largest RCUT value, its cutoff.
   */
  static void expand_grids(tadah::core::Context &c);
};
} // namespace models
} // namespace tadah
#endif

