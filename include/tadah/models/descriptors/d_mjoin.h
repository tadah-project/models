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
   *  counted after LIN/GEOM expansion), and when RCTYPE or RCUT has fewer
   *  values than the line has descriptors.
   */
  static std::vector<tadah::core::Context> parse_context(const tadah::core::Context &c, std::string type_str);

  static void expand_grids(tadah::core::Context &c, std::vector<tadah::core::Context> &contexts, std::string type_str);

  static void expand_grids(tadah::core::Context &c, std::string type_str);
  static void expand_grids(tadah::core::Context &c);
};
} // namespace models
} // namespace tadah
#endif

