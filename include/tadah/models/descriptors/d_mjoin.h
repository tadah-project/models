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
  static std::vector<tadah::core::Context> parse_context(tadah::core::Context &c, std::string type_str);

  static void expand_grids(tadah::core::Context &c, std::vector<tadah::core::Context> &contexts, std::string type_str);

  static void expand_grids(tadah::core::Context &c, std::string type_str);
  static void expand_grids(tadah::core::Context &c);
};
} // namespace models
} // namespace tadah
#endif

