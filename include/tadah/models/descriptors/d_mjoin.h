#ifndef D_MJOIN_H
#define D_MJOIN_H
#include <tadah/core/config.h>
#include <vector>

/** \brief Base class to be used by meta descriptors for joining
*/
class D_mJoin {
protected:
  struct D_Type {
    std::string type;
    std::string rc_type;
    double rcut;
    std::vector<std::string> params;
    Config c;
  };
  std::vector<D_Type> dt;
  std::string lab = "D_mJoin";
public:
  virtual ~D_mJoin();
  D_mJoin();
  virtual std::vector<Config> parse_config(Config &c, std::string type_str);
};
#endif

