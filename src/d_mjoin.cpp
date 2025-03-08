#include "tadah/core/registry.h"
#include <cstddef>
#include <iostream>
#include <string>
#include <tadah/models/descriptors/d_mjoin.h>
#include <tadah/models/descriptors/d_base.h>
#include <vector>
namespace tadah {
namespace models {


D_mJoin::D_mJoin() {}
D_mJoin::~D_mJoin() {}

std::vector<std::string> convert_to_strings(const tadah::core::v_type& vec, int precision) {
  std::vector<std::string> result;
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(precision);

  for (double num : vec) {
    oss.str("");   // Clear the stream
    oss << num;
    result.push_back(oss.str());
  }

  return result;
}

void D_mJoin::expand_grids(tadah::core::Config &c) {
  if (c.exist("TYPE2B")) {
    std::vector<tadah::core::Config> configs1 = parse_config(c, "TYPE2B");
    expand_grids(c,configs1,"TYPE2B");
  }
  if (c.exist("TYPEMB")) {
    std::vector<tadah::core::Config> configs2 = parse_config(c, "TYPEMB");
    expand_grids(c,configs2,"TYPEMB");
  }
}
void D_mJoin::expand_grids(tadah::core::Config &c, std::string type_str) {
  std::vector<tadah::core::Config> configs = parse_config(c, type_str);
  expand_grids(c,configs,type_str);
}
void D_mJoin::expand_grids(tadah::core::Config &c, std::vector<tadah::core::Config> &configs,std::string type_str) {
  if (type_str=="TYPE2B") {
    c.remove("CGRID2B"); c.remove("SGRID2B");
    for (auto &c1 : configs) {
      if (c1.exist("CGRID2B"))
        for (size_t i = 0; i < c1.size("CGRID2B"); ++i)
          c.add("CGRID2B", c1.get<double>("CGRID2B", i));
      if (c1.exist("SGRID2B"))
        for (size_t i = 0; i < c1.size("SGRID2B"); ++i)
          c.add("SGRID2B", c1.get<double>("SGRID2B", i));
    }
  }
  else if (type_str=="TYPEMB") {
    c.remove("CGRIDMB"); c.remove("SGRIDMB");
    c.remove("CEMBFUNC"); c.remove("SEMBFUNC");

    for (auto &c1 : configs) {
      if (c1.exist("CGRIDMB"))
        for (size_t i = 0; i < c1.size("CGRIDMB"); ++i)
          c.add("CGRIDMB", c1.get<double>("CGRIDMB", i));
      if (c1.exist("SGRIDMB"))
        for (size_t i = 0; i < c1.size("SGRIDMB"); ++i)
          c.add("SGRIDMB", c1.get<double>("SGRIDMB", i));
      if (c1.exist("CEMBFUNC"))
        for (size_t i = 0; i < c1.size("CEMBFUNC"); ++i)
          c.add("CEMBFUNC", c1.get<double>("CEMBFUNC", i));
      if (c1.exist("SEMBFUNC"))
        for (size_t i = 0; i < c1.size("SEMBFUNC"); ++i)
          c.add("SEMBFUNC", c1.get<double>("SEMBFUNC", i));
    }
  }
  else {
    std::cerr << "Cannot expand config grid: wrong keyword: " << type_str << std::endl;
  }
}
std::vector<tadah::core::Config> D_mJoin::parse_config(tadah::core::Config &c, std::string type_str) {
  std::string prefix = type_str=="TYPE2B" ? "D2_" : "DM_";
  std::string rctypekey = "RC" + type_str;
  std::string rcutkey = "RCUT" + type_str.substr(type_str.size() - 2);
  std::string cgridkey = "CGRID" + type_str.substr(type_str.size() - 2);
  std::string sgridkey = "SGRID" + type_str.substr(type_str.size() - 2);

  std::map<std::string, std::vector<std::string>::iterator> data_map;

  std::vector<std::string> types(c.size(type_str)); // including params TYPE2B D2_Blip 7 7
  c.get(type_str, types);
  //types.erase(types.begin()); // remove D2_mJoin or DM_mJoin
  std::string pattern="mJoin";
  types.erase(
    std::remove_if(types.begin(), types.end(), [&pattern](const std::string& s) {
      return s.find(pattern) != std::string::npos;
    }),
    types.end()
  );
  
  std::vector<std::string> rctypes(c.size(rctypekey));
  c.get(rctypekey, rctypes);
  data_map[rctypekey]=rctypes.begin();

  std::vector<std::string> rcut(c.size(rcutkey));
  c.get(rcutkey, rcut);
  data_map[rcutkey]=rcut.begin();

  size_t p = c.get<size_t>("OUTPREC");
  std::vector<std::string> cgrid;
  if (c.exist(cgridkey)) {
    cgrid.resize(c.size(cgridkey));
    c.get(cgridkey, cgrid);
    cgrid = convert_to_strings(D_Base::get_grid(cgrid),p); 
    data_map[cgridkey] = cgrid.begin();
  }
  std::vector<std::string> sgrid;
  if (c.exist(sgridkey)) {
    sgrid.resize(c.size(sgridkey));
    c.get(sgridkey, sgrid);
    sgrid = convert_to_strings(D_Base::get_grid(sgrid),p); 
    data_map[sgridkey] = sgrid.begin();
  }
  std::vector<std::string> cembgrid;
  if (c.exist("CEMBFUNC")) {
    cembgrid.resize(c.size("CEMBFUNC"));
    c.get("CEMBFUNC", cembgrid);
    cembgrid = convert_to_strings(D_Base::get_grid(cembgrid),p); 
    data_map["CEMBFUNC"] = cembgrid.begin();
  }
  std::vector<std::string> sembgrid;
  if (c.exist("SEMBFUNC")) {
    sembgrid.resize(c.size("SEMBFUNC"));
    c.get("SEMBFUNC", sembgrid);
    sembgrid = convert_to_strings(D_Base::get_grid(sembgrid),p); 
    data_map["SEMBFUNC"] = sembgrid.begin();
  }

  std::vector<tadah::core::Config> configs;
  for (auto it = types.begin(); it != types.end();) {
    tadah::core::Config c1 = c;
    c1.remove(type_str); c1.remove(rctypekey); c1.remove(rcutkey);
    c1.remove(cgridkey); c1.remove(sgridkey);
    const auto token = *it++;
    c1.add(type_str,token);
    c1.add(rctypekey,*(data_map[rctypekey]++));
    c1.add(rcutkey,*(data_map[rcutkey]++));
    D_Base *d = tadah::core::factory<D_Base>(token);
    size_t fparam = d->nparams - d->keys.size();
    for (size_t j=0; j<fparam; ++j) {
      c1.add(type_str,*it++);
    }

    for (const auto& k : d->keys) {
      int n = std::stoi(*it++);
      c1.add(type_str,n);
      for (int j = 0; j < n; ++j) {
        c1.add(k,*(data_map[k]++));
      } 
    }

    // read init atoms
    size_t found_init_atoms=0;
    while (it != types.end()) {
      if (it->find(prefix) == 0) {
        break;
      }
      c1.add(type_str,*it++);
      found_init_atoms++;
    }
    // if (it != types.end() && it != types.begin()) {
    //   --it;
    // }
    if (!found_init_atoms) {
      throw std::runtime_error(
        "Error: No element types specified. "
        "Please provide element types for this calculator. "
        "For example: '"+token+" Ti Ti Ti Nb' sets the calculator to compute Ti-Ti and Ti-Nb interactions. "
        "Token: " + token);
    }

    if (found_init_atoms % 2 != 0) {
      throw std::runtime_error(
        "Error: Element types must be provided in pairs. "
        "Ensure each element is paired correctly. "
        "For instance, '" + token + " Ti * Nb Zr' means that interactions of Ti with any atom, as well as Nb-Zr, will be computed. "
        "Token: " + token);
    }

    configs.push_back(c1);
    delete d;
  }
  return configs;
}

} // namespace models
} // namespace tadah

