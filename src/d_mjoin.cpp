#include "tadah/core/registry.h"
#include <cstddef>
#include <string>
#include <tadah/models/descriptors/d_mjoin.h>
#include <tadah/models/descriptors/d_base.h>
#include <vector>

D_mJoin::D_mJoin() {}
D_mJoin::~D_mJoin() {}

std::vector<std::string> convert_to_strings(const v_type& vec, int precision) {
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

std::vector<Config> D_mJoin::parse_config(Config &c, std::string type_str) {
  std::string prefix = type_str=="TYPE2B" ? "D2_" : "DM_";
  std::string rctypekey = "RC" + type_str;
  std::string rcutkey = "RCUT" + type_str.substr(type_str.size() - 2);
  std::string cgridkey = "CGRID" + type_str.substr(type_str.size() - 2);
  std::string sgridkey = "SGRID" + type_str.substr(type_str.size() - 2);

  std::map<std::string, std::vector<std::string>::iterator> data_map;

  std::vector<std::string> types(c.size(type_str)); // including params TYPE2B D2_Blip 7 7
  c.get(type_str, types);
  types.erase(types.begin()); // remove D2_mJoin or DM_mJoin
  
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

  std::vector<Config> configs;
  for (auto it = types.begin(); it != types.end();) {
    Config c1 = c;
    c1.remove(type_str); c1.remove(rctypekey); c1.remove(rcutkey);
    c1.remove(cgridkey); c1.remove(sgridkey);
    const auto token = *it++;
    c1.add(type_str,token);
    c1.add(rctypekey,*(data_map[rctypekey]++));
    c1.add(rcutkey,*(data_map[rcutkey]++));
    D_Base *d = CONFIG::factory<D_Base>(token);
    size_t fparam = d->nparams - d->keys.size();
    for (size_t j=0; j<fparam; ++j) {
      c1.add(type_str,*it++);
    }

    for (const auto& k : d->keys) {
      int n = std::stoi(*it++);
      for (int j = 0; j < n; ++j) {
        c1.add(k,*(data_map[k]++));
      } 
    }

    configs.push_back(c1);
    delete d;
  }
  return configs;
}
