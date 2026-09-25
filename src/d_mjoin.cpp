#include "tadah/core/registry.h"
#include <cstddef>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <tadah/core/number_format.h>
#include <tadah/models/descriptors/d_mjoin.h>
#include <tadah/models/descriptors/d_base.h>
#include <vector>
namespace tadah {
namespace models {

namespace {
// parse_context() splits the standard grid keys (and the declared
// data_keys, e.g. the knot-degree lists). A descriptor whose keys fall
// outside the supported set is left untouched by expand_grids; it reads its
// own grids straight from the context. mJoin meta-descriptors are always
// handled (their sub-descriptors are standard). Keyless descriptors (D2_EAM,
// DM_EAM) are trivially handled.
bool mjoin_can_handle(const tadah::core::Context &c, const std::string &type_str) {
  if (!c.contains(type_str) || c.size(type_str) == 0) return false;
  const std::string label = c.get<std::string>(type_str, 0);
  if (label.find("mJoin") != std::string::npos) return true;
  static const std::set<std::string> sup2b = {"CGRID2B", "SGRID2B",
                                              "KDEG2B"};
  static const std::set<std::string> supmb = {"CGRIDMB", "SGRIDMB",
                                              "CEMBFUNC", "SEMBFUNC",
                                              "AMPGRIDMB", "KDEGMB",
                                              "KDEGEMB"};
  const std::set<std::string> &sup = (type_str == "TYPE2B") ? sup2b : supmb;
  std::unique_ptr<D_Base> d(tadah::core::factory<D_Base>(label));
  for (const auto &k : d->keys)
    if (!sup.count(k)) return false;
  for (const auto &k : d->data_keys)
    if (!sup.count(k)) return false;
  return true;
}
} // namespace


D_mJoin::D_mJoin() {}
D_mJoin::~D_mJoin() {}

// Lossless: each value is written as its shortest round-trip text, so the
// sub-descriptor contexts see bit-identical grids. (A fixed NUMERIC-decimal
// rendering here used to round every expanded grid value, and flush values
// below 10^-NUMERIC to zero, on every mJoin construction and reinit.)
std::vector<std::string> convert_to_strings(const tadah::core::v_type& vec) {
  std::vector<std::string> result;
  result.reserve(vec.size());
  for (double num : vec) {
    result.push_back(tadah::core::format_roundtrip(num));
  }
  return result;
}

void D_mJoin::expand_grids(tadah::core::Context &c) {
  if (mjoin_can_handle(c, "TYPE2B")) {
    std::vector<tadah::core::Context> contexts1 = parse_context(c, "TYPE2B");
    expand_grids(c,contexts1,"TYPE2B");
  }
  if (mjoin_can_handle(c, "TYPEMB")) {
    std::vector<tadah::core::Context> contexts2 = parse_context(c, "TYPEMB");
    expand_grids(c,contexts2,"TYPEMB");
  }
}
void D_mJoin::expand_grids(tadah::core::Context &c, std::string type_str) {
  std::vector<tadah::core::Context> contexts = parse_context(c, type_str);
  expand_grids(c,contexts,type_str);
}
void D_mJoin::expand_grids(tadah::core::Context &c, std::vector<tadah::core::Context> &contexts,std::string type_str) {
  if (type_str=="TYPE2B") {
    c.remove("CGRID2B"); c.remove("SGRID2B");
    for (auto &c1 : contexts) {
      if (c1.contains("CGRID2B"))
        for (size_t i = 0; i < c1.size("CGRID2B"); ++i)
          c.add("CGRID2B", c1.get<double>("CGRID2B", i));
      if (c1.contains("SGRID2B"))
        for (size_t i = 0; i < c1.size("SGRID2B"); ++i)
          c.add("SGRID2B", c1.get<double>("SGRID2B", i));
    }
  }
  else if (type_str=="TYPEMB") {
    c.remove("CGRIDMB"); c.remove("SGRIDMB");
    c.remove("CEMBFUNC"); c.remove("SEMBFUNC");

    for (auto &c1 : contexts) {
      if (c1.contains("CGRIDMB"))
        for (size_t i = 0; i < c1.size("CGRIDMB"); ++i)
          c.add("CGRIDMB", c1.get<double>("CGRIDMB", i));
      if (c1.contains("SGRIDMB"))
        for (size_t i = 0; i < c1.size("SGRIDMB"); ++i)
          c.add("SGRIDMB", c1.get<double>("SGRIDMB", i));
      if (c1.contains("CEMBFUNC"))
        for (size_t i = 0; i < c1.size("CEMBFUNC"); ++i)
          c.add("CEMBFUNC", c1.get<double>("CEMBFUNC", i));
      if (c1.contains("SEMBFUNC"))
        for (size_t i = 0; i < c1.size("SEMBFUNC"); ++i)
          c.add("SEMBFUNC", c1.get<double>("SEMBFUNC", i));
    }
  }
  else {
    std::cerr << "Cannot expand context grid: wrong keyword: " << type_str << std::endl;
  }
}
std::vector<tadah::core::Context> D_mJoin::parse_context(const tadah::core::Context &c, std::string type_str) {
  std::string prefix = type_str=="TYPE2B" ? "D2_" : "DM_";
  std::string rctypekey = "RC" + type_str;
  std::string rcutkey = "RCUT" + type_str.substr(type_str.size() - 2);
  std::string cgridkey = "CGRID" + type_str.substr(type_str.size() - 2);
  std::string sgridkey = "SGRID" + type_str.substr(type_str.size() - 2);

  std::map<std::string, std::pair<const std::vector<std::string>*, size_t>> data_map;

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
  data_map[rctypekey]= { &rctypes, 0 };

  std::vector<std::string> rcut(c.size(rcutkey));
  c.get(rcutkey, rcut);
  data_map[rcutkey]={&rcut, 0};

  std::vector<std::string> cgrid;
  if (c.contains(cgridkey)) {
    cgrid.resize(c.size(cgridkey));
    c.get(cgridkey, cgrid);
    cgrid = convert_to_strings(D_Base::get_grid(cgrid));
    data_map[cgridkey] = {&cgrid, 0};
  }
  std::vector<std::string> sgrid;
  if (c.contains(sgridkey)) {
    sgrid.resize(c.size(sgridkey));
    c.get(sgridkey, sgrid);
    sgrid = convert_to_strings(D_Base::get_grid(sgrid));
    data_map[sgridkey] = {&sgrid, 0};
  }
  std::vector<std::string> cembgrid;
  if (c.contains("CEMBFUNC")) {
    cembgrid.resize(c.size("CEMBFUNC"));
    c.get("CEMBFUNC", cembgrid);
    cembgrid = convert_to_strings(D_Base::get_grid(cembgrid));
    data_map["CEMBFUNC"] = {&cembgrid, 0};
  }
  std::vector<std::string> sembgrid;
  if (c.contains("SEMBFUNC")) {
    sembgrid.resize(c.size("SEMBFUNC"));
    c.get("SEMBFUNC", sembgrid);
    sembgrid = convert_to_strings(D_Base::get_grid(sembgrid));
    data_map["SEMBFUNC"] = {&sembgrid, 0};
  }
  // additional splittable lists: density amplitudes (optimizable key) and
  // the knot-degree data_keys (split-only, never OPTIM)
  std::vector<std::string> extra_lists[4];
  const char *extra_names[4] = {"AMPGRIDMB", "KDEG2B", "KDEGMB", "KDEGEMB"};
  for (int e = 0; e < 4; ++e) {
    const std::string k = extra_names[e];
    if (c.contains(k)) {
      extra_lists[e].resize(c.size(k));
      c.get(k, extra_lists[e]);
      extra_lists[e] = convert_to_strings(D_Base::get_grid(extra_lists[e]));
      data_map[k] = {&extra_lists[e], 0};
    }
  }

  std::vector<tadah::core::Context> contexts;
  for (auto it = types.begin(); it != types.end();) {
    tadah::core::Context c1 = c;
    c1.remove(type_str); c1.remove(rctypekey); c1.remove(rcutkey);
    c1.remove(cgridkey); c1.remove(sgridkey);
    c1.remove("CEMBFUNC"); c1.remove("SEMBFUNC");
    c1.remove("AMPGRIDMB");
    c1.remove("KDEG2B"); c1.remove("KDEGMB"); c1.remove("KDEGEMB");
    const auto token = *it++;
    c1.add(type_str,token);
    c1.add(rctypekey, (*data_map[rctypekey].first)[ data_map[rctypekey].second++ ]);
    c1.add(rcutkey, (*data_map[rcutkey].first)[ data_map[rcutkey].second++ ]);
    D_Base *d = tadah::core::factory<D_Base>(token);
    size_t fparam = d->nparams - d->keys.size() - d->data_keys.size();
    for (size_t j=0; j<fparam; ++j) {
      c1.add(type_str,*it++);
    }
    std::vector<std::string> all_list_keys = d->keys;
    all_list_keys.insert(all_list_keys.end(), d->data_keys.begin(),
                         d->data_keys.end());
    for (const auto& k : all_list_keys) {
      int n = std::stoi(*it++);
      // check dimensions: ensure there are at least n remaining tokens
      if (data_map.find(k) != data_map.end() && data_map[k].second + n > data_map[k].first->size()) {
        throw std::runtime_error(
            "Error: Not enough data for key '" + k + "' in context. "
            "Expected " + std::to_string(n) + " but got " + std::to_string(data_map[k].first->size()) +
            ". Token: " + token);
      }
      c1.add(type_str, n);
      for (int j = 0; j < n; ++j) {
        c1.add(k, (*data_map[k].first)[ data_map[k].second++ ]);
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

    contexts.push_back(c1);
    delete d;
  }
  return contexts;
}

} // namespace models
} // namespace tadah

