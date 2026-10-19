#include "tadah/core/registry.h"
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <tadah/core/number_format.h>
#include <tadah/models/descriptors/d_mjoin.h>
#include <tadah/models/descriptors/d_base.h>
#include <tadah/models/descriptors/dm/embed_spec.h>
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
  if(type_str=="TYPEMB" && !c.blocks("DESCRIPTOR").empty()) return true;
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

// The grid of `key` (its values and LIN/GEOM generators) as values; a
// refusal of get_grid names the key, as D_Base::get_grid(Context) does.
static std::vector<std::string> expanded_grid(const std::vector<std::string> &g,
                                              const std::string &key) {
  try {
    return convert_to_strings(D_Base::get_grid(g));
  } catch (const std::invalid_argument &e) {
    throw std::invalid_argument(key + ": " + e.what());
  }
}

void D_mJoin::expand_grids(tadah::core::Context &c) {
  resolve_embed_blocks(c);
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
  if(type_str=="TYPEMB") resolve_embed_blocks(c);
}
std::vector<tadah::core::Context> D_mJoin::parse_context(const tadah::core::Context &c, std::string type_str) {
  if(type_str=="TYPEMB" && !c.blocks("DESCRIPTOR").empty()) return embed_contexts(c);
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
    cgrid = expanded_grid(cgrid, cgridkey);
    data_map[cgridkey] = {&cgrid, 0};
  }
  std::vector<std::string> sgrid;
  if (c.contains(sgridkey)) {
    sgrid.resize(c.size(sgridkey));
    c.get(sgridkey, sgrid);
    sgrid = expanded_grid(sgrid, sgridkey);
    data_map[sgridkey] = {&sgrid, 0};
  }
  std::vector<std::string> cembgrid;
  if (c.contains("CEMBFUNC")) {
    cembgrid.resize(c.size("CEMBFUNC"));
    c.get("CEMBFUNC", cembgrid);
    cembgrid = expanded_grid(cembgrid, "CEMBFUNC");
    data_map["CEMBFUNC"] = {&cembgrid, 0};
  }
  std::vector<std::string> sembgrid;
  if (c.contains("SEMBFUNC")) {
    sembgrid.resize(c.size("SEMBFUNC"));
    c.get("SEMBFUNC", sembgrid);
    sembgrid = expanded_grid(sembgrid, "SEMBFUNC");
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
      extra_lists[e] = expanded_grid(extra_lists[e], k);
      data_map[k] = {&extra_lists[e], 0};
    }
  }

  // The descriptors of the line, each with its tokens: its label, its
  // parameters and its element pairs, up to the next label.
  std::vector<std::vector<std::string>> lines;
  for (const auto &t : types) {
    if (lines.empty() || t.find(prefix) == 0) lines.emplace_back();
    lines.back().push_back(t);
  }

  // The descriptors of the line, for a message: "2 descriptors (D2_BP,
  // D2_BP)".
  std::string named = std::to_string(lines.size()) +
                      (lines.size() == 1 ? " descriptor (" : " descriptors (");
  for (std::size_t i = 0; i < lines.size(); ++i)
    named += (i ? ", " : "") + lines[i][0];
  named += ")";
  auto values = [](size_t v) {
    return std::to_string(v) + (v == 1 ? " value" : " values");
  };

  // One RCTYPE and one RCUT value per descriptor, in the order of the lines:
  // fewer would be read past their end, more would never be read.
  for (const std::string &key : {rctypekey, rcutkey}) {
    const std::size_t have = data_map[key].first->size();
    if (have == lines.size()) continue;
    throw std::runtime_error(
        type_str + " names " + named + " but " + key + " has " +
        values(have) + ": each descriptor takes one " + key +
        " value, in the order of the " + type_str + " lines");
  }

  std::vector<tadah::core::Context> contexts;
  // The list keys some descriptor of the line reads (SIZE 0 included), with
  // the values they take, for the check of leftover values below.
  std::map<std::string, size_t> read_keys;
  std::string last_text;
  for (size_t li = 0; li < lines.size(); ++li) {
    const auto &line = lines[li];
    const std::string &token = line[0];
    std::unique_ptr<D_Base> d(tadah::core::factory<D_Base>(token));
    // The line has the descriptor's parameters and its SIZE parameters are
    // integers (else an error naming the line), so they can be read below.
    const std::vector<std::string> elements = d->element_tokens(type_str, line);
    std::string text = type_str;
    for (const auto &t : line) text += " " + t;

    std::vector<std::string> all_list_keys = d->keys;
    all_list_keys.insert(all_list_keys.end(), d->data_keys.begin(),
                         d->data_keys.end());
    const size_t fparam = d->nparams - all_list_keys.size();
    std::vector<size_t> sizes;
    for (size_t k = 0; k < all_list_keys.size(); ++k) {
      const std::string &size_token = line[1 + fparam + k];
      try {
        sizes.push_back(std::stoul(size_token));
      } catch (const std::out_of_range &) {
        throw std::runtime_error(
            text + ": parameter " + std::to_string(fparam + k + 1) +
            ", the number of " + all_list_keys[k] + " values, is too large: " +
            size_token);
      }
    }

    // Every list key whose SIZE asks for values must be set.
    std::string takes, missing;
    size_t nmissing = 0;
    for (size_t k = 0; k < all_list_keys.size(); ++k) {
      const std::string &key = all_list_keys[k];
      if (sizes[k] == 0 || data_map.find(key) != data_map.end()) continue;
      ++nmissing;
      takes += (takes.empty() ? "" : " and ") + std::to_string(sizes[k]) +
               " " + key + (sizes[k] == 1 ? " value" : " values");
      missing += (missing.empty() ? "" : " and ") + key;
    }
    if (nmissing) {
      throw std::runtime_error(
          type_str + " " + token + " takes " + takes + ", but " + missing +
          (nmissing == 1 ? " is" : " are") + " not set: " + text);
    }

    tadah::core::Context c1 = c;
    c1.remove(type_str); c1.remove(rctypekey); c1.remove(rcutkey);
    c1.remove(cgridkey); c1.remove(sgridkey);
    c1.remove("CEMBFUNC"); c1.remove("SEMBFUNC");
    c1.remove("AMPGRIDMB");
    c1.remove("KDEG2B"); c1.remove("KDEGMB"); c1.remove("KDEGEMB");
    c1.add(type_str,token);
    c1.add(rctypekey, (*data_map[rctypekey].first)[ data_map[rctypekey].second++ ]);
    c1.add(rcutkey, (*data_map[rcutkey].first)[ data_map[rcutkey].second++ ]);
    for (size_t j=0; j<fparam; ++j) {
      c1.add(type_str, line[1 + j]);
    }
    for (size_t k = 0; k < all_list_keys.size(); ++k) {
      const std::string &key = all_list_keys[k];
      const size_t n = sizes[k];
      // The key (set when n > 0, checked above) must hold n more values:
      // those the descriptors before this one on the line took are theirs.
      read_keys[key] += n;
      const size_t given = n > 0 ? data_map.at(key).first->size() : 0;
      const size_t taken = n > 0 ? data_map.at(key).second : 0;
      if (n > 0 && taken + n > given) {
        std::string where = type_str + " " + token;
        if (lines.size() > 1)
          where += " (descriptor " + std::to_string(li + 1) + " of " +
                   std::to_string(lines.size()) + ")";
        const std::string has =
            taken == 0 ? values(given)
                       : values(given - taken) + " left of " +
                             std::to_string(given) +
                             " (the descriptors before it take " +
                             std::to_string(taken) + ")";
        throw std::runtime_error(where + " takes " + std::to_string(n) + " " +
                                 key + (n == 1 ? " value" : " values") +
                                 ", but " + key + " has " + has + ": " + text);
      }
      c1.add(type_str, n);
      for (size_t j = 0; j < n; ++j) {
        c1.add(key, (*data_map[key].first)[ data_map[key].second++ ]);
      }
    }

    // the element pairs
    if (elements.empty()) {
      throw std::runtime_error(
        "Error: No element types specified. "
        "Please provide element types for this calculator. "
        "For example: '"+token+" Ti Ti Ti Nb' sets the calculator to compute Ti-Ti and Ti-Nb interactions. "
        "Token: " + token);
    }

    if (elements.size() % 2 != 0) {
      throw std::runtime_error(
        "Error: Element types must be provided in pairs. "
        "Ensure each element is paired correctly. "
        "For instance, '" + token + " Ti * Nb Zr' means that interactions of Ti with any atom, as well as Nb-Zr, will be computed. "
        "Token: " + token);
    }
    for (const auto &e : elements) c1.add(type_str, e);

    contexts.push_back(c1);
    last_text = text;
  }

  // Every value of a list key the line reads belongs to a descriptor: values
  // past those the descriptors take would be dropped from the model and the
  // potential without a word.
  for (const auto &[key, taken] : read_keys) {
    const auto it = data_map.find(key);
    const size_t given = it == data_map.end() ? 0 : it->second.first->size();
    if (given <= taken) continue;
    const std::string takes = std::to_string(taken) + " " + key +
                              (taken == 1 ? " value" : " values");
    if (lines.size() == 1)
      throw std::runtime_error(type_str + " " + lines[0][0] + " takes " +
                               takes + ", but " + key + " has " +
                               values(given) + ": " + last_text);
    throw std::runtime_error(
        type_str + " names " + named + ", which take " + takes + ", but " +
        key + " has " + values(given) + ": each descriptor takes its " + key +
        " values, in the order of the " + type_str + " lines");
  }
  return contexts;
}

} // namespace models
} // namespace tadah

