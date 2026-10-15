#include "catch2/catch.hpp"
#include <tadah/core/core_user_config_parser.h>
#include <tadah/core/periodic_table.h>
#include <tadah/models/dc_selector.h>
#include <tadah/models/descriptors/d_basis_functions.h>
#include <memory>
#include <sstream>

using namespace tadah::core;
using namespace tadah::models;
namespace {
Context embed_config(std::string density = "Blip5", std::string embedding = "Gaussian",
                     std::string form = "NONNEG", std::string cutoff = "Cut_Poly",
                     std::string envelope = "OFF", bool widths = false) {
  PeriodicTable::initialize();
  const bool dk = density == "Knot5", ek = embedding == "Knot5";
  std::string text = "VERBOSE 0\nINIT2B false\nINITMB true\nATOM Ta\n"
      "MODEL M_BLR BF_Linear\nDESCRIPTOR ta DM_Embed\nPAIRS Ta Ta\n"
      "DENSITY\nBASIS " + density + "\nFORM " + form + "\nCENTRES 1 2\n" +
      (dk ? "DEGREES 3 5\n" : widths ? "WIDTHS 0.5 0.5\n" : "ETA 2 2\n") +
      "AMPLITUDES 0.4 0.2\nCUTOFF\nTYPE " + cutoff +
      "\nRADIUS 6\nENDCUTOFF\nENDDENSITY\nEMBEDDING\nBASIS " + embedding +
      "\nCENTRES 0 0.5 1\n" + (ek ? "DEGREES 5 3 5\n" :
      widths ? "WIDTHS 0.5 0.5 0.5\n" : "ETA 2 2 2\n") +
      "REFERENCE ZERO\nENVELOPE\nTYPE " + envelope + "\n" +
      (envelope == "OFF" ? "" : "CUTOFF 2\n") +
      "ENDENVELOPE\nENDEMBEDDING\nENDDESCRIPTOR\n";
  std::istringstream in(text);
  return CoreUserConfigParser().parse(in);
}
double energy(DM_Base &d, double r) {
  rho_type rho; d.init_rhoi(rho);
  d.calc_rho(73,73,r,r*r,Vec3d(r,0,0),rho);
  aed_type a(d.size()); a.set_zero(); d.calc_aed(rho,a);
  double e=0; for(std::size_t i=0;i<d.size();++i) e+=(i+1)*a(i);
  return e;
}
}

TEST_CASE("DM_Embed constructs through the descriptor selector", "[unified-dm]") {
  auto c = embed_config();
  DC_Selector s(c);
  REQUIRE(s.dmb->label() == "DM_Embed");
  REQUIRE(s.dmb->size() == 3);
  REQUIRE(s.dmb->eam_form());
  REQUIRE(s.dmb->get_rcut() == 6.0);
  rho_type rho; s.dmb->init_rhoi(rho);
  aed_type a(3); a.set_zero(); s.dmb->calc_aed(rho,a);
  for(std::size_t m=0;m<3;++m) REQUIRE(a(m) == 0.0);
}

TEST_CASE("DM_Embed independently selects density embedding cutoff and envelope",
          "[unified-dm][matrix]") {
  const std::vector<std::string> families{
    "Blip","Blip5","Gaussian","Knot5","WendlandC2","WendlandC4"};
  for(const auto &d:families) for(const auto &e:families)
    for(const auto &cut:{"Cut_Cos","Cut_CosS","Cut_Poly","Cut_PolyS","Cut_Tanh","Cut_PT","NONE"})
      for(const auto &env:{"OFF","SinInv","Smoothstep2","Smoothstep4"}) {
        CAPTURE(d,e,cut,env);
        auto c=embed_config(d,e,"NONNEG",cut,env);
        if(d=="Gaussian" && std::string(cut)=="NONE") {
          REQUIRE_THROWS(DC_Selector(c)); continue;
        }
        DC_Selector s(c);
        const double r=1.37,h=1e-5;
        const double fd=(energy(*s.dmb,r+h)-energy(*s.dmb,r-h))/(2*h);
        rho_type rho; s.dmb->init_rhoi(rho);
        s.dmb->calc_rho(73,73,r,r*r,Vec3d(r,0,0),rho);
        aed_type a(3); a.set_zero(); s.dmb->calc_aed(rho,a);
        fd_type f(3); f.set_zero();
        s.dmb->calc_dXijdri(73,73,r,r*r,Vec3d(r,0,0),rho,f);
        double value=0; for(std::size_t m=0;m<3;++m) value+=(m+1)*f(m,0);
        REQUIRE(value == Approx(fd).margin(2e-7));
      }
}

TEST_CASE("DM_Embed force buffers are additive and invalid reinit preserves state",
          "[unified-dm]") {
  auto c=embed_config(); DC_Selector s(c);
  const double before=energy(*s.dmb,1.3);
  auto bad=c;
  auto blocks=bad.blocks();
  for(auto &line:blocks[0].lines) if(line.find("AMPLITUDES")!=std::string::npos)
    line="AMPLITUDES -1 1";
  bad.replace_blocks(blocks);
  REQUIRE_THROWS(s.dmb->reinit(bad));
  REQUIRE(energy(*s.dmb,1.3) == before);
  rho_type rho; s.dmb->init_rhoi(rho);
  s.dmb->calc_rho(73,73,1.3,1.69,Vec3d(1.3,0,0),rho);
  aed_type a(3); a.set_zero(); s.dmb->calc_aed(rho,a);
  fd_type f(3),g(3); f.set_zero(); g.set_zero();
  s.dmb->calc_dXijdri(73,73,1.3,1.69,Vec3d(1.3,0,0),rho,f);
  for(std::size_t m=0;m<3;++m) g(m,0)=7;
  s.dmb->calc_dXijdri(73,73,1.3,1.69,Vec3d(1.3,0,0),rho,g);
  for(std::size_t m=0;m<3;++m) REQUIRE(g(m,0)==Approx(7+f(m,0)));
}
