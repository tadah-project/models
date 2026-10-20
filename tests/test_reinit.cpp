/* Regression tests for D_Base::reinit() on meta-descriptors and for
 * M_Core::reset_trained().
 *
 * Bug 1 (HPO predict mismatch): D2_mJoin::reinit() and DM_mJoin::reinit()
 * updated each child's fcut but left the container's own outer fcut
 * (D2_mJoin / DM_mJoin base-class fcut) stale.  DescriptorsCalc::calc()
 * uses d2.get_rcut() and dm.get_rcut() as a pre-filter before calling into
 * the descriptor.  When RCUT2BMAX / RCUTMBMAX grew during HPO, the prefilter
 * was frozen at the initial value, so atom pairs beyond the original cutoff
 * were silently excluded from training and validation — but included by
 * tadah-predict (which reconstructs everything fresh), giving a descriptor
 * mismatch and catastrophically wrong predictions.
 *
 * Bug 2 (re-train guard): M_Core::trained is set to true after the first
 * training call.  Without an explicit reset, a later train() of the same
 * model is refused, so a model retrained in place must call
 * model->reset_trained() first.  (D_Base::reinit() and DC_Selector::reinit()
 * remain the in-place update of descriptor parameters; tadah hpo no longer
 * uses them: it builds a Trainer per evaluation.)
 */

#include "catch2/catch.hpp"
#include <tadah/core/core_types.h>
#include <tadah/models/descriptors/d_all.h>
#include <tadah/models/cutoffs.h>
#include <tadah/models/dc_selector.h>
#include <tadah/models/m_core.h>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace tadah::core;
using namespace tadah::models;

struct ReinitFixture {
    ReinitFixture() { PeriodicTable::initialize(); }
};

/* ------------------------------------------------------------------
 * Helpers
 * ------------------------------------------------------------------ */

// Context for DM_mJoin with two DM_Blip children.
// child1: Lmax=lmax1, cgrid/sgrid count=n1
// child2: Lmax=lmax2, cgrid/sgrid count=n2
// Grid values are 0.5, 1.0, ..., 0.5*(n1+n2) for cgrid; all 0.5 for sgrid.
// rcut: cutoff applied to both children (RCUTMBMAX = rcut).
static Context make_dm_mjoin_ctx(int lmax1, int n1, int lmax2, int n2,
                                  double rcut = 10.0)
{
    Context ctx;
    // RCUTMBMAX is computed automatically from RCUTMB by Context::add().
    ctx.add("RCTYPEMB", "Cut_Cos");
    ctx.add("RCTYPEMB", "Cut_Cos");
    ctx.add("RCUTMB", rcut);
    ctx.add("RCUTMB", rcut);
    ctx.remove("INITMB");
    ctx.add("INITMB", "true");
    // DM_mJoin inline format: mJoin marker, then per-child specs.
    // DM_Blip has nparams=3, keys={"CGRIDMB","SGRIDMB"}, fparam=1 (Lmax).
    // Per-child: "DM_Blip" <Lmax> <n_cgrid> <n_sgrid> <atom> <atom>
    ctx.add("TYPEMB", "DM_mJoin");
    ctx.add("TYPEMB", "DM_Blip");
    ctx.add("TYPEMB", lmax1); ctx.add("TYPEMB", n1); ctx.add("TYPEMB", n1);
    ctx.add("TYPEMB", "H"); ctx.add("TYPEMB", "H");
    ctx.add("TYPEMB", "DM_Blip");
    ctx.add("TYPEMB", lmax2); ctx.add("TYPEMB", n2); ctx.add("TYPEMB", n2);
    ctx.add("TYPEMB", "H"); ctx.add("TYPEMB", "H");
    for (int i = 0; i < n1 + n2; ++i) ctx.add("CGRIDMB", 0.5 * (i + 1));
    for (int i = 0; i < n1 + n2; ++i) ctx.add("SGRIDMB", 0.5);
    return ctx;
}

// Context for D2_mJoin with two D2_Blip children.
// D2_Blip has nparams=2, keys={"CGRID2B","SGRID2B"}, fparam=0.
// Per-child: "D2_Blip" <n_cgrid> <n_sgrid> <atom> <atom>
// rcut: cutoff applied to both children (RCUT2BMAX = rcut).
static Context make_d2_mjoin_ctx(int n1, int n2, double rcut = 10.0)
{
    Context ctx;
    // RCUT2BMAX is computed automatically from RCUT2B by Context::add().
    ctx.add("RCTYPE2B", "Cut_Cos");
    ctx.add("RCTYPE2B", "Cut_Cos");
    ctx.add("RCUT2B", rcut);
    ctx.add("RCUT2B", rcut);
    ctx.remove("INIT2B");
    ctx.add("INIT2B", "true");
    ctx.add("TYPE2B", "D2_mJoin");
    ctx.add("TYPE2B", "D2_Blip");
    ctx.add("TYPE2B", n1); ctx.add("TYPE2B", n1);
    ctx.add("TYPE2B", "H"); ctx.add("TYPE2B", "H");
    ctx.add("TYPE2B", "D2_Blip");
    ctx.add("TYPE2B", n2); ctx.add("TYPE2B", n2);
    ctx.add("TYPE2B", "H"); ctx.add("TYPE2B", "H");
    for (int i = 0; i < n1 + n2; ++i) ctx.add("CGRID2B", 0.5 * (i + 1));
    for (int i = 0; i < n1 + n2; ++i) ctx.add("SGRID2B", 0.5);
    return ctx;
}


TEST_CASE_METHOD(ReinitFixture,
    "DM_mJoin::reinit calc_aed result matches fresh construction", "[DM][reinit]")
{
    auto ctx1 = make_dm_mjoin_ctx(0, 2, 1, 3);
    DM_mJoin dm(ctx1);
    dm.set_fidx(0);

    auto ctx2 = make_dm_mjoin_ctx(0, 2, 1, 3);  // same grid as reinit does not change grid dims
    dm.reinit(ctx2);
    // dm.set_fidx(0);

    // Fresh DM_mJoin built directly with ctx2.
    DM_mJoin dm_fresh(ctx2);
    dm_fresh.set_fidx(0);

    REQUIRE(dm.size() == dm_fresh.size());

    // Exercise calc_rho + calc_aed for both; results must agree.
    // Using a single pair of H atoms at known positions.
    Vec3d del(2.1, 3.3, 4.4);
    const double rsq = del * del;
    const double r   = std::sqrt(rsq);

    auto eval = [&](DM_Base& d) {
        rho_type rho;
        d.init_rhoi(rho);
        d.calc_rho(1, 1, r, rsq, del, rho);
        aed_type aed(d.size());
        aed.set_zero();
        d.calc_aed(rho, aed);
        return std::vector<double>(aed.data(), aed.data() + d.size());
    };

    const auto v_reinit = eval(dm);
    const auto v_fresh  = eval(dm_fresh);
    REQUIRE_THAT(v_reinit, Catch::Matchers::Approx(v_fresh).margin(1e-12));
}

TEST_CASE_METHOD(ReinitFixture,
    "DM_mJoin::reinit returns false and preserves values when sizes unchanged",
    "[DM][reinit]")
{
    auto ctx1 = make_dm_mjoin_ctx(1, 3, 0, 2);
    DM_mJoin dm(ctx1);
    dm.set_fidx(0);

    const size_t s_before   = dm.size();
    const size_t rho_before = dm.rhoi_size();

    // Same grid counts, same context → no structural change.
    auto ctx2 = make_dm_mjoin_ctx(1, 3, 0, 2);
    const bool changed = dm.reinit(ctx2);

    REQUIRE(changed == false);
    REQUIRE(dm.size()      == s_before);
    REQUIRE(dm.rhoi_size() == rho_before);
}

/* ------------------------------------------------------------------
 * D2_mJoin: aggregate size is updated after children change
 * ------------------------------------------------------------------ */

TEST_CASE_METHOD(ReinitFixture,
    "D2_mJoin::reinit updates aggregate size", "[D2][reinit]")
{
    // child1: 2 grids → s=2 | child2: 3 grids → s=3  | total: 5
    auto ctx1 = make_d2_mjoin_ctx(2, 3);
    D2_mJoin dm(ctx1);
    dm.set_fidx(0);

    REQUIRE(dm.size() == 2u + 3u);

    // Reinit: child1 gains a grid point (2 → 3).
    auto ctx2 = make_d2_mjoin_ctx(3, 3);
    const bool changed = dm.reinit(ctx2);
    dm.set_fidx(0);

    REQUIRE(changed == true);
    REQUIRE(dm.size() == 3u + 3u);
}

TEST_CASE_METHOD(ReinitFixture,
    "D2_mJoin::reinit calc_aed result matches fresh construction", "[D2][reinit]")
{
    auto ctx1 = make_d2_mjoin_ctx(2, 3);
    D2_mJoin dm(ctx1);
    dm.set_fidx(0);

    auto ctx2 = make_d2_mjoin_ctx(3, 3);
    dm.reinit(ctx2);
    dm.set_fidx(0);

    D2_mJoin dm_fresh(ctx2);
    dm_fresh.set_fidx(0);

    REQUIRE(dm.size() == dm_fresh.size());

    const double r   = 2.5;
    const double rsq = r * r;

    auto eval_d2 = [&](D2_Base& d) {
        aed_type aed(d.size());
        aed.set_zero();
        d.calc_aed(1, 1, r, rsq, aed);
        return std::vector<double>(aed.data(), aed.data() + d.size());
    };

    const auto v_reinit = eval_d2(dm);
    const auto v_fresh  = eval_d2(dm_fresh);
    REQUIRE_THAT(v_reinit, Catch::Matchers::Approx(v_fresh).margin(1e-12));
}

/* ------------------------------------------------------------------
 * Bug 1 regression: outer fcut must be updated when RCUT changes
 *
 * DescriptorsCalc::calc() uses d2.get_rcut() (the D2_mJoin container's
 * own fcut, not any child fcut) as the 2B neighbour pre-filter.
 * Before the fix, reinit() left this outer fcut frozen at the value
 * set during the first Trainer construction.  When RCUT2BMAX grew during
 * HPO, pairs between the old and new rcut were silently excluded from
 * training/validation but included by tadah-predict, causing the saved
 * weights to mismatch the saved descriptor parameters.
 * ------------------------------------------------------------------ */

TEST_CASE_METHOD(ReinitFixture,
    "D2_mJoin::reinit updates outer fcut when RCUT2B changes", "[D2][reinit]")
{
    // Initial context: RCUT2BMAX = 3.0
    auto ctx1 = make_d2_mjoin_ctx(3, 3, 3.0);
    D2_mJoin dm(ctx1);

    REQUIRE(dm.get_rcut() == Approx(3.0));

    // Reinit with a larger cutoff: RCUT2BMAX = 7.0
    auto ctx2 = make_d2_mjoin_ctx(3, 3, 7.0);
    dm.reinit(ctx2);

    // The outer container fcut must reflect the new value.
    // This is what DescriptorsCalc reads as the neighbour pre-filter.
    REQUIRE(dm.get_rcut() == Approx(7.0));
}

TEST_CASE_METHOD(ReinitFixture,
    "D2_mJoin::reinit calc_aed at extended range matches fresh construction",
    "[D2][reinit]")
{
    // r = 5.0 lies between old rcut=3.0 and new rcut=8.0.
    // After reinit the descriptor must equal a freshly-constructed object.
    auto ctx1 = make_d2_mjoin_ctx(3, 3, 3.0);
    D2_mJoin dm(ctx1);
    dm.set_fidx(0);

    auto ctx2 = make_d2_mjoin_ctx(3, 3, 8.0);
    dm.reinit(ctx2);
    dm.set_fidx(0);

    D2_mJoin dm_fresh(ctx2);
    dm_fresh.set_fidx(0);

    const double r   = 5.0;
    const double rsq = r * r;

    auto eval = [&](D2_Base& d) {
        aed_type aed(d.size());
        aed.set_zero();
        d.calc_aed(1, 1, r, rsq, aed);
        return std::vector<double>(aed.data(), aed.data() + d.size());
    };

    const auto v_fresh = eval(dm_fresh);

    // Sanity: r=5.0 is inside rcut=8.0, so values must be non-zero.
    bool any_nonzero = false;
    for (auto val : v_fresh) if (std::abs(val) > 1e-15) { any_nonzero = true; break; }
    REQUIRE(any_nonzero);

    const auto v_reinit = eval(dm);
    REQUIRE_THAT(v_reinit, Catch::Matchers::Approx(v_fresh).margin(1e-12));
}

TEST_CASE_METHOD(ReinitFixture,
    "DM_mJoin::reinit updates outer fcut when RCUTMB changes", "[DM][reinit]")
{
    // Initial context: RCUTMBMAX = 3.0
    auto ctx1 = make_dm_mjoin_ctx(0, 3, 1, 3, 3.0);
    DM_mJoin dm(ctx1);

    REQUIRE(dm.get_rcut() == Approx(3.0));

    // Reinit with a larger cutoff: RCUTMBMAX = 7.0
    auto ctx2 = make_dm_mjoin_ctx(0, 3, 1, 3, 7.0);
    dm.reinit(ctx2);

    // The outer container fcut must reflect the new value.
    // This is what DescriptorsCalc reads as the MB neighbour pre-filter.
    REQUIRE(dm.get_rcut() == Approx(7.0));
}

TEST_CASE_METHOD(ReinitFixture,
    "DM_mJoin::reinit calc_rho at extended range matches fresh construction",
    "[DM][reinit]")
{
    // del = (2.1, 3.3, 4.4)  =>  r ~= 5.89, between old rcut=3.0 and new rcut=8.0.
    auto ctx1 = make_dm_mjoin_ctx(0, 3, 1, 3, 3.0);
    DM_mJoin dm(ctx1);
    dm.set_fidx(0);

    auto ctx2 = make_dm_mjoin_ctx(0, 3, 1, 3, 8.0);
    dm.reinit(ctx2);
    dm.set_fidx(0);

    DM_mJoin dm_fresh(ctx2);
    dm_fresh.set_fidx(0);

    Vec3d del(2.1, 3.3, 4.4);
    const double rsq = del * del;
    const double r   = std::sqrt(rsq);

    auto eval_rho = [&](DM_Base& d) {
        rho_type rho;
        d.init_rhoi(rho);
        d.calc_rho(1, 1, r, rsq, del, rho);
        return std::vector<double>(rho.data(), rho.data() + rho.size());
    };

    const auto v_fresh = eval_rho(dm_fresh);

    // Sanity: r ~= 5.89 is inside rcut=8.0, so density values must be non-zero.
    bool any_nonzero = false;
    for (auto val : v_fresh) if (std::abs(val) > 1e-15) { any_nonzero = true; break; }
    REQUIRE(any_nonzero);

    const auto v_reinit = eval_rho(dm);
    REQUIRE_THAT(v_reinit, Catch::Matchers::Approx(v_fresh).margin(1e-12));
}

/* ------------------------------------------------------------------
 * DC_Selector::reinit: the many-body block follows the two-body block
 * ------------------------------------------------------------------ */

// A DC_Selector context: D2_Blip of n2b functions, then DM_Blip (Lmax 0,
// 2 functions), H only, both cutoffs 5.0 (Cut_Cos), no BIAS.
static Context make_selector_ctx(int n2b)
{
    Context ctx;
    ctx.add("RCTYPE2B", "Cut_Cos");
    ctx.add("RCUT2B", 5.0);
    ctx.remove("INIT2B");
    ctx.add("INIT2B", "true");
    ctx.add("TYPE2B", "D2_Blip");
    ctx.add("TYPE2B", n2b); ctx.add("TYPE2B", n2b);
    ctx.add("TYPE2B", "H"); ctx.add("TYPE2B", "H");
    for (int i = 0; i < n2b; ++i) ctx.add("CGRID2B", 0.5 * (i + 1));
    for (int i = 0; i < n2b; ++i) ctx.add("SGRID2B", 1.0);
    ctx.add("RCTYPEMB", "Cut_Cos");
    ctx.add("RCUTMB", 5.0);
    ctx.remove("INITMB");
    ctx.add("INITMB", "true");
    ctx.add("TYPEMB", "DM_Blip");
    ctx.add("TYPEMB", 0); ctx.add("TYPEMB", 2); ctx.add("TYPEMB", 2);
    ctx.add("TYPEMB", "H"); ctx.add("TYPEMB", "H");
    ctx.add("CGRIDMB", 1.0); ctx.add("CGRIDMB", 2.0);
    ctx.add("SGRIDMB", 1.0); ctx.add("SGRIDMB", 1.0);
    return ctx;
}

/* Pins: DC_Selector::reinit that changes the two-body size (D2_Blip of 2
 * functions to 3) moves the many-body block after it, as init() places it:
 * dmb's first index is 3 (bias 0 + the 3 two-body values), dsize is 3 + the
 * many-body size, and the descriptor of an H pair (r = 2.0) that the
 * selector's descriptors write is the one a DC_Selector built from the new
 * context writes, value for value.
 * Before the fix: reinit updated dsize but left dmb's first index at 2 (the
 * old two-body size), so dmb wrote its 2 values from the third value on:
 * the first was added onto the third two-body value, and the last value
 * stayed 0 ({..., 0.4972, 0.4284, 0} for {..., 0.4704, 0.0268, 0.4284}); a
 * reinit that shrank the two-body block made dmb write past dsize. */
TEST_CASE_METHOD(ReinitFixture,
    "DC_Selector::reinit moves the many-body block after a two-body size "
    "change", "[DC_Selector][reinit]")
{
    DC_Selector dc(make_selector_ctx(2));
    REQUIRE(dc.d2b->size() == 2);
    REQUIRE(dc.dmb->get_fidx() == 2);
    const size_t mb = dc.dmb->size();
    REQUIRE(mb > 0);

    const Context ctx3 = make_selector_ctx(3);
    CHECK(dc.reinit(ctx3));
    CHECK(dc.d2b->size() == 3);
    CHECK(dc.dsize == 3 + mb);
    CHECK(dc.dmb->get_fidx() == 3);

    DC_Selector fresh(ctx3);
    REQUIRE(fresh.dsize == 3 + mb);
    REQUIRE(fresh.dmb->get_fidx() == 3);

    // The descriptor of one H-H pair, as the selector's descriptors write it.
    const Vec3d del(2.0, 0.0, 0.0);
    const double rsq = del * del;
    const double r = std::sqrt(rsq);
    auto descriptor = [&](DC_Selector &s) {
        aed_type aed(s.dsize);
        aed.set_zero();
        s.d2b->calc_aed(1, 1, r, rsq, aed);
        rho_type rho;
        s.dmb->init_rhoi(rho);
        s.dmb->calc_rho(1, 1, r, rsq, del, rho);
        s.dmb->calc_aed(rho, aed);
        return std::vector<double>(aed.data(), aed.data() + s.dsize);
    };
    const std::vector<double> expected = descriptor(fresh);
    REQUIRE(expected.back() != 0.0);
    CHECK(descriptor(dc) == expected);
}

// What a selector's descriptors write for one H-H pair (r = 2.0 along x),
// into buffers longer than any of the sizes compared (n descriptor values,
// nr density values), so that a value written past the selector's sizes is
// compared too: the descriptor (aed), the density array after calc_aed
// (rho_i, then F'(rho_i)) and the force descriptor (fd, x direction).
struct SelectorValues {
    std::vector<double> aed, rho, fd;
    std::size_t init_rhoi_size = 0;  // what init_rhoi gives the array
};
static SelectorValues selector_values(DC_Selector &s, std::size_t n,
                                      std::size_t nr)
{
    const Vec3d del(2.0, 0.0, 0.0);
    const double rsq = del * del;
    const double r = std::sqrt(rsq);
    SelectorValues v;
    rho_type sized;
    s.dmb->init_rhoi(sized);
    v.init_rhoi_size = sized.size();
    aed_type aed(n);
    aed.set_zero();
    fd_type fd(n);
    s.d2b->calc_all(1, 1, r, rsq, aed, fd);
    rho_type rho(nr);
    rho.set_zero();
    s.dmb->calc_rho(1, 1, r, rsq, del, rho);
    s.dmb->calc_aed(rho, aed);
    s.dmb->calc_dXijdri(1, 1, r, rsq, del, rho, fd);
    v.aed.assign(aed.data(), aed.data() + n);
    v.rho.assign(rho.data(), rho.data() + nr);
    for (std::size_t i = 0; i < n; ++i) v.fd.push_back(fd(i, 0));
    return v;
}

// A DC_Selector context: D2_Blip of 2 functions, then DM_mJoin of two
// DM_Blip (Lmax 0, n1 and n2 functions), H only, Cut_Cos at 5.0, no BIAS.
static Context make_selector_dm_mjoin_ctx(int n1, int n2)
{
    Context ctx;
    ctx.add("RCTYPE2B", "Cut_Cos");
    ctx.add("RCUT2B", 5.0);
    ctx.remove("INIT2B");
    ctx.add("INIT2B", "true");
    ctx.add("TYPE2B", "D2_Blip");
    ctx.add("TYPE2B", 2); ctx.add("TYPE2B", 2);
    ctx.add("TYPE2B", "H"); ctx.add("TYPE2B", "H");
    ctx.add("CGRID2B", 0.5); ctx.add("CGRID2B", 1.0);
    ctx.add("SGRID2B", 1.0); ctx.add("SGRID2B", 1.0);
    ctx.add("RCTYPEMB", "Cut_Cos"); ctx.add("RCTYPEMB", "Cut_Cos");
    ctx.add("RCUTMB", 5.0); ctx.add("RCUTMB", 5.0);
    ctx.remove("INITMB");
    ctx.add("INITMB", "true");
    ctx.add("TYPEMB", "DM_mJoin");
    for (int n : {n1, n2}) {
        ctx.add("TYPEMB", "DM_Blip");
        ctx.add("TYPEMB", 0); ctx.add("TYPEMB", n); ctx.add("TYPEMB", n);
        ctx.add("TYPEMB", "H"); ctx.add("TYPEMB", "H");
    }
    for (int i = 0; i < n1 + n2; ++i) ctx.add("CGRIDMB", 0.5 * (i + 1));
    for (int i = 0; i < n1 + n2; ++i) ctx.add("SGRIDMB", 1.0);
    return ctx;
}

/* Pins: DC_Selector::reinit that grows a child of a DM_mJoin many-body
 * descriptor (the first DM_Blip from 2 functions to 3) leaves the selector
 * as a DC_Selector built from the new context: dsize 7 (2 + 5) and rsize
 * 10, and the descriptor, the density array and the force descriptor of an
 * H pair are the fresh selector's, value for value, also past dsize and
 * rsize (none written there); init_rhoi gives the array rsize values.
 * Before the fix: DM_mJoin::reinit kept its total size (4), its density
 * sizes and the density offsets of its children, so dsize stayed 6, rsize
 * and init_rhoi 8, the second child's density overlapped the first's, and
 * the descriptor was {0.0205, 0.1636, 0.0004, 0.0268, 1.2655, 1.2655,
 * 0.2213, 0, 0} against the fresh {0.0205, 0.1636, 0.0004, 0.0268, 0.2213,
 * 0.4284, 0.2213, 0, 0}: its value 6 written past dsize 6; the density
 * array and the force descriptor differed as well. */
TEST_CASE_METHOD(ReinitFixture,
    "DC_Selector::reinit that grows a DM_mJoin child gives the selector a "
    "fresh one has", "[DC_Selector][reinit][DM]")
{
    DC_Selector dc(make_selector_dm_mjoin_ctx(2, 2));
    REQUIRE(dc.dsize == 6);
    REQUIRE(dc.rsize == 8);

    const Context ctx3 = make_selector_dm_mjoin_ctx(3, 2);
    CHECK(dc.reinit(ctx3));
    DC_Selector fresh(ctx3);
    REQUIRE(fresh.dsize == 7);
    REQUIRE(fresh.rsize == 10);
    CHECK(dc.dsize == fresh.dsize);
    CHECK(dc.rsize == fresh.rsize);

    const std::size_t n = fresh.dsize + 2, nr = fresh.rsize + 2;
    const SelectorValues expected = selector_values(fresh, n, nr);
    REQUIRE(expected.aed[fresh.dsize - 1] != 0.0);
    REQUIRE(expected.init_rhoi_size == fresh.rsize);
    const SelectorValues got = selector_values(dc, n, nr);
    CHECK(got.init_rhoi_size == expected.init_rhoi_size);
    CHECK(got.aed == expected.aed);
    CHECK(got.rho == expected.rho);
    CHECK(got.fd == expected.fd);
}

// A DC_Selector context: BIAS true, ATOM H, D2_mJoin of two D2_Blip (n1
// and n2 functions), Cut_Cos at 5.0, then DM_Blip (Lmax 0, 2 functions).
static Context make_selector_d2_mjoin_ctx(int n1, int n2)
{
    Context ctx;
    ctx.remove("BIAS");
    ctx.add("BIAS", "true");
    ctx.add("ATOM", "H");
    ctx.add("RCTYPE2B", "Cut_Cos"); ctx.add("RCTYPE2B", "Cut_Cos");
    ctx.add("RCUT2B", 5.0); ctx.add("RCUT2B", 5.0);
    ctx.remove("INIT2B");
    ctx.add("INIT2B", "true");
    ctx.add("TYPE2B", "D2_mJoin");
    for (int n : {n1, n2}) {
        ctx.add("TYPE2B", "D2_Blip");
        ctx.add("TYPE2B", n); ctx.add("TYPE2B", n);
        ctx.add("TYPE2B", "H"); ctx.add("TYPE2B", "H");
    }
    for (int i = 0; i < n1 + n2; ++i) ctx.add("CGRID2B", 0.5 * (i + 1));
    for (int i = 0; i < n1 + n2; ++i) ctx.add("SGRID2B", 1.0);
    ctx.add("RCTYPEMB", "Cut_Cos");
    ctx.add("RCUTMB", 5.0);
    ctx.remove("INITMB");
    ctx.add("INITMB", "true");
    ctx.add("TYPEMB", "DM_Blip");
    ctx.add("TYPEMB", 0); ctx.add("TYPEMB", 2); ctx.add("TYPEMB", 2);
    ctx.add("TYPEMB", "H"); ctx.add("TYPEMB", "H");
    ctx.add("CGRIDMB", 1.0); ctx.add("CGRIDMB", 2.0);
    ctx.add("SGRIDMB", 1.0); ctx.add("SGRIDMB", 1.0);
    return ctx;
}

/* Pins: with BIAS true (ATOM H: one intercept column, 0) a D2_mJoin
 * two-body descriptor starts at column 1: its first index (get_fidx) is 1
 * in a new selector and after a reinit that grows its first D2_Blip from 2
 * functions to 3, and the descriptor, the density array and the force
 * descriptor of an H pair after that reinit are those of a DC_Selector
 * built from the new context, value for value (column 0, the intercept,
 * left 0 by the descriptors).
 * Before the fix: D2_mJoin::set_fidx placed its children but kept its own
 * first index 0, so get_fidx() was 0 in both selectors, and its reinit
 * placed the children again from 0: the two-body values moved onto the
 * intercept column and the last one stayed 0 ({0.0205, 0.1636, 0.4704,
 * 0.6545, 0.4704, 0, 0.0268, 0.4284, 0, 0} against {0, 0.0205, 0.1636,
 * 0.4704, 0.6545, 0.4704, 0.0268, 0.4284, 0, 0}), in the force descriptor
 * too. */
TEST_CASE_METHOD(ReinitFixture,
    "DC_Selector::reinit with BIAS true keeps a D2_mJoin after the "
    "intercept column", "[DC_Selector][reinit][D2]")
{
    DC_Selector dc(make_selector_d2_mjoin_ctx(2, 2));
    REQUIRE(dc.bias == 1);
    CHECK(dc.d2b->get_fidx() == 1);
    CHECK(dc.dmb->get_fidx() == 5);

    const Context ctx3 = make_selector_d2_mjoin_ctx(3, 2);
    CHECK(dc.reinit(ctx3));
    DC_Selector fresh(ctx3);
    REQUIRE(fresh.dsize == 1 + 5 + 2);
    CHECK(fresh.d2b->get_fidx() == 1);
    CHECK(dc.dsize == fresh.dsize);
    CHECK(dc.d2b->get_fidx() == 1);
    CHECK(dc.dmb->get_fidx() == 6);

    const std::size_t n = fresh.dsize + 2, nr = fresh.rsize + 2;
    const SelectorValues expected = selector_values(fresh, n, nr);
    REQUIRE(expected.aed[0] == 0.0);
    REQUIRE(expected.aed[1] != 0.0);
    const SelectorValues got = selector_values(dc, n, nr);
    CHECK(got.aed == expected.aed);
    CHECK(got.rho == expected.rho);
    CHECK(got.fd == expected.fd);
}

// A DC_Selector context: D2_Blip of 2 functions, then DM_Blip (Lmax lmax,
// n functions), H only, Cut_Cos at 5.0, no BIAS.
static Context make_selector_dm_blip_ctx(int lmax, int n)
{
    Context ctx = make_selector_ctx(2);
    for (const char *k : {"TYPEMB", "CGRIDMB", "SGRIDMB"}) ctx.remove(k);
    ctx.add("TYPEMB", "DM_Blip");
    ctx.add("TYPEMB", lmax); ctx.add("TYPEMB", n); ctx.add("TYPEMB", n);
    ctx.add("TYPEMB", "H"); ctx.add("TYPEMB", "H");
    for (int i = 0; i < n; ++i) ctx.add("CGRIDMB", 1.0 + 0.5 * i);
    for (int i = 0; i < n; ++i) ctx.add("SGRIDMB", 1.0);
    return ctx;
}

/* Pins: a reinit to a context whose BIAS changes (ATOM H, false to true,
 * and true to false; the descriptors unchanged) returns true and gives the
 * selector the intercept columns of the new context, as a fresh one: bias
 * 1 and 0, dsize, the first index of the two-body (D2_mJoin) and of the
 * many-body block, and the descriptor of an H pair.
 * Before the fix: reinit kept the bias of the context the selector was
 * built from and returned false (the descriptors' sizes were unchanged):
 * bias 0 for 1 (dsize 6 for 7, the descriptors from column 0 instead of 1)
 * and 1 for 0 (dsize 7 for 6, from column 1 instead of 0). */
TEST_CASE_METHOD(ReinitFixture,
    "DC_Selector::reinit follows a change of BIAS", "[DC_Selector][reinit]")
{
    auto with_bias = [](bool bias) {
        Context ctx = make_selector_d2_mjoin_ctx(2, 2);
        ctx.remove("BIAS");
        ctx.add("BIAS", bias ? "true" : "false");
        return ctx;
    };
    for (bool from : {false, true}) {
        DYNAMIC_SECTION("BIAS " << from << " to " << !from) {
            DC_Selector dc(with_bias(from));
            REQUIRE(dc.bias == (from ? 1 : 0));
            const Context ctx2 = with_bias(!from);
            CHECK(dc.reinit(ctx2));
            DC_Selector fresh(ctx2);
            CHECK(dc.bias == fresh.bias);
            CHECK(dc.dsize == fresh.dsize);
            CHECK(dc.d2b->get_fidx() == static_cast<size_t>(fresh.bias));
            CHECK(dc.dmb->get_fidx() == fresh.dmb->get_fidx());
            const std::size_t n = 1 + 4 + 2 + 2, nr = fresh.rsize + 2;
            const SelectorValues expected = selector_values(fresh, n, nr);
            const SelectorValues got = selector_values(dc, n, nr);
            CHECK(got.aed == expected.aed);
            CHECK(got.fd == expected.fd);
        }
    }
}

/* Pins: a reinit that changes a DM_Blip many-body descriptor of a selector
 * (Lmax 0: 2 functions to 3 and 3 to 2; Lmax 1 with 2 functions to Lmax 0
 * with 4, the same 4 descriptor values from 8 density values to 4) returns
 * true and leaves the selector as a fresh one: dsize, rsize, the first
 * index of F'(rho_i) (get_drfidx), and the density array (rho_i, then
 * F'(rho_i) right after it), the descriptor and the force descriptor of
 * an H pair, value for value, also past dsize and rsize.
 * Before the fix: DM_Blip::reinit updated its density size but not the
 * first index of F'(rho_i), which stayed at the old density size, and the
 * selector updated rsize only when the descriptor size changed. From 2 to 3
 * functions the index was 2 (fresh 3): F'(rho_i) overwrote the third
 * density and the last descriptor value was 0.1071 for 0.4284 (the density
 * array and the force descriptor differed too); from 3 to 2 it was 3
 * (fresh 2): F'(rho_i) was written at 3 and 4, past rsize 4 (the
 * descriptor and the force descriptor agreed); from Lmax 1 to Lmax 0
 * reinit returned false, rsize stayed 16 (fresh 8) and the index 8 (fresh
 * 4): F'(rho_i) was written at 8 to 11, past rsize 8. */
TEST_CASE_METHOD(ReinitFixture,
    "DC_Selector::reinit that resizes a DM_Blip gives the selector a fresh "
    "one has", "[DC_Selector][reinit][DM]")
{
    struct Change { int lmax0, n0, lmax1, n1; };
    for (const Change &c : {Change{0, 2, 0, 3}, Change{0, 3, 0, 2},
                            Change{1, 2, 0, 4}}) {
        DYNAMIC_SECTION("DM_Blip Lmax " << c.lmax0 << ", " << c.n0
                        << " functions to Lmax " << c.lmax1 << ", " << c.n1
                        << " functions") {
            DC_Selector dc(make_selector_dm_blip_ctx(c.lmax0, c.n0));
            const Context ctx2 = make_selector_dm_blip_ctx(c.lmax1, c.n1);
            CHECK(dc.reinit(ctx2));
            DC_Selector fresh(ctx2);
            REQUIRE(fresh.dmb->get_drfidx() == fresh.dmb->rhoi_size());
            CHECK(dc.dsize == fresh.dsize);
            CHECK(dc.rsize == fresh.rsize);
            CHECK(dc.dmb->get_drfidx() == fresh.dmb->get_drfidx());
            const std::size_t n = 2 + 4 + 2, nr = 16 + 2;
            const SelectorValues expected = selector_values(fresh, n, nr);
            const SelectorValues got = selector_values(dc, n, nr);
            CHECK(got.aed == expected.aed);
            CHECK(got.rho == expected.rho);
            CHECK(got.fd == expected.fd);
        }
    }
}

/* Pins: DM_Blip::reinit and DM_EAD::reinit that change the number of
 * functions (Lmax 0, 2 to 3) place F'(rho_i) right after the new density,
 * as the constructor places it: get_drfidx() is rhoi_size() (3), that of a
 * descriptor built from the new context, and the descriptor of an H pair
 * equals the fresh one's.
 * Before the fix: get_drfidx() stayed 2 (the old density size): F'(rho_i)
 * overwrote the third density, and the descriptor's last value differed
 * (DM_Blip 0.1071, DM_EAD 0.2319, for 0.4284). */
TEST_CASE_METHOD(ReinitFixture,
    "DM_Blip and DM_EAD reinit place F'(rho) after the new density",
    "[DM][reinit]")
{
    for (const std::string name : {"DM_Blip", "DM_EAD"}) {
        DYNAMIC_SECTION(name) {
            auto ctx_of = [&](int n) {
                Context ctx = make_selector_dm_blip_ctx(0, n);
                ctx.remove("TYPEMB");
                ctx.add("TYPEMB", name);
                ctx.add("TYPEMB", 0); ctx.add("TYPEMB", n);
                ctx.add("TYPEMB", n);
                ctx.add("TYPEMB", "H"); ctx.add("TYPEMB", "H");
                return ctx;
            };
            std::unique_ptr<DM_Base> dm(
                tadah::core::factory<DM_Base, const Context &>(name,
                                                                ctx_of(2)));
            const Context ctx3 = ctx_of(3);
            CHECK(dm->reinit(ctx3));
            std::unique_ptr<DM_Base> fresh(
                tadah::core::factory<DM_Base, const Context &>(name, ctx3));
            REQUIRE(fresh->rhoi_size() == 3);
            REQUIRE(fresh->get_drfidx() == 3);
            CHECK(dm->rhoi_size() == 3);
            CHECK(dm->get_drfidx() == fresh->get_drfidx());
            auto descriptor = [](DM_Base &d) {
                const Vec3d del(2.0, 0.0, 0.0);
                const double rsq = del * del;
                rho_type rho;
                d.init_rhoi(rho);
                d.calc_rho(1, 1, std::sqrt(rsq), rsq, del, rho);
                aed_type aed(d.size());
                aed.set_zero();
                d.calc_aed(rho, aed);
                return std::vector<double>(aed.data(), aed.data() + d.size());
            };
            const std::vector<double> expected = descriptor(*fresh);
            REQUIRE(expected.back() != 0.0);
            CHECK(descriptor(*dm) == expected);
        }
    }
}

/* ------------------------------------------------------------------
 * Bug 2 regression: M_Core::reset_trained() must clear the trained flag
 *
 * M_Core::trained is set to true after the first model::train() call.
 * A model trained again in place must be reset before training, or the
 * training call throws "This object is already trained!".  tadah hpo did
 * so through Trainer::reinit() (MLIP, removed: every evaluation now builds
 * its own Trainer); reset_trained() stays the way to retrain a model, and
 * this test ensures the flag is cleared correctly.
 * ------------------------------------------------------------------ */

namespace {

// Minimal concrete M_Core subclass used only by the test below.
struct ConcreteModel : tadah::models::M_Core {
    double predict(const tadah::core::aed_type&) const override { return 0.0; }
};

} // anonymous namespace

TEST_CASE("M_Core::reset_trained clears the trained flag", "[model][reinit]")
{
    ConcreteModel m;
    REQUIRE(m.is_trained() == false);   // default: not yet trained

    m.trained = true;
    REQUIRE(m.is_trained() == true);

    m.reset_trained();
    REQUIRE(m.is_trained() == false);   // regression: must clear the flag
}
