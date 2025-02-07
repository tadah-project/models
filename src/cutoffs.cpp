#include <tadah/models/cutoffs.h>
#include <cmath>
#include <stdexcept>


template<> CONFIG::Registry<Cut_Base,double>::Map CONFIG::Registry<Cut_Base,double>::registry{};
CONFIG::Registry<Cut_Base,double>::Register<Cut_Cos> Cut_Cos_1( "Cut_Cos" );
CONFIG::Registry<Cut_Base,double>::Register<Cut_Cos_S> Cut_Cos_S_1( "Cut_Cos_S" );
CONFIG::Registry<Cut_Base,double>::Register<Cut_Tanh> Cut_Tanh_1( "Cut_Tanh" );
CONFIG::Registry<Cut_Base,double>::Register<Cut_Poly2> Cut_Poly_1( "Cut_Poly2" );
CONFIG::Registry<Cut_Base,double>::Register<Cut_Dummy> Cut_Dummy_1( "Cut_Dummy" );
CONFIG::Registry<Cut_Base,double>::Register<Cut_PT> Cut_PT_1( "Cut_PT" );


Cut_Base::~Cut_Base() {}
void Cut_Base::test_rcut(const double r)
{
    if (r<0.0)
        throw std::runtime_error("Cutoff distance cannot be negative.");
}

Cut_Cos::Cut_Cos() {}
Cut_Cos::Cut_Cos(double rcut)
{
    set_rcut(rcut);
    test_rcut(rcut);
}

std::string Cut_Cos::label() {
    return lab;
}

void Cut_Cos::set_rcut(const double r) {
    test_rcut(r);
    rcut=r;
    rcut_sq=r*r;
    rcut_inv= r<=0 ? 0.0 : 1.0/r;
}
double Cut_Cos::get_rcut() {
    return rcut;
}

double Cut_Cos::get_rcut_sq() {
    return rcut_sq;
}

double Cut_Cos::calc(double r) {
    if (r>=rcut) return 0.0;
    return 0.5 * (std::cos(M_PI*r*rcut_inv) + 1.0);
}
double Cut_Cos::calc_prime(double r) {
    if (r>=rcut) return 0.0;
    double t = M_PI*rcut_inv;
    return -0.5 * t*sin(t*r);

}

Cut_Dummy::Cut_Dummy() {}
Cut_Dummy::Cut_Dummy(double rcut)
{
    set_rcut(rcut),
    test_rcut(rcut);
}

std::string Cut_Dummy::label() {
    return lab;
}

void Cut_Dummy::set_rcut(const double r) {
    test_rcut(r);
    rcut=r;
    rcut_sq=r*r;
    rcut_inv= r<=0 ? 0.0 : 1.0/r;
}
double Cut_Dummy::get_rcut() {
    return rcut;
}

double Cut_Dummy::get_rcut_sq() {
    return rcut_sq;
}

double Cut_Dummy::calc(double r) {
    if (r>=rcut) return 0.0;
    return 1.0;
}
double Cut_Dummy::calc_prime(double) {
    return 0.0;
}

Cut_Tanh::Cut_Tanh() {}
Cut_Tanh::Cut_Tanh(double rcut)
{
    set_rcut(rcut);
    test_rcut(rcut);
}

std::string Cut_Tanh::label() {
    return lab;
}

void Cut_Tanh::set_rcut(const double r) {
    test_rcut(r);
    rcut=r;
    rcut_sq=r*r;
    rcut_inv= r<=0 ? 0.0 : 1.0/r;
}
double Cut_Tanh::get_rcut() {
    return rcut;
}

double Cut_Tanh::get_rcut_sq() {
    return rcut_sq;
}

double Cut_Tanh::calc(double r) {
    if (r>=rcut) return 0.0;
    double t = std::tanh(1.0 - r*rcut_inv);
    return t*t*t;
}
double Cut_Tanh::calc_prime(double r) {
    if (r>=rcut) return 0.0;
    double a = (1.0-r*rcut_inv);
    double t = std::tanh(a);
    double s = 1.0/(std::cosh(a));
    return -3.0*t*t*rcut_inv*s*s;
}

Cut_Poly2::Cut_Poly2() {}
Cut_Poly2::Cut_Poly2(double rcut)
{
    if (rcut<=1.0)
        throw std::runtime_error(
                "Cut_Poly2 cutoff distance cannot be smaller than 1.0.");
    set_rcut(rcut);
    test_rcut(rcut);
}

std::string Cut_Poly2::label() {
    return lab;
}

void Cut_Poly2::set_rcut(const double r) {
    test_rcut(r);
    rcut=r;
    rcut_sq=r*r;
    rcut_inv= r<=0 ? 0.0 : 1.0/r;
    rcut_inner=rcut-1.0;
}
double Cut_Poly2::get_rcut() {
    return rcut;
}

double Cut_Poly2::get_rcut_sq() {
    return rcut_sq;
}

double Cut_Poly2::calc(double r) {
    if (r>=rcut) return 0.0;
    else if (r<= rcut_inner) return 1.0;
    double rs=r-rcut_inner;
    return rs*rs*rs*(rs*(15.0-6.0*rs)-10.0)+1;
}
double Cut_Poly2::calc_prime(double r) {
    if (r>=rcut) return 0.0;
    else if (r<= rcut_inner) return 0.0;
    double rs=r-rcut_inner;
    return -30.0*(rs-1.0)*(rs-1.0)*rs*rs;
}

Cut_Cos_S::Cut_Cos_S() {}
Cut_Cos_S::Cut_Cos_S(double rcut)
{
    set_rcut(rcut);
    test_rcut(rcut);
}

std::string Cut_Cos_S::label() {
    return lab;
}

void Cut_Cos_S::set_rcut(const double r) {
    test_rcut(r);
    rcut=r;
    rcut_sq=r*r;
    rcut_inv= r<=0 ? 0.0 : 1.0/r;
    rcut_inner=rcut-1.0;
}
double Cut_Cos_S::get_rcut() {
    return rcut;
}

double Cut_Cos_S::get_rcut_sq() {
    return rcut_sq;
}

double Cut_Cos_S::calc(double r) {
  if (r>=rcut) return 0.0;
  else if (r<= rcut_inner) return 1.0;
  double rs = (r-rcut_inner)/(rcut-rcut_inner);
  return 0.5*(1+std::cos(M_PI*rs));
}
double Cut_Cos_S::calc_prime(double r) {
  if (r>=rcut || r<= rcut_inner) return 0.0;
  else {
    double rs = (r - rcut_inner) / (rcut - rcut_inner);
    double drs_dr = 1.0 / (rcut - rcut_inner);
    return -0.5 * M_PI * std::sin(M_PI * rs) * drs_dr;
  }
}

Cut_PT::Cut_PT() {}
Cut_PT::Cut_PT(double rcut) {
  set_rcut(rcut);
  test_rcut(rcut);
}
std::string Cut_PT::label() {
  return lab;
}

// Set the outer cutoff radius
void Cut_PT::set_rcut(const double r) {
    test_rcut(r);
    rcut=r;
    rcut_sq=r*r;
    rcut_inv= r<=0 ? 0.0 : 1.0/r;
    rcut_inner=-rcut/2;
}
double Cut_PT::get_rcut() {
  return rcut;
}
double Cut_PT::get_rcut_sq() {
  return rcut_sq;
}
double Cut_PT::calc(double r) {
  if (r>=rcut) return 0.0;
  double u = 0.5 * ((rcut / (r - rcut_inner)) + (rcut / (r - rcut)));
  return 0.5 + 0.5 * tanh(u);
}
double Cut_PT::calc_prime(double r) {
  if (r>=rcut) return 0.0;
  double u = 0.5 * ((rcut / (r - rcut_inner)) + (rcut / (r - rcut)));
  double cosh_u = cosh(u);
  double sech2_u = 1.0 / (cosh_u * cosh_u);

  double term1 = 1.0 / ((r - rcut_inner) * (r - rcut_inner));
  double term2 = 1.0 / ((r - rcut) * (r - rcut));
  double sum_terms = term1 + term2;

  return - (rcut / 4.0) * sech2_u * sum_terms;
}
