#include <cassert>
#include <sstream>
#include <iomanip>
#include <functional>

#include "Vec.hh"
#include "particle.h"
#include "fcache.h"
#include "root_finding.h"
#include "newton_solver.h"
#include "nlopt.hpp"

using namespace std;
using namespace Unit;

double particle::angle = 45 * deg;
particle::HCSFORM particle::hcsform = Jokipii_Thomas;

particle::particle() {
  mass = 0.93827 * GeV;
  B0 = 5 * nT;
  polarity = -1;
  angle = 15 * deg;
  D = 5 * 1e22 * cm * cm / sec;
  indexA = 2;

  theta = pi / 2.0 + 1e-6;
  Vs_eq = Wind();

  Bn = B0 * AU * AU / 1.35883;
}

particle::particle(const map<string, docopt::value>& args)
{
  auto fargs = [&](const string& key) -> double {
    return stod(args.at(key).asString());
  };
  mass = fargs("--mass") * GeV;
  B0 = fargs("--B0") * nT;
  polarity = args.at("--polarity").asLong();
  angle = fargs("--angle") * deg;
  D = fargs("--D") * 1e22 * cm * cm / sec;
  indexA = fargs("--indexA");

  theta = pi / 2.0 + 1e-6;
  Vs_eq = Wind();

  Bn = B0 * AU * AU / 1.35883;
  // Bn = 4.7*1e13;

  //hcsform = Jokipii_Thomas;
//   Ek = mass;
//   std::cout << "Mass:  " << mass/GeV << "   " << Ek/GeV << std::endl;
//   getchar();
}

particle::~particle() {}

double particle::Wind(double r, double theta, double phi, double angle) const{
  double value;
    if (0. <= theta && theta <= pi / 2.) {
      value = (1.475 - 0.4 * tanh(6.8 * ((theta - pi / 2) + (15 * deg + angle)))) * (3.5 / 5. - 1.5 / 5. * tanh((r - 120 * AU) / 1.2 / AU));
    } else if (pi / 2. < theta && theta <= pi) {
      value = (1.475 + 0.4 * tanh(6.8 * ((theta - pi / 2) - (15 * deg + angle)))) * (3.5 / 5. - 1.5 / 5. * tanh((r - 120 * AU) / 1.2 / AU));
    }
    return  value * 400 * (km / sec);
}

double particle::Wind() const {
  return Wind(r, theta, phi, angle);
}

void particle::r_bound(double r, double phi, double phi0, double& rlow, double& rup) const {
  double Tr = 2 * pi * Vs_eq / Omega;

  rlow = (phi0 - phi + Omega * (t - t0)) * Vs_eq / Omega;
  rlow += floor((r - rlow) / Tr) * Tr;
  rup = rlow + Tr;
}

double particle::Theta_S_Jokipii_Thomas(double phi0) const {
  static fcache theta_jokipii_thomas([](double x) { return pi / 2 - asin(sin(angle) * sin(x)); }, 1000000);

  return theta_jokipii_thomas(phi0);
}
double particle::Phi0_S_Jokipii_Thomas(double theta) const {
  double ratio = sin(pi / 2 - theta) / sin(angle);
  if (ratio >= 1) return pi / 2;
  else if (ratio <= -1) return -pi / 2;

  return asin(ratio);
}

double particle::Theta_S_Kota_Jokipii(double phi0) const {
  static fcache theta_kota_jokipii([](double x) { return pi / 2 - atan(tan(angle) * sin(x)); }, 1000000);

  return theta_kota_jokipii(phi0);
}
double particle::Phi0_S_Kota_Jokipii(double theta) const {
  double ratio = tan(pi / 2 - theta) / tan(angle);
  if (ratio >= 1) return pi / 2;
  else if (ratio <= -1) return -pi / 2;

  return asin(ratio);
}

double particle::Theta_S(double r, double phi) const {
  if (hcsform == Jokipii_Thomas)
    return Theta_S_Jokipii_Thomas(r, phi);
  else if (hcsform == Kota_Jokipii)
    return Theta_S_Kota_Jokipii(r, phi);

  assert(false && "hcsform not supported");
  return 0;
}
extern "C" double Theta_S_C(double r, double phi) {
  return particle().Theta_S(r * AU, phi);
}
double particle::Phi0_S(double theta) const {
  if (hcsform == Jokipii_Thomas)
    return Phi0_S_Jokipii_Thomas(theta);
  else if (hcsform == Kota_Jokipii)
    return Phi0_S_Kota_Jokipii(theta);

  assert(false && "hcsform not supported");
  return 0;
}

const double particle::Heav() {
  double theta_s = Theta_S(r, phi);
  double value;
  if (theta <= theta_s)
    value = 1.;
  else if (theta_s < theta)
    value = -1.;

  return value;
}

const double particle::B_r(const double& r, const double& theta, const double& phi, const double& heaviside) {
  // double r0 = 1 * AU;
  // std::cout << "D:   " << heaviside << "   " << polarity << "   " << r << "
  // " << r0 << std::endl;
  //     getchar();

  return Bn * heaviside * polarity * pow(1. / r, 2.);
}

const double particle::B_p(const double& r, const double& theta, const double& phi, const double& heaviside) {
  // std::cout << "D:   " << Omega << "   " << sin(theta) << "   " << Vs << " "
  // << polarity<< std::endl; getchar();
  // double r0 = 1 * AU;

  //  std::cout << "D:   " << "  " << Omega * sin(theta) * heaviside * polarity / Vs * r  << std::endl; getchar();
  // std::cout << Bn << "  " << r << "  " << sin(theta ) << "  " << a << std::endl;
  // getchar();
  return -1. * Bn / r * Omega * sin(theta) * heaviside * polarity / Vs;
}

const double particle::K_rr(const double& r, const double& theta, const double& phi, const double& psi, const double& B, const double& B0, const double& M_p, const double& M_p0, const double& V_p) {
  double lamda0 = 0.05 * AU;
  double lamda_p = lamda0 * B0 / B * M_p / M_p0;
  double kx = 1. / 3. * lamda_p * V_p;
  double ky = 0.02 * kx;

  double kr = kx * pow(cos(psi), 2.) + ky * pow(sin(psi), 2.);

  return kr;
}

const double particle::K_tt(const double& r, const double& theta, const double& phi, const double& psi, const double& B, const double& B0, const double& M_p, const double& M_p0, const double& V_p) {
  double lamda0 = 0.05 * AU;
  double lamda_p = lamda0 * B0 / B * M_p / M_p0;
  double kx = 1. / 3. * lamda_p * V_p;
  double ky = 0.02 * kx;

  double kz = ky;
  // if (theta < pi / 2. && 0 <= theta)
  //   kz = ky * (2. - 1. * tanh(8 * ((theta + (-90. + 35.) * deg))));
  // else if (pi / 2. <= theta && theta <= pi)
  //   kz = ky * (2. + 1. * tanh(8 * ((theta + (-90. - 35.) * deg))));

  return kz;
}

const double particle::K_pp(const double& r, const double& theta, const double& phi, const double& psi, const double& B, const double& B0, const double& M_p, const double& M_p0, const double& V_p) {
  double lamda0 = 0.05 * AU;
  double lamda_p = lamda0 * B0 / B * M_p / M_p0;
  double kx = 1. / 3. * lamda_p * V_p;
  double ky = 0.02 * kx;

  double kp = -1. * kx * pow(sin(psi), 2.) + ky * pow(cos(psi), 2.);

  return kp;
}

void particle::HCS_rphi(const double &r, const double &phi, double& x, double& y, double& z) const {
  double cs0 = Theta_S(r, phi);
  x = r * sin(cs0) * cos(phi);
  y = r * sin(cs0) * sin(phi);
  z = r * cos(cs0);
}

struct nlopt_info {
  double x, y, z;
  double r, phi, theta;
  const particle *p;
};
static int iter_d = 0;
double distance_to_point(const std::vector<double> &x, std::vector<double> &grad, void *voidp) {
  iter_d++;
  nlopt_info *info = reinterpret_cast<nlopt_info*>(voidp);

  auto distance = [&](double r, double phi) {
    // r = 53.01 * AU;
    // phi = 1e-10;
    // std::cout << "Find:  " << r / AU << "  " << phi ;
    // getchar();
    double vx, vy, vz;
    info->p->HCS_rphi(r, phi, vx, vy, vz);
    return sqrt((vx - info->x) * (vx - info->x) + (vy - info->y) * (vy - info->y) + (vz - info->z) * (vz - info->z));
  };

  double d = distance(x[0], x[1]);
  double ddr = distance(x[0] + 1e-10 * AU, x[1]);
  double ddp = distance(x[0], x[1] + 1e-10);

  grad = { (ddr - d) / (1e-10 * AU), (ddp - d) / (1e-10) };

  return d;
}

double particle::get_HCS_distance_old(double ftol_abs) const {
  nlopt_info info = { r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta), r, phi, theta, this };
  nlopt::opt opt(nlopt::LD_MMA, 2);
  opt.set_min_objective(distance_to_point, &info);
//   opt.set_xtol_rel(1e-8);
  opt.set_ftol_abs(ftol_abs);

  double dr = pi / Omega * Vs_eq / 2;
  vector<double> position_hcs, position_hcs_low, position_hcs_up;
  auto assign_region = [&](double v) {
    position_hcs = { fmax(1e-5, r + v), phi };
    position_hcs_low = { fmax(0, position_hcs[0] - dr), position_hcs[1] - pi };
    position_hcs_up = { position_hcs[0] + dr, position_hcs[1] + pi };
    opt.set_lower_bounds(position_hcs_low);
    opt.set_upper_bounds(position_hcs_up);
  };

  double dlow, dmid, dup;
  assign_region(-dr);
  opt.optimize(position_hcs, dlow);

  //cout << "pbest: " << position_hcs[0] / AU << " " << Theta_S(position_hcs[0], position_hcs[1]) / deg << " " << position_hcs[1] / deg << endl;

  assign_region(0);
  opt.optimize(position_hcs, dmid);
  //cout << "pbest: " << position_hcs[0] / AU << " " << Theta_S(position_hcs[0], position_hcs[1]) / deg << " " << position_hcs[1] / deg << endl;

  assign_region(dr);
  opt.optimize(position_hcs, dup);
  //cout << "pbest: " << position_hcs[0] / AU << " " << Theta_S(position_hcs[0], position_hcs[1]) / deg << " " << position_hcs[1] / deg << endl;

  double vx, vy, vz;
  HCS_rphi(r, phi, vx, vy, vz);
  double sign = info.z > vz ? 1 : -1;

  //cout << "dlow: " << dlow / AU << " dmid: " << dmid / AU << " dup: " << dup / AU << endl;

  return sign * fmin(fmin(dlow, dup), dmid);
}

static int iter_d_fix = 0;
double distance_to_point_fix_phi(const std::vector<double> &x, std::vector<double> &grad, void *voidp) {
  nlopt_info *info = reinterpret_cast<nlopt_info*>(voidp);

  auto distance = [&](double r) {
    double vx, vy, vz;
    info->p->HCS_rphi(r, info->phi, vx, vy, vz);
    return sqrt((vx - info->x) * (vx - info->x) + (vy - info->y) * (vy - info->y) + (vz - info->z) * (vz - info->z));
  };


  double d = distance(x[0]);
  double ddr = distance(x[0] + 1e-10 * AU);

  grad = { (ddr - d) / (1e-10 * AU) };
  iter_d_fix++;
  return d;
}
double dtheta_to_point(const std::vector<double> &x, std::vector<double> &grad, void *voidp) {
  nlopt_info *info = reinterpret_cast<nlopt_info*>(voidp);

  double d = fabs(info->p->Theta_S(x[0], info->phi) - info->theta);
  double ddr = fabs(info->p->Theta_S(x[0] + 1e-13 * AU, info->phi) - info->theta);

  grad = { (ddr - d) / (1e-13 * AU) };
  iter_d_fix++;
  return d;
}

class SpiralVdot {
  public:
  const Vec p_cs0, target_point;
  const double r_cs0, phi_cs0;
  double ov, theta_cs, phi0, ctheta, stheta;
  SpiralVdot(double ov_, const Vec& p_cs, const Vec& target_point_) :
    p_cs0(p_cs), target_point(target_point_),
    r_cs0(p_cs.len()), phi_cs0(p_cs.phi()),
    ov(ov_), theta_cs(p_cs.theta()), phi0(phi_cs0 + r_cs0 * ov), ctheta(p_cs.z / r_cs0), stheta(sqrt(1 - ctheta * ctheta)) {}

  Vec tangent_vec(double r, const Vec& point) const {
    double rphi = r * stheta;
    double sphi = point.y / rphi, cphi = point.x / rphi;

    Vec dr(stheta * cphi, stheta * sphi, ctheta);
    Vec dphi(-rphi * sphi, rphi * cphi, 0);

    /************************************************************
     * Delta_phi + Omega / V * Delta_r = 0
     * Delta_phi = - Omega / V * Delta_r
     * dr_phi = Delta_phi * dphi + Delta_r * dr
     *        ~ - Omega / V * dphi + dr
     ************************************************************/
    Vec dv = - ov * dphi + dr;
    dv.normalize();
    return dv;
  }

  double operator()(double r) const {
    Vec p_cs;
    if (r == r_cs0) p_cs = p_cs0;
    else {
      double phi = phi0 - r * ov;
      p_cs.set_spherical(r, theta_cs, phi);
    }

    Vec dl = target_point - p_cs;
    return dl.dot(tangent_vec(r, p_cs)) / dl.len();
  }
};

bool particle::spiral_iterate(const Vec& target_point, Vec& p_cs, double& diter) const {
  double ov = Omega / Vs_eq;

  SpiralVdot vdot(ov, p_cs, target_point);

  double vdot0 = vdot(vdot.r_cs0);
  if (vdot0 == 0) return false;

  double r1;
  double dangle = - asin(vdot0) * fmin(diter / vdot.r_cs0, 1); // the dangle should be smaller when the distance is much small than r_cs

  int id = 0;
  do {
    id++;
    r1 = (vdot.phi0 - vdot.phi_cs0 - id * dangle) / ov;
  } while (vdot(r1) * vdot0 > 0);
  //if (r1 < 0) r1 = 1e-5 * AU;

  double rh = ridders_method(vdot, vdot.r_cs0, r1, 1e-3);

  p_cs.set_spherical(rh, vdot.theta_cs, vdot.phi0 -  rh * ov);
  //cout << "++++ " << rh / AU << " " << vdot(rh) << " | " << (target_point - vdot.p_cs0).len() / AU << " " << (target_point - p_cs).len() / AU << endl;

  double diter_next = (target_point - p_cs).len();
  //cout << "diter -> next: " << diter / AU << " " << diter_next / AU << endl;
  assert(rh > 0 && "the spiral_iterate should not give a negative radius");
  assert(diter_next < diter * (1 + 1e-8) && "the spiral_iterate should decrease the distance to the target point");
  diter = diter_next;
  return true;
}

class WaveVdot {
  public:
  const Vec p_cs0, target_point;
  const double r_cs0, theta_cs0;
  double ov, phi_cs, cphi, sphi;
  const particle* p;

  WaveVdot(double ov_, const Vec& p_cs, const Vec& target_point_, const particle* p_) :
    ov(ov_), p_cs0(p_cs), target_point(target_point_), r_cs0(p_cs.len()), theta_cs0(p_cs.theta()), phi_cs(p_cs.phi()), p(p_)
  {
    double rphi = sqrt(p_cs.x * p_cs.x + p_cs.y * p_cs.y);
    cphi = p_cs.x / rphi;
    sphi = p_cs.y / rphi;
  }

  Vec tangent_vec(double r, const Vec& point) const {
    double ctheta = point.z / r, stheta = sqrt(1 - ctheta * ctheta);
    double rphi = stheta * r;

    Vec dr(stheta * cphi, stheta * sphi, ctheta);
    Vec dtheta(r * ctheta * cphi, r * ctheta * sphi, - r * stheta);

    /************************************************************
     * Jokipii_Thomas:  sin(pi / 2 - theta) = sin(angle) * sin(phi + ov * r)
     * Jokipii_Thomas:  - cos(pi / 2 - theta) * Delta_theta = sin(angle) * cos(phi + ov * r) * Delta_r * ov
     * Jokipii_Thomas:  Delta_theta = - sin(angle) * cos(phi + ov * r) / cos(pi / 2 - theta) * Delta_r * ov
     * 
     * Kota_Jokipii:    tan(pi / 2 - theta) =  tan(angle) * sin(phi + ov * r)
     * Kota_Jokipii:    - cos(pi / 2 - theta)^-2 * Delta_theta = tan(angle) * cos(phi + ov * r) * Delta_r * ov
     * Kota_Jokipii:    Delta_theta = - tan(angle) * cos(phi + ov * r) * cos(pi / 2 - theta)^2 * Delta_r * ov
     * 
     *       dtheta_phi = Delta_theta * dtheta + Delta_r * dr
     * Jokipii_Thomas:  ~ - sin(angle) * cos(phi + ov * r) / cos(pi / 2 - theta) * ov * dtheta + dphi
     * Kota_Jokipii:    ~ - tan(angle) * cos(phi + ov * r) * cos(pi / 2 - theta)^2 * ov * dtheta + dphi
     ************************************************************/

    double phi0 = phi_cs + r * ov;

    Vec dv;
    if (p->hcsform == particle::Jokipii_Thomas) dv = - sin(p->angle) * cos(phi0) / stheta * ov * dtheta + dr;
    else if (p->hcsform == particle::Kota_Jokipii) dv = - tan(p->angle) * cos(phi0) * stheta * stheta * ov * dtheta + dr;
    else assert(false && "Unsuported hcsform");

    dv.normalize();
    return dv;
  }

  double operator()(double r, double *vdot_r = NULL) const {
    Vec p_cs;
    if (r == r_cs0) p_cs = p_cs0;
    else
      p_cs.set_spherical(r, p->Theta_S(r, phi_cs), phi_cs);

    //cout << "-- " << r / AU << " " << (target_point - p_cs).dot(tangent_vec(r, p_cs)) / AU << endl;
    Vec dl = target_point - p_cs;
    Vec vtangent = tangent_vec(r, p_cs);
    //auto show_vec = [](const Vec& v) {
    //  cout << " [" << sqrt(v.x * v.x + v.y * v.y) << "," << v.z << "]" << " " << v.phi() / deg;
    //};
    //cout << "r p_cs dl vtangent vdot: " << r / AU << " | ";
    //show_vec(target_point / AU);
    //show_vec(p_cs / AU);
    //show_vec(dl / AU);
    //show_vec(vtangent);
    //cout << " " << dl.dot(tangent_vec(r, p_cs)) / AU << endl;
    if (vdot_r != NULL) *vdot_r = fabs(p_cs.dot(vtangent) / p_cs.len());

    return dl.dot(vtangent) / dl.len();
  }
};

bool particle::wave_iterate(const Vec& target_point, Vec& p_cs, double& diter) const {
  double ov = Omega / Vs_eq;

  WaveVdot vdot(ov, p_cs, target_point, this);

  double vdot_r;
  double vdot0 = vdot(vdot.r_cs0, &vdot_r);
  if (vdot0 == 0) return false;

  double vdr = vdot0 * (target_point - vdot.p_cs0).len() * vdot_r;
  if (fabs(vdr) > pi / 2 / ov) vdr = pi / 2 / ov * (vdr > 0 ? 1 : -1);

  int ir = 0;
  double r1;
  do {
    ir++;
    r1 = vdot.r_cs0 + ir * vdr;
  } while (vdot(r1) * vdot0 > 0);
  //double vd1 = vdot(vdot.r_cs0);
  //double vd2 = vdot(vdot.r_cs0 + vdr);
  //double vd3 = vdot(vdot.r_cs0 + 2 * vdr);
  //cout << "ir = " << ir << " " << vdr / AU << " " << vdot.r_cs0 / AU << " " << (vdot.r_cs0 + vdr) / AU << " " << r1 / AU << " | " << vd1 << " " << vd2 << " " << vd3 << endl;

  //for (double vr = 0.1 * AU; vr < 5.2 * AU; vr += 0.05 * AU)
  //  cout << "-- " << vr / AU << " " << Theta_S(vr, vdot.phi_cs) / deg << " " << vdot(vr) << endl;
  double rh = ridders_method(vdot, vdot.r_cs0, r1, 1e-3);
  //cout << vdot.r_cs0 / AU << " " << rh / AU << " " << r1 / AU << " | " << vdot0 << " " << vdot(rh) << " " << vdot(r1) << endl;
  //cout << "p_cs before: " << p_cs.len() / AU << " " << p_cs.theta() / deg << " " << Theta_S(p_cs.len(), vdot.phi_cs) / deg << " " << vdot.phi_cs / deg << " " << p_cs.phi() / deg << " " << p_cs / AU << endl;
  //cout << "set with: " << rh / AU << " " << Theta_S(rh, vdot.phi_cs) / deg << " " << vdot.phi_cs / deg << endl;
  //p_cs.set_spherical(rh, Theta_S(rh, vdot.phi_cs), vdot.phi_cs);
  //cout << "set over: " << p_cs.len() / AU << " " << p_cs.theta() / deg << " " << p_cs.phi() / deg << endl;
  //cout << "p_cs after:  " << p_cs.len() / AU << " " << p_cs.theta() / deg << " " << Theta_S(p_cs.len(), vdot.phi_cs) / deg << " " << p_cs / AU << endl;
  //cout << "==== " << rh / AU << " " << vdot(rh) << " | " << (target_point - vdot.p_cs0).len() / AU << " " << (target_point - p_cs).len() / AU << endl;

  double diter_next = (target_point - p_cs).len();
  assert(diter_next <= diter * (1 + 1e-8) && "the wave_iterate should decrease the distance to the target point");
  diter = diter_next;
  return true;
}

Vec particle::norm_vec(const Vec& p_cs) const {
  double r_cs = p_cs.len();
  double phi_cs = p_cs.phi();

  double rphi = sqrt(p_cs.x * p_cs.x + p_cs.y * p_cs.y);
  double ctheta = p_cs.z / r_cs, stheta = sqrt(1 - ctheta * ctheta);
  double sphi = p_cs.y / rphi, cphi = p_cs.x / rphi;
  Vec dr(stheta * cphi, stheta * sphi, ctheta);
  Vec dtheta(r_cs * ctheta * cphi, r_cs * ctheta * sphi, - r_cs * stheta);
  Vec dphi(-r_cs * stheta * sphi, r_cs * stheta * cphi, 0);

  double ov = Omega / Vs_eq;
  /************************************************************
   * Delta_phi + Omega / V * Delta_r = 0
   * Delta_phi = - Omega / V * Delta_r
   * dr_phi = Delta_phi * dphi + Delta_r * dr
   *        ~ - Omega / V * dphi + dr
   ************************************************************/
  Vec dr_phi = - ov * dphi + dr;
  //cout << "phi = (" << dphi << ")" << endl;
  //cout << "r = (" << dr << ")" << endl;
  //cout << "theta = (" << dtheta << ")" << endl;

  /************************************************************
   * Jokipii_Thomas:  sin(pi / 2 - theta) = sin(angle) * sin(phi + ov * r)
   * Jokipii_Thomas:  - cos(pi / 2 - theta) * Delta_theta = sin(angle) * cos(phi + ov * r) * Delta_phi
   * Jokipii_Thomas:  Delta_theta = - sin(angle) * cos(phi + ov * r) / cos(pi / 2 - theta) * Delta_phi
   * 
   * Kota_Jokipii:    tan(pi / 2 - theta) =  tan(angle) * sin(phi + ov * r)
   * Kota_Jokipii:    - cos(pi / 2 - theta)^-2 * Delta_theta = tan(angle) * cos(phi + ov * r) * Delta_phi
   * Kota_Jokipii:    Delta_theta = - tan(angle) * cos(phi + ov * r) * cos(pi / 2 - theta)^2 * Delta_phi
   * 
   *       dtheta_phi = Delta_theta * dtheta + Delta_phi * dphi
   * Jokipii_Thomas:  ~ - sin(angle) * cos(phi + ov * r) / cos(pi / 2 - theta) * dtheta + dphi
   * Kota_Jokipii:    ~ - tan(angle) * cos(phi + ov * r) * cos(pi / 2 - theta)^2 * dtheta + dphi
   ************************************************************/
  Vec dtheta_phi;
  if (hcsform == Jokipii_Thomas) dtheta_phi = - sin(angle) * cos(phi0(r_cs, phi_cs)) / stheta * dtheta + dphi;
  else if (hcsform == Kota_Jokipii) dtheta_phi = - tan(angle) * cos(phi0(r_cs, phi_cs)) * stheta * stheta * dtheta + dphi;
  else assert(false && "Unsuported hcsform");

  Vec dh = dr_phi.cross(dtheta_phi);
  dh.normalize();
  return dh;
}

bool particle::point_iterate(const Vec& target_point, Vec& p_cs, Vec& dh, double& diter) const {
  Vec dl = p_cs - target_point;

  Vec pnext = target_point + dh * dh.dot(dl);
  Vec dv = pnext - p_cs;

  if (dv.len() > 1 * AU) {
    dv *= 1 * AU / dv.len();
    pnext = p_cs + dv;
  }

  Vec dh_next;

  do {
    double r_cs = pnext.len();
    double phi_cs = pnext.phi();
    double theta_cs = Theta_S(r_cs, phi_cs);
    pnext.set_spherical(r_cs, theta_cs, phi_cs);

    if (dv.len() < 1) break;
    dh_next = norm_vec(pnext);

    if ((pnext - target_point).len() < dl.len() && dl.cross(dh).dot((pnext - target_point).cross(dh_next)) > 0) break;

    dv *= 0.5;
    pnext = p_cs + dv;
  } while (true);

  p_cs = pnext;
  dh = dh_next;
  diter = (target_point - p_cs).len();
  return true;
}

double particle::get_HCS_distance_polygon(double Rg2, Polygon polygon) const {
    static map<Polygon, vector<vector<double>>> polygon_coordinates;
    if (polygon_coordinates.size() == 0) {
        polygon_coordinates[Polygon::Dodecahedron] = {
            {-8.49106388e-01,  4.09765508e-01},
            { 0.00000000e+00,  0.00000000e+00},
            { 9.27795731e-01, -1.67583849e-01},
            { 2.06010245e-01,  6.34039244e-01},
            {-6.09029588e-01, -7.19703286e-01},
            { 7.79418897e-01,  5.30466663e-01},
            {-6.52100048e-01, -1.38607428e-01},
            {-3.18764701e-01,  8.87287139e-01},
            { 6.96854372e-02, -9.40230521e-01},
            { 4.46086072e-01, -4.95429100e-01}
        };
        polygon_coordinates[Polygon::Icosahedron] = {
            {-8.94427401e-01,  0.00000000e+00},
            { 0.00000000e+00,  0.00000000e+00},
            { 7.23606342e-01, -5.25731196e-01},
            { 7.23606342e-01,  5.25731196e-01},
            {-2.76392839e-01, -8.50650862e-01},
            {-2.76392839e-01,  8.50650862e-01}
        };
        polygon_coordinates[Polygon::test] = {
            {-0.78947368,  0.52631579},
            { 0.71956327, -0.07106798},
            {-0.83085071,  0.16411866},
            { 0.05134629, -0.87225298},
            { 0.76631844, -0.6391243 },
            {-0.35677185,  0.70473452},
            {-0.59843068, -0.79522101},
            {-0.29036692, -0.49708905},
            {-0.66244178, -0.14539189},
            { 0.73033085,  0.25005566},
            { 0.9901128 ,  0.0237064 },
            { 0.2401403 , -0.23190504},
            {-0.24356815, -0.79728833},
            {-0.78495416, -0.60489685},
            { 0.17998138,  0.74263921},
            {-0.52734902,  0.42989847},
            {-0.24666497, -0.39700991},
            { 0.73528321, -0.01291032},
            { 0.53520874,  0.60655655},
            { 0.32039536,  0.94311487},
            {-0.55824441,  0.80374943},
            { 0.95502034,  0.2029875 },
            {-0.51220465, -0.22983242},
            {-0.90234999, -0.37391205},
            { 0.20163116, -0.92320135},
            {-0.7894378 ,  0.39952827},
            { 0.56505887, -0.27118692},
            {-0.28194829,  0.94581223},
            {-0.43702741, -0.62483872},
            { 0.0487374 , -0.14547679},
            {-0.92085544,  0.0107791 },
            { 0.52079749, -0.56779835},
            { 0.59814782, -0.75851087},
            {-0.07508506, -0.94131328},
            {-0.15800437,  0.23119067},
            { 0.25872507,  0.71775799},
            { 0.83117746, -0.43318884},
            { 0.00607969,  0.98422672},
            {-0.79853414, -0.59199361},
            { 0.81465343,  0.35671526},
            {-0.12407075,  0.4502174 },
            {-0.14144438,  0.95669616},
            { 0.93709725, -0.32477502},
            { 0.43284022,  0.46132957},
            {-0.83461502, -0.21058945},
            { 0.04541209, -0.8260639 },
            {-0.80879467,  0.52416916},
            { 0.93998817, -0.29865455},
            {-0.78571549,  0.12876106},
            { 0.11032759,  0.80168824},
            { 0.49907824, -0.86263992},
            {-0.46441649,  0.34216349},
            { 0.37691676, -0.47943857},
            {-0.5896296 ,  0.80577937},
            {-0.1382676 , -0.83054393},
            {-0.0959989 , -0.41577574},
            { 0.21666526, -0.66737826},
            { 0.90665782, -0.29265911},
            {-0.7387881 , -0.60218003},
            {-0.28743024, -0.51365067},
            {-0.97281803,  0.23016991},
            {-0.70024794, -0.02732258},
            { 0.66804122,  0.46605879},
            { 0.2368295 , -0.64993051},
            {-0.21691327, -0.76790176},
            {-0.76026287, -0.60747272},
            { 0.44161257,  0.87878276},
            {-0.75907322,  0.3325662 },
            {-0.11739113, -0.38933563},
            { 0.1267021 ,  0.2283998 },
            { 0.25664534,  0.87508575},
            {-0.43450207,  0.64269533},
            {-0.8623405 ,  0.07277181},
            { 0.93368285, -0.2999191 },
            {-0.21057294, -0.01173747},
            { 0.44864425, -0.89331679},
            {-0.5394668 ,  0.50279678},
            { 0.71255482, -0.45151713},
            {-0.15179426,  0.97324577},
            { 0.9226174 ,  0.28720156},
            { 0.54176207,  0.06106739},
            {-0.97505653,  0.22165457},
            { 0.06215722, -0.86824944},
            {-0.44214501, -0.88984543},
            {-0.65153555, -0.31256437},
            {-0.27878992, -0.11928233},
            { 0.35266608,  0.25712442},
            { 0.54410163, -0.26724821},
            {-0.22216102, -0.7881555 },
            { 0.8883393 ,  0.45352427},
            { 0.01496563,  0.83456093},
            {-0.4907137 ,  0.65846601},
            {-0.94923198, -0.13756806},
            {-0.74769062, -0.41482637},
            { 0.60544967, -0.02834844},
            { 0.06754417,  0.46308946},
            { 0.57396225,  0.74940962},
            {-0.87993516,  0.41561846},
            {-0.13639258,  0.73738229},
            { 0.79549249,  0.29128793}
        };
    }

    // --- //

    nlopt_info info = { r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta), r, phi, theta, this };
    double dphi = Rg2 / r;
    if (dphi > 0.5) {
        dphi = asin(dphi);
    }

    auto distance = [&](double r, double phi) {
      double vx, vy, vz;
      info.p->HCS_rphi(r, phi, vx, vy, vz);
      return sqrt((vx - info.x) * (vx - info.x) + (vy - info.y) * (vy - info.y) + (vz - info.z) * (vz - info.z));
    };

    double d_min = __DBL_MAX__;
    for (auto pc : polygon_coordinates[polygon]) {
        double next_d = distance(r + Rg2 * pc[0], phi + dphi * pc[1]);
        if (next_d < d_min) d_min = next_d;
    }
    return d_min;
}
    
double particle::get_HCS_distance() const {
  Vec target, p_cs;
  target.set_spherical(r, theta, phi);

  double rlow, rup;
  double phi0 = Phi0_S(theta);
  r_bound(r, phi, phi0, rlow, rup);
  if (rlow < 0) rlow = 1e-2*AU; // Avoid negative radius

  double phi_cs = phi, theta_cs = theta;
  if (theta_cs < pi / 2 - angle) theta_cs = pi / 2 - angle;
  else if (theta_cs > pi / 2 + angle) theta_cs = pi / 2 + angle;

  auto distance_iter = [&](Vec& point, int& iter) -> double {
    double diter = 1e5 * AU,
           diter_last = 1e5 * AU;

    point.set_spherical(point.len(), Theta_S(point.len(), point.phi()), point.phi()); // initialize the point to HCS.
    //cout << "begin:   " << r / AU << " " << theta / deg << " " << phi / deg << " | " << point.len() / AU << " " << point.theta() / deg << " " << point.phi() / deg << " -> " << diter / AU << endl;
    Vec dh = norm_vec(point);
    int viter = 0;
    while (diter == 1e5 * AU || diter_last == 1e5 * AU || fabs(diter_last - diter) / diter_last > 1e-4) {
      diter_last = diter;
      if (fabs(pi / 2 - theta) > angle - 0.5 * deg && fabs(pi / 2 - point.theta()) > angle - 0.5 * deg) {
        spiral_iterate(target, point, diter);
        //cout << "diter_s: " << r / AU << " " << theta / deg << " " << phi / deg << " | " << point.len() / AU << " " << point.theta() / deg << " " << point.phi() / deg << " -> " << diter / AU << endl;
        wave_iterate(target, point, diter);
        //cout << "diter_w: " << r / AU << " " << theta / deg << " " << phi / deg << " | " << point.len() / AU << " " << point.theta() / deg << " " << point.phi() / deg << " -> " << diter / AU << endl;
        dh = norm_vec(point);
        //cout << dh << endl;
      } else {
        point_iterate(target, point, dh, diter);
        //cout << "diter: " << r / AU << " " << theta / deg << " " << phi / deg << " | " << point.len() / AU << " " << point.theta() / deg << " " << point.phi() / deg << " -> " << diter / AU << endl;
      }
      //cout << diter_last / AU << " " << diter / AU << " " << (diter_last - diter) / diter_last << endl;
   }
    return (target - point).len();
  };

  int ilow = 0, imid = 0, iup = 0;

  //cout << "----------------------low------------------------" << endl;
  double dlow = 1e5 * AU;
  if ((rup - r) / (r - rlow) > 0.35) {
    p_cs.set_spherical(rlow, theta_cs, phi_cs);
    dlow = distance_iter(p_cs, ilow);
  }

  //cout << "----------------------mid------------------------" << endl;
  double dmid = 1e5 * AU;
  if (fabs(pi / 2 - theta) < angle) {
    p_cs.set_spherical(r, Theta_S(r, phi), phi);
    dmid = distance_iter(p_cs, imid);
  }

  //cout << "----------------------up------------------------" << endl;
  double dup = 1e5 * AU;
  if ((r - rlow) / (rup - r) > 0.35) {
    p_cs.set_spherical(rup, theta_cs, phi_cs);
    dup = distance_iter(p_cs, iup);
  }

  double sign = Theta_S(r, phi) < theta ? -1 : 1;
  return sign * fmin(fmin(dlow, dup), dmid);
}

void particle::step(const string& logname) {
  if (!fix_seed) {
    random_device rd;
    seed = rd();
  }
  mt19937 gen(seed);
  double mean = 0.0;
  double dev = 1.0;
  std::normal_distribution<double> dist(mean, dev);
  std::uniform_real_distribution<> dis(0, pi);

  double record_T = 0.;
  hcsform = Kota_Jokipii;

  ostringstream osname;
  if (!logname.empty()) osname << "s" << seed << "_" << logname;
  std::ofstream logfile(osname.str());
  if (logfile.is_open())
     logfile << "t[month],r[AU],theta[rad],phi[rad],Ek[GeV],drift,Vdr_gc[km/s]" << endl;

  auto r2V_r = [&](double r) {
    return r * r * Wind(r, theta, phi, angle);
  };

  auto dvdx = [&](double x, const function<double(double)>& v) {
    return (v(x * (1 + 1e-3)) - v(x)) / (1e-3 * x);
  };

  // theta = 1e-3;
  double Dt = 0;
  double M_p0 = sqrt(Ek * (Ek + 2. * mass));
  // theta = 0.00;
  // r = 5*AU;
  theta = Theta_S(r, phi)*(1-1e-2/2.);
  // double theta10 = theta;
  while (r<boundary) {//theta<pi/2
    Dt += dt;
    M_p = sqrt(Ek * (Ek + 2. * mass));
    rigidity = A / (Z * e) * M_p;
    V_p = M_p / (Ek + mass) * c_speed;
    Vs = Wind();

    // std::cout << "Mp:  " << r << "  " << M_p/GeV << "  " << Ek/GeV << "  " << mass/GeV << std::endl;
    // getchar();

    double cs = Theta_S(r, phi);
    // theta = cs*(1 + 1e-3);
    heaviside = Heav();
    Br = B_r(r, theta, phi, heaviside);
    Bp = B_p(r, theta, phi, heaviside);
    psi = atan(fabs(Bp / Br));
    double B = sqrt(Br * Br + Bp * Bp);

    double Br1 = B_r(r*(1+1e-3), theta*(1+1e-3), phi, heaviside);
    double Bp1 = B_p(r*(1+1e-3), theta*(1+1e-3), phi, heaviside);
    double psi1 = atan(fabs(Bp1/Br1));
    double B1 = sqrt(Br1 * Br1 + Bp1 * Bp1);

    k_rr = K_rr(r, theta, phi, psi, B, B0, M_p, M_p0, V_p);
    k_tt = K_tt(r, theta, phi, psi, B, B0, M_p, M_p0, V_p);
    k_pp = K_pp(r, theta, phi, psi, B, B0, M_p, M_p0, V_p);
        
    double k_rr1 = K_rr(r*(1+1e-3), theta*(1+1e-3), phi, psi1, B1, B0, M_p, M_p0, V_p);
    double k_tt1 = K_tt(r*(1+1e-3), theta*(1+1e-3), phi, psi1, B1, B0, M_p, M_p0, V_p);
    double k_pp1 = K_pp(r*(1+1e-3), theta*(1+1e-3), phi, psi1, B1, B0, M_p, M_p0, V_p);

    double dr2V_dr = dvdx(r, r2V_r);

    double gamma = tan(psi);
    double drift = 2 * M_p * V_p * r / (3 * Z * e * c_speed * Bn );
    double Vdr_gc =  drift / pow(1 + gamma * gamma, 2.) * (-1. * gamma ) / fabs(tan(theta));
    double Vdt_gc = A * polarity * drift / pow(1 + gamma * gamma, 2.) *  (2. + gamma * gamma) * gamma * heaviside;
    double Vdp_gc = A * polarity * drift / pow(1 + gamma * gamma, 2.) * gamma * gamma / fabs(tan(theta)) * heaviside;//;

    double delta = (Theta_S(r+0.1, phi) - Theta_S(r, phi)) / 0.1;
    if(delta<0) delta = -1.;
    else if(0<(delta)) delta = 1.;
    double beta = atan(Omega * r * sqrt(fabs(sin(angle) * sin(angle) - cos(cs) * cos(cs))) / (Vs * sin(psi) * sin(cs))) * delta;
    if(Z*polarity<0) beta = pi + beta;
    else if(0<Z*polarity) beta = beta;

    //if (fabs(r - 66.76 * AU) / r < 0.01) logfile << setprecision(15);
    //r = 66.7606607635679 * AU;
    //theta = 2.12431446451745;
    //phi = 2.4543874540613;
    double Rg = M_p / (B * Z * e * c_speed);
    double d_HCS = fabs(get_HCS_distance());
    //exit(0);
    double Vns = 0.;
    if(d_HCS<2.*Rg) Vns = (0.457 - 0.412 * d_HCS / Rg + 0.0915 * d_HCS * d_HCS / Rg / Rg) * V_p * A_drift;//

    double zonal = dis(gen);
    double Vdr_HCS = Vns * cos(beta) * sin(zonal) * A ;
    double Vdt_HCS = Vns * sin(beta) * A ;
    double Vdp_HCS = Vns * cos(beta) * cos(zonal) * A ;

    double Vdr = Vdr_gc + Vdr_HCS;
    double Vdt = Vdt_gc + Vdt_HCS;
    double Vdp = Vdp_gc + Vdp_HCS;

    double r0 = r;
    double theta0 = theta;
    double phi0 = phi;

    double dwr = dist(gen);
    double dwp = dist(gen);
    double dwt = dist(gen);
    if(3<fabs(dwr)) dwr = dist(gen);
    if(3<fabs(dwp)) dwp = dist(gen);
    if(3<fabs(dwt)) dwt = dist(gen);

    // r += 0.01 * AU;
    // r += -1. * Vs * dt;- Vdr_HCS +

// std::cout << "V:  " << Vs << "  " << V1 << std::endl;
// getchar();
    // logfile << Dt << "  " << r/AU << "  " << theta << "  " << Vs*dt << "  " << Vdr_gc*dt << "  " << 2. * k_rr / r *dt << "  " << sqrt(2. * fabs(k_rr) * dt) * dwr << "  " << k_rr << std::endl;
// 
// logfile << Dt << "  " << r/AU << "  " << theta << "  " << Vs*dt << "  " << Vdr_gc*dt << "  " << 1. / r / r * (r * r * k_rr - r10 * r10 * k_rr10) / (r - r10) << std::endl;
//  
    if(r<4.*AU) dwr = fabs(dwr);
    r += (- Vs - Vdr
          + 1. / r / r * (r*(1+1e-3) * r*(1+1e-3) * k_rr1 - r * r * k_rr) / (r*1e-3)) * dt +
         sqrt(2. * fabs(k_rr) * dt) * dwr;

    // r += - Vdr_HCS * dt;

    theta += ( - Vdt / r
              + 1. / r / r * (sin(theta * (1+1e-3)) * k_tt1 - sin(theta) * k_tt) / (theta * 1e-3)
             ) *dt 
              + 1. / r * sqrt(2. * fabs(k_tt) * dt) * dwt ;

    // theta += - Vdt_HCS / r * dt;
    // theta10 += - Vdt_HCS / r * dt;

    phi += - Vdp / (r * sin(theta)) * dt +
           sqrt(2. * fabs(k_pp) * dt) * dwp / (r * sin(theta));

    double r1 = r0 * (1 + 1e-3);
    double V1 = Wind(r1, theta0, phi0, angle);
    double E = Ek + mass;
    double p2 = E * E - mass * mass;

    Ek += dr2V_dr / 3 / r0 / r0 * p2 / E * dt;

    //Ek += p2 / E * 1. / (3 * r0 * r0) * (r1 * r1 * V1 - r0 * r0 * Vs) / (r1 - r0) * dt;
    //Ek += 2. / 3. * Vs / r0 * p2 / E * dt;
    if (r < 0.) {
      r = 0.;
      break;
    }

    if (theta < 0.) {
      theta = fabs(theta);
      phi += pi;
    } 
    else if (pi < theta) {
      theta = 2. * pi - theta;
      phi += pi;
    }
    
    if(theta < 5*1e-3){
      theta = 5*1e-3;
    }
    else if(3.14-5*1e-3<theta){
      theta = 3.14-5*1e-3;
    }

    if (phi < 0 || 2 * pi < phi)
      phi -= floor(phi / (2 * pi)) * 2 * pi;

    // logfile << r / AU << "  " << theta << "  " << 1. * Vdt_gc / r*dt << "  " << Vdt_HCS / r*dt << "  " << 1. / (r * r * sin(theta)) * cos(theta) * k_tt * dt << "  " << 1. / r * sqrt(2. * fabs(k_tt) * dt) * dwt << std::endl;
    // logfile << r/AU << "  " << M_p / GeV << "  " << Ek/GeV << std::endl;
    // if(60*60*24*365*1.5<record_T) break;
    // logfile << r/AU << "  " << theta << "  " << Vdt_gc / r * dt << "  " << 1. / (r * r * sin(theta)) * cos(theta) * k_tt * dt << "  " << k_tt << "  " << 1. / r * sqrt(2. * fabs(k_tt) * dt) * dwt << std::endl;
    // logfile << r/AU << "  " << theta << "  " << Vdr_HCS << "  " << Vdr_gc << "  " << Vdt_HCS << "  " << Vdt_gc << "  " << Vs << "  " << drift << std::endl;
    if (logfile.is_open())
       logfile << Dt/day/30 << "," << r/AU << "," << theta << "," << phi  << "," << Ek / GeV << "," << drift << "," << Vdr_gc << endl;

    //if (r > 120 * AU) {
    //  cout << seed << " " << r / AU << " " << theta / deg << " " << phi / deg << endl;
    //}
    //assert(r < 120 * AU && "The particle should not get out of the solar system too much");
    // if(theta<pi*10./180. || pi*170./180.<theta) break;
    // if(24*30*120*3600.<Dt){
    //   std::cout << "Out time" << std::endl;
    //   break;
    // }
  }
  // if(90 < r/AU)  {
  //   logfile << r/AU << "  " << theta << "  " << phi << "  " << Ek / GeV << std::endl;
  // }
  if (logfile.is_open()) logfile.close();
}
