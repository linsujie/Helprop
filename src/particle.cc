#include <cassert>
#include <iomanip>

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

particle::particle() {}

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

const double particle::Wind() {
  double value;
  if (0. <= theta && theta <= pi / 2.) {
    value = (1.475 - 0.4 * tanh(6.8 * ((theta - pi / 2) + (15 * deg + angle)))) *
                        (3.5 / 5. - 1.5 / 5. * tanh((r - 95 * AU) / 1.2 / AU));
  } else if (pi / 2. < theta && theta <= pi) {
    value = (1.475 + 0.4 * tanh(6.8 * ((theta - pi / 2) - (15 * deg + angle)))) *
                        (3.5 / 5. - 1.5 / 5. * tanh((r - 95 * AU) / 1.2 / AU));
  }

  //     for(int i=0;i<9;i++){
  //         double a = (5+10*i)*pi/180;
  //         cout << a << "   " << 1.475 - 0.4 * tanh(6.8*((a - pi/2.) +
  //         (15.*deg + angle))) * (3.5/5. - 1.5/5.*tanh((r-95*AU)/1.2/AU)) <<
  //         endl;
  //     }
  //     for(int i=9;i<18;i++){
  //         double a = (5+10*i)*pi/180;
  //         cout << a << "   " << 1.475 + 0.4 * tanh(6.8*((a - pi/2.) -
  //         (15.*deg + angle))) * (3.5/5. - 1.5/5.*tanh((r-95*AU)/1.2/AU)) <<
  //         endl;
  //     }
  // getchar();
  return value * 400 * (km / sec);
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

  return Theta_S_Jokipii_Thomas(r, phi);
  assert(false && "hcsform not supported");
  return 0;
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

const double particle::B_r(const double& r, const double& heaviside) {
  // double r0 = 1 * AU;
  // std::cout << "D:   " << heaviside << "   " << polarity << "   " << r << "
  // " << r0 << std::endl;
  //     getchar();

  return Bn * heaviside * polarity * pow(1. / r, 2.);
}

const double particle::B_p(const double& r, const double& theta, const double& heaviside) {
  // std::cout << "D:   " << Omega << "   " << sin(theta) << "   " << Vs << " "
  // << polarity<< std::endl; getchar();
  // double r0 = 1 * AU;

  //  std::cout << "D:   " << "  " << Omega * sin(theta) * heaviside * polarity / Vs * r  << std::endl; getchar();
  // std::cout << Bn << "  " << r << "  " << sin(theta ) << "  " << a << std::endl;
  // getchar();
  return -1. * Bn / r * Omega * sin(theta) * heaviside * polarity / Vs;
}

const double particle::K_rr() {
  double kx = D * pow(Ek / GeV, indexA);
  double ky = 0.02 * kx;
  double kr = kx * pow(cos(psi), 2.) + ky * pow(sin(psi), 2.);

  // cout << "VD :  " << kx << "  " << D << "  " << ky << "  " <<  pow(Ek/GeV,
  // indexA) << "   " << kr << endl; getchar();

  return kr;
}

const double particle::K_tt(double theta) {
  double kx = D * pow(Ek / GeV, indexA);
  double ky = 0.02 * kx;
  double kz;
  if (theta < pi / 2. && 0 <= theta)
    kz = ky * (2. - 1. * tanh(8 * ((theta + (-90. + 35.) * deg))));
  else if (pi / 2. <= theta && theta <= pi)
    kz = ky * (2. + 1. * tanh(8 * ((theta + (-90. - 35.) * deg))));

  // cout << "VD :  " << ky<< "  " << kz << "  " << (2. + 1.*tanh(8*((theta +
  // (- 90. - 35.)*pi/180.)))) << endl;

  //         getchar();
  return kz;
}

const double particle::K_pp() {
  double kx = D * pow(Ek / GeV, indexA);
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

double particle::get_HCS_distance_old() const {
  nlopt_info info = { r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta), r, phi, theta, this };
  nlopt::opt opt(nlopt::LD_MMA, 2);
  opt.set_min_objective(distance_to_point, &info);
  opt.set_xtol_rel(1e-8);

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

double particle::spiral_iterate(const Vec& target_point, Vec& p_cs) const {
  double ov = Omega / Vs_eq;

  SpiralVdot vdot(ov, p_cs, target_point);

  double vdot0 = vdot(vdot.r_cs0);
  if (vdot0 == 0) return 0;

  double r1;
  double dangle = - asin(vdot0);

  do {
    dangle *= 2;
    r1 = (vdot.phi0 - vdot.phi_cs0 - dangle * 2) / ov;
  } while (vdot(r1) * vdot0 > 0);

  double rh = ridders_method(vdot, vdot.r_cs0, r1, 1e-3);

  p_cs.set_spherical(rh, vdot.theta_cs, vdot.phi0 -  rh * ov);
  //cout << "++++ " << rh / AU << " " << vdot(rh) << " | " << (target_point - vdot.p_cs0).len() / AU << " " << (target_point - p_cs).len() / AU << endl;
  return (target_point - p_cs).len();
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

  double operator()(double r) const {
    Vec p_cs;
    if (r == r_cs0) p_cs = p_cs0;
    else
      p_cs.set_spherical(r, p->Theta_S(r, phi_cs), phi_cs);

    //cout << "-- " << r / AU << " " << (target_point - p_cs).dot(tangent_vec(r, p_cs)) / AU << endl;
    Vec dl = target_point - p_cs;
    return dl.dot(tangent_vec(r, p_cs)) / dl.len();
  }
};

double particle::wave_iterate(const Vec& target_point, Vec& p_cs) const {
  double ov = Omega / Vs_eq;

  WaveVdot vdot(ov, p_cs, target_point, this);

//  Vec p_cs0 = p_cs;
//  double r0 = p_cs.len();
//  double phi_cs = p_cs.phi();
//
//  bool update_p_cs = false;
//
//  double rphi0 = sqrt(p_cs.x * p_cs.x + p_cs.y * p_cs.y);
//  double sphi = p_cs.y / rphi0, cphi = p_cs.x / rphi0;
//
//  int iter = 0;
//  auto tangent_vec = [&](const Vec& point) -> Vec {
//    double r = point.len();
//    double phi0_prime = phi_cs + r * ov;
//    double rphi = sqrt(point.x * point.x + point.y * point.y);
//    double stheta = rphi / r, ctheta = point.z / r;
//
//    Vec dr(stheta * cphi, stheta * sphi, ctheta);
//    Vec dtheta(r * ctheta * cphi, r * ctheta * sphi, - r * stheta);
//
//    /************************************************************
//     * Jokipii_Thomas:  sin(pi / 2 - theta) = sin(angle) * sin(phi + ov * r)
//     * Jokipii_Thomas:  - cos(pi / 2 - theta) * Delta_theta = sin(angle) * cos(phi + ov * r) * Delta_r * ov
//     * Jokipii_Thomas:  Delta_theta = - sin(angle) * cos(phi + ov * r) / cos(pi / 2 - theta) * Delta_r * ov
//     * 
//     * Kota_Jokipii:    tan(pi / 2 - theta) =  tan(angle) * sin(phi + ov * r)
//     * Kota_Jokipii:    - cos(pi / 2 - theta)^-2 * Delta_theta = tan(angle) * cos(phi + ov * r) * Delta_r * ov
//     * Kota_Jokipii:    Delta_theta = - tan(angle) * cos(phi + ov * r) * cos(pi / 2 - theta)^2 * Delta_r * ov
//     * 
//     *       dtheta_phi = Delta_theta * dtheta + Delta_r * dr
//     * Jokipii_Thomas:  ~ - sin(angle) * cos(phi + ov * r) / cos(pi / 2 - theta) * ov * dtheta + dphi
//     * Kota_Jokipii:    ~ - tan(angle) * cos(phi + ov * r) * cos(pi / 2 - theta)^2 * ov * dtheta + dphi
//     ************************************************************/
//
//    Vec dv;
//    if (hcsform == Jokipii_Thomas) dv = - sin(angle) * cos(phi0_prime) / stheta * ov * dtheta + dr;
//    else if (hcsform == Kota_Jokipii) dv = - tan(angle) * cos(phi0_prime) * stheta * stheta * ov * dtheta + dr;
//    else assert(false && "Unsuported hcsform");
//
//    dv.normalize();
//    iter++;
//    return dv;
//  };
//
//  map<double, double> vdot_values;
//  auto vdot = [&](double rprime) -> double {
//    if (vdot_values.find(rprime) != vdot_values.end())
//      return vdot_values.at(rprime);
//
//    if (update_p_cs)
//      p_cs.set_spherical(rprime, Theta_S(rprime, phi_cs), phi_cs);
//
//    Vec dv = tangent_vec(p_cs);
//    double res = (target_point - p_cs).dot(dv);
//
//    vdot_values[rprime] = res;   
//    return res;
//  };

  double vdot0 = vdot(vdot.r_cs0);
  if (vdot0 == 0) return 0;

  double vdr = vdot0 * (target_point - vdot.p_cs0).len();
  if (fabs(vdr) > pi / 2 / ov) vdr = pi / 2 / ov * (vdr > 0 ? 1 : -1);

  int ir = 0;
  double r1;
  do {
    ir++;
    r1 = vdot.r_cs0 + ir * vdr;
  } while (vdot(r1) * vdot0 > 0);

  double rh = ridders_method(vdot, vdot.r_cs0, r1, 1e-3);
  p_cs.set_spherical(rh, Theta_S(rh, vdot.phi_cs), vdot.phi_cs);
  //cout << "==== " << rh / AU << " " << vdot(rh) << " | " << (target_point - vdot.p_cs0).len() / AU << " " << (target_point - p_cs).len() / AU << endl;
  return (target_point - p_cs).len();
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

double particle::point_iterate(const Vec& target_point, Vec& p_cs, Vec& dh) const {
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
  return (target_point - p_cs).len();
}

double particle::get_HCS_distance() const {
  Vec target, p_cs;
  target.set_spherical(r, theta, phi);

  double rlow, rup;
  double phi0 = Phi0_S(theta);
  r_bound(r, phi, phi0, rlow, rup);

  double phi_cs = phi, theta_cs = theta;
  if (theta_cs < pi / 2 - angle) theta_cs = pi / 2 - angle;
  else if (theta_cs > pi / 2 + angle) theta_cs = pi / 2 + angle;

  auto distance_iter = [&](Vec& point, int& iter) -> double {
    double diter = 1e5 * AU,
           diter_last = 1e5 * AU;

    Vec dh = norm_vec(point);
    int viter = 0;
    while (diter == 1e5 * AU || diter_last == 1e5 * AU || fabs(diter_last - diter) / diter_last > 1e-4) {
      diter_last = diter;
      diter = point_iterate(target, point, dh);
      //cout << "diter: " << r / AU << " " << theta / deg << " " << phi /deg << " | " << point.len() / AU << " " << point.theta() / deg << " " << point.phi() / deg << " -> " << diter / AU << endl;
      if (fabs(pi / 2 - theta) > angle - 0.5 * deg && fabs(pi / 2 - point.theta()) > angle - 0.5 * deg) {
        diter = spiral_iterate(target, point);
        //cout << "diter_s: " << r / AU << " " << theta / deg << " " << phi /deg << " | " << point.len() / AU << " " << point.theta() / deg << " " << point.phi() / deg << " -> " << diter / AU << endl;
        diter = wave_iterate(target, point);
        //cout << "diter_w: " << r / AU << " " << theta / deg << " " << phi /deg << " | " << point.len() / AU << " " << point.theta() / deg << " " << point.phi() / deg << " -> " << diter / AU << endl;
        dh = norm_vec(point);
        //cout << dh << endl;
      }
      //cout << diter_last / AU << " " << diter / AU << " " << (diter_last - diter) / diter_last << endl;
   }
    return (target - point).len();
  };

  int ilow = 0, imid = 0, iup = 0;
  //cout << "----------------------------------------------" << endl;
  double dlow = 1e5 * AU;
  if ((rup - r) / (r - rlow) > 0.35) {
    p_cs.set_spherical(rlow, theta_cs, phi_cs);
    dlow = distance_iter(p_cs, ilow);
  }

  //cout << "----------------------------------------------" << endl;
  double dmid = 1e5 * AU;
  if (fabs(pi / 2 - theta) < angle) {
    p_cs.set_spherical(r, Theta_S(r, phi), phi);
    dmid = distance_iter(p_cs, imid);
  }

  //cout << "----------------------------------------------" << endl;
  double dup = 1e5 * AU;
  if ((r - rlow) / (rup - r) > 0.35) {
    p_cs.set_spherical(rup, theta_cs, phi_cs);
    dup = distance_iter(p_cs, iup);
  }

  double sign = Theta_S(r, phi) < theta ? -1 : 1;
  return sign * fmin(fmin(dlow, dup), dmid);
}


double wind1(double r, double theta, double phi, double angle){
  double value;
    if (0. <= theta && theta <= pi / 2.) {
      value = (1.475 - 0.4 * tanh(6.8 * ((theta - pi / 2) + (15 * deg + angle)))) * (3.5 / 5. - 1.5 / 5. * tanh((r - 95 * AU) / 1.2 / AU));
    } else if (pi / 2. < theta && theta <= pi) {
      value = (1.475 + 0.4 * tanh(6.8 * ((theta - pi / 2) - (15 * deg + angle)))) * (3.5 / 5. - 1.5 / 5. * tanh((r - 95 * AU) / 1.2 / AU));
    }
    return  value * 400 * (km / sec);
}

void particle::step() {
  random_device rd;
  mt19937 gen(rd());
  double mean = 0.0;
  double dev = 1.0;
  std::normal_distribution<double> dist(mean, dev);
  std::uniform_real_distribution<> dis(0, pi);

  double record_T = 0.;
  hcsform = Kota_Jokipii;

  std::ofstream f2("file2");  //, std::ios::app
  // f2.open("file2");
  // theta = 1e-3;
  double Dt = 0;
  double M_p0 = sqrt(Ek * (Ek + 2. * mass));
  // theta = 0.00;
  // r = 5*AU;
  while (r<boundary) {//theta<pi/2
  // theta += 0.001;

  // r = 63.4163*AU;
  // theta = 3.135;
    Dt += dt;
    M_p = sqrt(Ek * (Ek + 2. * mass));
    rigidity = A / (Z * e) * M_p;
    V_p = M_p / (Ek + mass) * c_speed;
    Vs = Wind();

    // std::cout << "Mp:  " << r << "  " << M_p/GeV << "  " << Ek/GeV << "  " << mass/GeV << std::endl;
    // getchar();

    double cs = Theta_S(r, phi);
    // theta = cs + 1e-3;
    heaviside = Heav();
    Br = B_r(r, heaviside);
    Bp = B_p(r, theta, heaviside);
    psi = atan(fabs(Bp / Br));

    double B = sqrt(Br * Br + Bp * Bp);

    // k_rr = K_rr();
    // k_tt = K_tt(theta);
    // k_pp = K_pp();

    // Br1 = B_r(r+1e-5, heaviside);
    // Bp1 = B_p(r+1e-5, theta+1e-5, heaviside);
    // psi1 = atan(fabs(Bp1 / Br1));

    double lamda0 = 0.05*AU;
    double lamda_p = lamda0 * B0 / B * M_p / M_p0;
    double kx = lamda_p * V_p / 3.;
    double ky = 0.02 * kx;
    double k_rr = kx * pow(cos(psi), 2.) + ky * pow(sin(psi), 2.);
    double k_tt = ky;


    // std::cout << "k:  " << k_rr << "  " << k_rr1 << "  " << r << "  " << Ek << std::endl;
    // getchar();

    // double lamda_p1 = lamda0 * B1 / B0 * M_p / M_p0;
    // double kx1 = lamda_p1 * V_p / 3.;
    // k_rr1 = kx1 * pow(cos(psi1), 2.) + ky * pow(sin(psi1), 2.);


    double gamma = tan(psi);//r * Omega * sin(theta) / Vs;

    double drift = 2 * M_p * V_p * r / (3 * Z * e * c_speed * Bn );//B0 * r0 * r0
    // double drift = 2 * M_p * V_p * r / (3 * e);
    Vdr_gc = drift / pow(1 + gamma * gamma, 2.) * (-1. * gamma ) / tan(theta);//; * heaviside polarity * 
    Vdt_gc = A * polarity * drift / pow(1 + gamma * gamma, 2.) *  (2. + gamma * gamma) * gamma * heaviside;
    Vdp_gc = polarity * drift / pow(1 + gamma * gamma, 2.) * gamma * gamma / tan(theta) * heaviside;//;
    // std::cout << "Z:  " << Z << "  " << polarity << "  " << Vdt_gc << "  " << heaviside << "  " << drift / pow(1 + gamma * gamma, 2.) *  (2. + gamma * gamma) * gamma << std::endl;
    // getchar();

    double delta = (Theta_S(r+0.1, phi) - Theta_S(r, phi)) / 0.1;
    if(delta<0) delta = -1.;
    else if(0<(delta)) delta = 1.;
    double beta = atan(Omega * r * sqrt(fabs(sin(angle) * sin(angle) - cos(cs) * cos(cs))) / (Vs * sin(psi) * sin(cs))) * delta;
    if(Z*polarity<0) beta = pi + beta;
    else if(0<Z*polarity) beta = beta;
    // else if(0<Z*polarity && delta<0) beta = -1. * beta;
    // else if(0<Z*polarity && 0<delta) beta = 1. * beta;
    double Rg = M_p / (B * Z * e * c_speed);
// std::cout << "begin   " << r/AU << "  " << theta << "  " << phi << std::endl;
    double d_HCS = fabs(get_HCS_distance());
// std::cout << "begin1   " << std::endl;
    double Vns = 0.;
    if(d_HCS<2.*Rg) Vns = (0.457 - 0.412 * d_HCS / Rg + 0.0915 * d_HCS * d_HCS / Rg / Rg) * V_p / 2. ;//

    double zonal = dis(gen);
    Vdr_HCS = Vns * cos(beta)  * Z ;//* e * Bn* sin(zonal)
    Vdt_HCS = Vns * sin(beta) * Z ;//* e * Bn
    Vdp_HCS = Vns * cos(beta) * cos(zonal) * Z ;//* e * Bn
// std::cout << "drift:  " << drift << "  " << Bn << "  " << M_p << "  " << V_p << "  " << r << "  " << e << "  " << c_speed << std::endl;
// getchar();

    // std::cout << "V:  " << V_p << "  " << Vns << "  " << Vdr_HCS << std::endl;
    // getchar();

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
    // f2 << Dt << "  " << r/AU << "  " << theta << "  " << Vs*dt << "  " << Vdr_gc*dt << "  " << 2. * k_rr / r *dt << "  " << sqrt(2. * fabs(k_rr) * dt) * dwr << "  " << k_rr << std::endl;
// 
    r += (-1. * Vs + fabs(Vdr_gc) - Vdr_HCS + 2. * k_rr / r) * dt +
         sqrt(2. * fabs(k_rr) * dt) * dwr;

    theta += ( - Vdt_gc / r - Vdt_HCS / r
              + 1. / (r * r * sin(theta)) * cos(theta) * k_tt
             ) *dt 
              + 1. / r * sqrt(2. * fabs(k_tt) * dt) * dwt;

    // phi += (-1. * Vdp_gc - Vdp_HCS) / (r * sin(theta)) * dt +
    //        sqrt(2. * fabs(k_pp) * dt) * dwp / (r * sin(theta));

//     std::cout << "test:  " << (-1. * Vdp_gc - Vdp_HCS) / (r * sin(theta)) * dt << "  " << sqrt(2. * fabs(k_pp) * dt) * dwp / (r * sin(theta)) << std::endl;
// getchar();
    // M_p += 2. * Vs / (3. * r) * M_p * M_p / (Ek + mass)  * dt;
    // Ek = sqrt(M_p * M_p + mass * mass) - mass;

    double r1 = r;
    double V1 = wind1(r1, theta0, phi0, angle);

    Ek += 1. / (3 * r0 * r0) * (r1 * r1 * V1 - r0 * r0 * Vs) / (r1 - r0) * Ek * dt;
    // Ek += 2. / 3. * Vs / r0 * Ek * dt;
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

    if (phi < 0.) phi = 2. * pi + phi;
    else if (2. * pi < phi) phi -= 2. * pi;

    // if(10*30*24*3600<Dt) {
    //   std::cout << "break" << std::endl;make
    //   continue;
    // }

    // f2 << r/AU << "   " << theta << "   " << Vdr_HCS << "  " << Vdt_HCS << "  " << Vdp_HCS << "   " << heaviside << "  " << beta << "  " << delta << std::endl;
    // f2 << r / AU << "  " << theta << "  " << Vs*dt << "  " << Vdr_gc*dt << "  " << Vdr_HCS*dt << "  " << sqrt(2. * fabs(k_rr) * dt) * dwr << std::endl;
    // f2 << r / AU << "  " << theta << "  " << 1. * Vdt_gc / r*dt << "  " << Vdt_HCS / r*dt << "  " << 1. / (r * r * sin(theta)) * cos(theta) * k_tt * dt << "  " << 1. / r * sqrt(2. * fabs(k_tt) * dt) * dwt << std::endl;
    // f2 << r/AU << "  " << M_p / GeV << "  " << Ek/GeV << std::endl;
    // if(60*60*24*365*1.5<record_T) break;
    // f2 << r/AU << "  " << theta << "  " << Vdt_gc / r * dt << "  " << 1. / (r * r * sin(theta)) * cos(theta) * k_tt * dt << "  " << k_tt << "  " << 1. / r * sqrt(2. * fabs(k_tt) * dt) * dwt << std::endl;
    // f2 << r/AU << "  " << theta << "  " << Vdr_HCS << "  " << Vdr_gc << "  " << Vdt_HCS << "  " << Vdt_gc << "  " << Vs << "  " << drift << std::endl;
    // f2 << r/AU << "  " << theta << "  " << Vdr_gc*dt << "  " << Vdr_HCS*dt << "  " << 2. * k_rr / r*dt << "  " << sqrt(2. * fabs(k_rr) * dt) << std::endl;
    // f2 << r/AU << "  " << k_rr << "  " << psi << std::endl;
    f2 << r/AU << "  " << theta << "  " << Ek / GeV << "  " << drift << "  " << Vdr_gc << std::endl;
    // getchar();
  }
  // if(50 < r/AU)  {
  //   f2 << r/AU << "  " << theta << "  " << Ek / GeV << std::endl;
  // }
  std::cout << "get One particle." << std::endl;
  f2.close();
  getchar();
}