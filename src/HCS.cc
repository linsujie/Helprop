#include <fstream>
#include <vector>
#include <map>
#include <iomanip>
#include <cassert>
#include <iomanip>
#include <filesystem>
#include <regex>
#include <mutex>

#include "HCS.h"
#include "Unit.h"
#include "root_finding.h"
#include "newton_solver.h"
#include "nlopt.hpp"

using namespace std;
using namespace Unit;

double HCS::angle = 45 * Unit::deg;
HCS::HCSFORM HCS::hcsform = Jokipii_Thomas;
const double HCS::Omega = 2*Unit::pi/27.5/Unit::day;
std::vector<double> HCS::angle_axis;
std::vector<std::string> HCS::table_names;
std::vector<KDInterp*> HCS::tables;

HCS::HCS(double Vs_eq_, const std::string& table_dir) : Vs_eq(Vs_eq_), interp(filesystem::exists(table_dir)) {
  if (interp && tables.empty()) init_tables(table_dir);
}

HCS::~HCS() {}

void HCS::init_tables(const std::string& table_dir) {
  ifstream namelist((table_dir + "/namelist.txt").c_str());
  string line;

  double alow = -1, aup = -1, last_up = -1;
  angle_axis.clear();
  while (getline(namelist, line)) {
    table_names.push_back(line);

    regex r(".+/([0-9.]+)_([0-9.]+).bson");
    smatch m;
    assert(regex_match(line, m, r));

    alow = atof(m[1].str().c_str());
    aup = atof(m[2].str().c_str());

    assert(last_up == -1 || last_up == alow);
    last_up = aup;

    angle_axis.push_back(alow * deg);
  }

  angle_axis.push_back(aup * deg);

  tables.resize(table_names.size(), NULL);
}

std::string doubleToBinaryString(double value) {
    // 创建一个大小为 sizeof(double) 的字节数组
    unsigned char* bytes = reinterpret_cast<unsigned char*>(&value);
    std::string binaryString;

    // 遍历每个字节，并将其转换为8位二进制字符串
    for (int i = sizeof(double) - 1; i >= 0; --i) {
        // 将每个字节转换为8位二进制字符串并拼接到结果中
        binaryString += std::bitset<8>(bytes[i]).to_string();
    }

    return binaryString;
}

double binaryToDouble(const std::string& binaryString) {
    // 确保输入的二进制字符串长度为64位
    if (binaryString.length() != 64) {
        throw std::invalid_argument("Binary string must be 64 bits long.");
    }

    // 将二进制字符串转换为字节数组
    unsigned char bytes[sizeof(double)];
    for (size_t i = 0; i < sizeof(double); ++i) {
        // 提取每8位并转换为字节
        bytes[sizeof(double) - 1 - i] = static_cast<unsigned char>(std::bitset<8>(binaryString.substr(i * 8, 8)).to_ulong());
    }

    // 将字节数组重新解释为double
    return *reinterpret_cast<double*>(bytes);
}

struct nlopt_info {
  double x, y, z;
  double r, phi, theta;
  const HCS *h;
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
    info->h->rphi(r, phi, vx, vy, vz);
    return sqrt((vx - info->x) * (vx - info->x) + (vy - info->y) * (vy - info->y) + (vz - info->z) * (vz - info->z));
  };

  double d = distance(x[0], x[1]);
  double ddr = distance(x[0] + 1e-10 * AU, x[1]);
  double ddp = distance(x[0], x[1] + 1e-10);

  grad = { (ddr - d) / (1e-10 * AU), (ddp - d) / (1e-10) };

  return d;
}

double HCS::get_distance_old(double r, double theta, double phi, double ftol_abs) const {
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
  rphi(r, phi, vx, vy, vz);
  double sign = info.z > vz ? 1 : -1;

  //cout << "dlow: " << dlow / AU << " dmid: " << dmid / AU << " dup: " << dup / AU << endl;

  return sign * fmin(fmin(dlow, dup), dmid);
}

void HCS::rphi(const double &r, const double &phi, double& x, double& y, double& z) const {
  double cs0 = Theta_S(r, phi);
  x = r * sin(cs0) * cos(phi);
  y = r * sin(cs0) * sin(phi);
  z = r * cos(cs0);
}

void HCS::r_bound(double r, double phi, double phi0, double& rlow, double& rup) const {
  double Tr = 2 * pi * Vs_eq / Omega;

  rlow = (phi0 - phi + Omega * (t - t0)) * Vs_eq / Omega;
  rlow += floor((r - rlow) / Tr) * Tr;
  rup = rlow + Tr;
}

double HCS::Theta_S_Jokipii_Thomas(double phi0) {
  return pi / 2 - asin(sin(angle) * sin(phi0));
//  static fcache theta_jokipii_thomas([](double x) { return pi / 2 - asin(sin(angle) * sin(x)); }, 1000000);
//
//  return theta_jokipii_thomas(phi0);
}
double HCS::Phi0_S_Jokipii_Thomas(double theta) {
  double ratio = sin(pi / 2 - theta) / sin(angle);
  if (ratio >= 1) return pi / 2;
  else if (ratio <= -1) return -pi / 2;

  return asin(ratio);
}

double HCS::Theta_S_Kota_Jokipii(double phi0) {
  return pi / 2 - atan(tan(angle) * sin(phi0));
//  static fcache theta_kota_jokipii([](double x) { return pi / 2 - atan(tan(angle) * sin(x)); }, 1000000);
//
//  return theta_kota_jokipii(phi0);
}
double HCS::Phi0_S_Kota_Jokipii(double theta) {
  double ratio = tan(pi / 2 - theta) / tan(angle);
  if (ratio >= 1) return pi / 2;
  else if (ratio <= -1) return -pi / 2;

  return asin(ratio);
}

double HCS::Theta_S(double r, double phi) const {
  if (hcsform == Jokipii_Thomas)
    return Theta_S_Jokipii_Thomas(r, phi);
  else if (hcsform == Kota_Jokipii)
    return Theta_S_Kota_Jokipii(r, phi);

  assert(false && "hcsform not supported");
  return 0;
}
extern "C" double Theta_S_C(double r, double phi) {
  return HCS(430 * km / sec, "").Theta_S(r * AU, phi);
}
double HCS::Phi0_S(double theta) {
  if (hcsform == Jokipii_Thomas)
    return Phi0_S_Jokipii_Thomas(theta);
  else if (hcsform == Kota_Jokipii)
    return Phi0_S_Kota_Jokipii(theta);

  assert(false && "hcsform not supported");
  return 0;
}

class SpiralVdot {
  public:
  const Vec p_cs0, target_point;
  const double r_cs0, phi_cs0;
  double ov, theta_cs, phi0, ctheta, stheta;
  bool pflag = false;
  SpiralVdot(double ov_, const Vec& p_cs, const Vec& target_point_) :
    p_cs0(p_cs), target_point(target_point_),
    r_cs0(p_cs.len()), phi_cs0(p_cs.phi()),
    ov(ov_), theta_cs(p_cs.theta()), phi0(phi_cs0 + r_cs0 * ov), ctheta(p_cs.z / r_cs0), stheta(sqrt(1 - ctheta * ctheta)) {}

  Vec tangent_vec(double r, const Vec& point) const {
    double rphi = fabs(r) * stheta;
    double sphi = point.y / rphi, cphi = point.x / rphi;

    Vec dr(stheta * cphi, stheta * sphi, r < 0 ? - ctheta : ctheta);
    Vec dphi(-rphi * sphi, rphi * cphi, 0);

    /************************************************************
     * Delta_phi + Omega / V * Delta_r = 0
     * Delta_phi = - Omega / V * Delta_r
     * dr_phi = Delta_phi * dphi + Delta_r * dr
     *        ~ - Omega / V * dphi + dr
     ************************************************************/
    Vec dv = - ov * dphi + dr;
    dv.normalize();
    return r > 0 ? dv : - dv;
  }

  double operator()(double r) const {
    Vec p_cs;
    if (r == r_cs0) p_cs = p_cs0;
    else {
      double phi = phi0 - fabs(r * ov);
      p_cs.set_spherical(r, theta_cs, phi);
    }

    Vec dl = target_point - p_cs;
    if (pflag) {
      cout << "--SVdot " << p_cs / AU << " " << tangent_vec(r, p_cs) / 10 << endl;
    }
    //cout << fabs(r) / AU << " " << p_cs.theta() << " " << p_cs.phi() << " " << dl.dot(tangent_vec(r, p_cs)) / dl.len() << endl;
    return dl.dot(tangent_vec(r, p_cs)) / dl.len();
  }
};

bool HCS::spiral_iterate(const Vec& target_point, Vec& p_cs, double& diter) const {
  double ov = Omega / Vs_eq;

  SpiralVdot vdot(ov, p_cs, target_point);

  double vdot0 = vdot(vdot.r_cs0);
  if (vdot0 == 0) return false;

  double r1;
  double dangle = - asin(vdot0) * fmin(diter / vdot.r_cs0, 1); // the dangle should be smaller when the distance is much small than r_cs

  if (fabs(dangle) > 5 * deg) dangle = dangle > 0 ? 5 * deg : - 5 * deg; // Avoid the dangle too large

  int id = 0;
  do {
    id++;
    r1 = (vdot.phi0 - vdot.phi_cs0 - id * dangle) / ov;
  } while (vdot(r1) * vdot0 > 0);

  if (fabs(r1 - vdot.r_cs0) / vdot.r_cs0 < 1e-9) return false;

  double rh = ridders_method(vdot, vdot.r_cs0, r1, 1e-3);

  p_cs.set_spherical(rh, vdot.theta_cs, vdot.phi0 -  fabs(rh * ov));

  double diter_next = (target_point - p_cs).len();
//  if (diter_next > diter) {
//    cout << endl;
//    cout << setprecision(15) << target_point / AU << endl;
//    vdot.pflag = true;
//    double dr = (r1 - vdot.r_cs0) / 50;
//    cout << "| " << r1 << " " << vdot.r_cs0 << " " << dr << endl;
//    for (int ir = 0; ir < 100; ir += 1)
//      vdot(vdot.r_cs0 + ir * dr);
//    vdot.pflag = false;
//    cout << "diter diter_next: " << diter / AU << " " << diter_next / AU << endl;
//    cout << "thetas | rhs:     " << vdot.theta_cs << " " << p_cs.theta() << " | " << rh / AU << " " << p_cs.len() / AU << endl;
//    cout << "r0 rh r1:         " << vdot.r_cs0 / AU << " " << rh / AU << " " << r1 / AU << endl;
//    cout << "vdots:            " << vdot0 << " " << vdot(rh) << " " << vdot(r1) << endl;
//
//    p_cs.set_spherical(vdot.r_cs0, vdot.theta_cs, vdot.phi0 -  fabs(vdot.r_cs0 * ov));
//    cout << "dist: r0 r1 rh    " << (target_point - p_cs).len() / AU;
//    p_cs.set_spherical(r1, vdot.theta_cs, vdot.phi0 -  fabs(r1 * ov));
//    cout << " " << (target_point - p_cs).len() / AU;
//    p_cs.set_spherical(rh, vdot.theta_cs, vdot.phi0 -  fabs(rh * ov));
//    cout << " " << (target_point - p_cs).len() / AU << endl;
//
//    //double rhh = ridders_method(vdot, vdot.r_cs0, rh, 1e-3);
//    //cout << "rhh: " << rhh / AU << " " << vdot(rhh) << endl;
//  }
  assert(diter_next < diter * (1 + 1e-8) && "the spiral_iterate should decrease the distance to the target point");
  diter = diter_next;
  return rh > 0 ? true : false; // if the best fit point goes to negative radius, we should set spiral_iterate to fail and let the system change to point_iterate.
}

class WaveVdot {
  public:
  const Vec p_cs0, target_point;
  const double r_cs0, theta_cs0;
  double ov, phi_cs, cphi, sphi;
  const HCS* h;
  bool pflag = false;

  WaveVdot(double ov_, const Vec& p_cs, const Vec& target_point_, const HCS* h_) :
    ov(ov_), p_cs0(p_cs), target_point(target_point_), r_cs0(p_cs.len()), theta_cs0(p_cs.theta()), phi_cs(p_cs.phi()), h(h_)
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
    if (r < 0) dtheta *= -1;

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

    double phi0 = phi_cs + fabs(r * ov);

    Vec dv;
    if (HCS::hcsform == HCS::Jokipii_Thomas) dv = - sin(HCS::angle) * cos(phi0) / stheta * ov * dtheta + dr;
    else if (HCS::hcsform == HCS::Kota_Jokipii) dv = - tan(HCS::angle) * cos(phi0) * stheta * stheta * ov * dtheta + dr;
    else assert(false && "Unsuported hcsform");

    dv.normalize();
    return dv;
  }

  double operator()(double r, double *vdot_r = NULL) const {
    Vec p_cs;
    if (r == r_cs0) p_cs = p_cs0;
    else
      p_cs.set_spherical(r, h->Theta_S(fabs(r), phi_cs), phi_cs);

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


    if (pflag) {
      cout << "--WVdot " << p_cs / AU << " " << vtangent / 50 << endl;
    }
    return dl.dot(vtangent) / dl.len();
  }
};

bool HCS::wave_iterate(const Vec& target_point, Vec& p_cs, double& diter) const {
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

  if (fabs(r1 - vdot.r_cs0) / vdot.r_cs0 < 1e-9) return false;

  double rh = ridders_method(vdot, vdot.r_cs0, r1, 1e-5);
  p_cs.set_spherical(rh, Theta_S(fabs(rh), vdot.phi_cs), vdot.phi_cs);

  double diter_next = (target_point - p_cs).len();
  if (diter < diter_next && diter_next < diter + resolution)
    return false;

  if (diter < diter_next && (diter_next - diter) / diter < 1e-4)
    return false;
  if (diter < diter_next) {
    cout << endl;
    cout << target_point / AU << endl;
    vdot.pflag = true;
    double dr = (r1 - vdot.r_cs0) / 50;
    cout << "| " << r1 << " " << vdot.r_cs0 << " " << dr << endl;
    for (double vr = vdot.r_cs0 - 50 * dr; fabs(vr - r1) / fabs(r1) > 1e-5; vr += dr)
      vdot(vr);
    vdot.pflag = false;
    cout << diter / AU << " " << diter_next / AU
      << " " << vdot.phi_cs << " " << p_cs.phi() << " " << rh / AU << " " << p_cs.len() / AU << endl;
    cout << vdot.r_cs0 / AU << " " << rh / AU << " " << r1 / AU << " | " << vdot0 << " " << vdot(rh) << " " << vdot(r1) << endl;
    vdot.pflag = true;
    p_cs.set_spherical(vdot.r_cs0, Theta_S(fabs(vdot.r_cs0), vdot.phi_cs), vdot.phi_cs);
    cout << (target_point - p_cs).len() / AU << endl;
    p_cs.set_spherical(r1, Theta_S(fabs(r1), vdot.phi_cs), vdot.phi_cs);
    cout << (target_point - p_cs).len() / AU << endl;
    p_cs.set_spherical(rh, Theta_S(fabs(rh), vdot.phi_cs), vdot.phi_cs);
    cout << (target_point - p_cs).len() / AU << endl;
  }
  assert(diter_next <= diter * (1 + 1e-8) && "the wave_iterate should decrease the distance to the target point");
  diter = diter_next;
  return true;
}

Vec HCS::norm_vec(const Vec& p_cs) const {
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

bool HCS::point_iterate(const Vec& target_point, Vec& p_cs, Vec& dh, double& diter) const {
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

double HCS::get_distance_polygon(double r, double theta, double phi, double Rg2, Polygon polygon) const {
    static map<Polygon, vector<vector<double>>> polygon_coordinates;
    if (polygon_coordinates.size() == 0) {
        polygon_coordinates[Polygon::Icosahedron] = {
            {-8.94427401e-01,  0.00000000e+00},
            { 0.00000000e+00,  0.00000000e+00},
            { 7.23606342e-01, -5.25731196e-01},
            { 7.23606342e-01,  5.25731196e-01},
            {-2.76392839e-01, -8.50650862e-01},
            {-2.76392839e-01,  8.50650862e-01}
        };
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
        polygon_coordinates[Polygon::pseudorandom] = {
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

    auto distance = [&](double r, double phi) {
      double vx, vy, vz;
      info.h->rphi(r, phi, vx, vy, vz);
      return sqrt((vx - info.x) * (vx - info.x) + (vy - info.y) * (vy - info.y) + (vz - info.z) * (vz - info.z));
    };

    double d_min = 1.2 * Rg2;
    for (auto pc : polygon_coordinates[polygon]) {
        double dphi = d_min / r;
        if (dphi > 0.5) {
            dphi = asin(dphi);
        }
        double next_d = distance(r + d_min * pc[0], phi + dphi * pc[1]);
        if (next_d < d_min) d_min = next_d;
    }
    return d_min;
}

inline void show_line(const string& label) {
  //cout << "-------------------------" << label << "-------------------------" << endl;
}
inline void show_log(const string& title, const Vec& target, double diter, double resolution, const Vec& point) {
  //cout << setprecision(13) << title << ": " << target.len() / AU << " " << target.theta() / deg << " " << target.phi() / deg << " | " << point.len() / AU << " " << point.theta() / deg << " " << point.phi() / deg << " -> " << diter / AU << " +- " << resolution / AU << endl;
}

double HCS::get_distance_from_point(const Vec& target, Vec& point) const {
  double diter = 1e5 * AU,
         diter_last = 1e5 * AU;

  point.set_spherical(point.len(), Theta_S(point.len(), point.phi()), point.phi()); // initialize the point to HCS.
  Vec dh = norm_vec(point);
  bool spiral_available = true;
  while (diter == 1e5 * AU || diter_last == 1e5 * AU
         || //(fabs(diter_last - diter) / diter_last > 1e-4 &&
         fabs(diter_last - diter) > resolution) {
    diter_last = diter;
    if (fabs(pi / 2 - target.theta()) > angle - 0.5 * deg && fabs(pi / 2 - point.theta()) > angle - 0.5 * deg && spiral_available) {
      spiral_available = spiral_iterate(target, point, diter);
      show_log("diter_s: ", target, diter, resolution, point);
      if (wave_iterate(target, point, diter) == false) break;
      show_log("diter_w: ", target, diter, resolution, point);
      dh = norm_vec(point);
    } else {
      point_iterate(target, point, dh, diter);
      show_log("diter: ", target, diter, resolution, point);
    }
    //cout << diter_last / AU << " " << diter / AU << " " << (diter_last - diter) / diter_last << endl;
  }
  return (target - point).len();
}

std::mutex mtx;
double HCS::get_distance_intp(double r, double theta, double phi) {
  int i = upper_bound(angle_axis.begin(), angle_axis.end(), HCS::angle) - angle_axis.begin() - 1;
  assert(0 <= i && i < angle_axis.size() - 2 && "The HCS::angle should be in the range of angle_axis.");

  mtx.lock();
  if (!tables[i])
    tables[i] = new KDInterp(table_names[i]);
  mtx.unlock();

  return hcs_interp_eval(HCS::angle, r, theta, phi, tables[i], *this);
}

double HCS::get_distance(double r, double theta, double phi) {
  if (interp) return get_distance_intp(r, theta, phi);

  Vec p_cs;
  return get_distance(r, theta, phi, p_cs);
}

double HCS::get_distance(double r, double theta, double phi, Vec& p_cs) const {
  Vec target;
  target.set_spherical(r, theta, phi);

  double rlow, rup;
  double phi0 = Phi0_S(theta);
  r_bound(r, phi, phi0, rlow, rup);
  if (rlow < 0) rlow = 1e-2*AU; // Avoid negative radius

  double phi_cs = phi, theta_cs = theta;
  if (theta_cs < pi / 2 - angle) theta_cs = pi / 2 - angle;
  else if (theta_cs > pi / 2 + angle) theta_cs = pi / 2 + angle;

  Vec p_cs_l, p_cs_m, p_cs_u;
  show_line("low");
  double dlow = 1e5 * AU;
  if ((rup - r) / (r - rlow) > 0.1) {
    p_cs_l.set_spherical(rlow, theta_cs, phi_cs);
    dlow = get_distance_from_point(target, p_cs_l);
  }

  show_line("mid");
  double dmid = 1e5 * AU;
  if (fabs(pi / 2 - theta) < angle) {
    p_cs_m.set_spherical(r, Theta_S(r, phi), phi);
    dmid = get_distance_from_point(target, p_cs_m);
  }

  show_line("up");
  double dup = 1e5 * AU;
  if ((r - rlow) / (rup - r) > 0.1) {
    p_cs_u.set_spherical(rup, theta_cs, phi_cs);
    dup = get_distance_from_point(target, p_cs_u);
  }

  double d = dlow;
  p_cs = p_cs_l;

  if (dmid < d) {
    d = dmid;
    p_cs = p_cs_m;
  }

  if (dup < d) {
    d = dup;
    p_cs = p_cs_u;
  }

  return sign(r, theta, phi) * d;
}

double HCS::get_raw_distance(double r, double theta) const {
  if (fabs(pi / 2 - theta) < angle) return 0;

  double dangle = theta < pi / 2 ? pi / 2 - angle - theta : theta - pi / 2 - angle;
  return r * sin(dangle);
}
