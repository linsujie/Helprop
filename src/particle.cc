#include "particle.h"
#include "nlopt.hpp"

using namespace std;
using namespace Unit;

particle::particle() {}

particle::particle(const map<string, docopt::value>& args) : hcsform(Jokipii_Thomas) {
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

//   Ek = mass;
//   std::cout << "Mass:  " << mass/GeV << "   " << Ek/GeV << std::endl;
//   getchar();
}

particle::~particle() {}

const double particle::Wind() {
  double value;
  if (0. < theta && theta < pi / 2.) {
    value = 1.475 - 0.4 * tanh(6.8 * ((theta - pi / 2) + (15 * deg + angle))) *
                        (3.5 / 5. - 1.5 / 5. * tanh((r - 95 * AU) / 1.2 / AU));
  } else if (pi / 2. < theta && theta < pi) {
    value = 1.475 + 0.4 * tanh(6.8 * ((theta - pi / 2) - (15 * deg + angle))) *
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

double particle::Theta_S_Jokipii_Thomas(double r, double phi) const {
  return pi / 2 - asin(sin(angle) * sin(phi + Omega * r / Vs_eq));
}
double particle::Theta_S_Kota_Jokipii(double r, double phi) const {
  double phi0 = phi + r * Omega / Vs_eq - Omega * (t - t0);
  double alpha = angle;
  return pi / 2 - atan(tan(alpha) * sin(phi0));
}
double particle::Theta_S(double r, double phi) const {
  if (hcsform == Jokipii_Thomas)
    return Theta_S_Jokipii_Thomas(r, phi);
  else if (hcsform == Kota_Jokipii)
    return Theta_S_Kota_Jokipii(r, phi);

  // return 0;
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

const double particle::B_r(const double& heaviside) {
  // double r0 = 1 * AU;
  // std::cout << "D:   " << heaviside << "   " << polarity << "   " << r << "
  // " << r0 << std::endl;
  //     getchar();

  return Bn * heaviside * polarity * pow(1. / r, 2.);
}

const double particle::B_p(const double& heaviside) {
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

const double particle::K_tt() {
  double kx = D * pow(Ek / GeV, indexA);
  double ky = 0.02 * kx;
  double kz;
  if (theta < pi / 2. && 0 < theta)
    kz = ky * (2. - 1. * tanh(8 * ((theta + (-90. + 35.) * deg))));
  else if (pi / 2. < theta && theta < pi)
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
  const particle *p;
};
double distance_to_point(const std::vector<double> &x, std::vector<double> &grad, void *voidp) {
  nlopt_info *info = reinterpret_cast<nlopt_info*>(voidp);

  auto distance = [&](double r, double phi) {
    // r = 53.01 * AU;
    // phi = 1e-10;
    // std::cout << "Find:  " << r / AU << "  " << phi ;
    // getchar();
    double vx, vy, vz;
    info->p->HCS_rphi(r, phi, vx, vy, vz);
    // std::cout << r / AU << " " << phi << " "
    //  << info->p->Theta_S(r, phi) << " "
    //  << "  " << sqrt((vx - info->x) * (vx - info->x) + (vy - info->y) * (vy - info->y) + (vz - info->z) * (vz - info->z)) << std::endl;
    // std:: cout << "  " << info->x << "  " << info->y << "  " << info->z << std::endl;
    // std::cout << vx << "  " << vy << "  " << vz << std::endl;
    // std::cout << "Test:  " << r << "  " << phi  << endl;
    return sqrt((vx - info->x) * (vx - info->x) + (vy - info->y) * (vy - info->y) + (vz - info->z) * (vz - info->z));
  };


  double d = distance(x[0], x[1]);
  double ddr = distance(x[0] + 1e-10 * AU, x[1]);
  double ddp = distance(x[0], x[1] + 1e-10);

  // std::cout << "Test:  " << x[0] << "  " << x[1] << "  " << x[2] << endl;

  grad = { (ddr - d) / (1e-10 * AU), (ddp - d) / (1e-10) };

  return d;
}

double particle::get_HCS_distance() const {
  nlopt_info info = { r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta), this };
  // std::cout << "befor:  " << r/AU << "  " << theta << "  " << phi << "  " << r * sin(theta) * cos(phi) << "  " << r * sin(theta) * sin(phi) << "  " << r * cos(theta) << std::endl;
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
  double theta_hcs = Theta_S(r, phi);
  assign_region(-dr);
  opt.optimize(position_hcs, dlow);

  assign_region(0);
  opt.optimize(position_hcs, dmid);

  assign_region(dr);
  opt.optimize(position_hcs, dup);

  double vx, vy, vz;
  HCS_rphi(r, phi, vx, vy, vz);
  double sign = info.z > vz ? 1 : -1;

  return sign * fmin(fmin(dlow, dup), dmid);
}

void particle::step() {
    std::ofstream f2;
    f2.open("file2");

  random_device rd;
  mt19937 gen(rd());
  double mean = 0.0;
  double dev = 1.0;
  normal_distribution<double> dist(mean, dev);

  double record_T = 0.;
  hcsform = Kota_Jokipii;

  // theta = 1e-3;
  while (r<boundary) {

    M_p = sqrt(Ek * (Ek + 2. * mass));
    rigidity = A / (Z * e) * M_p;
    V_p = M_p / (Ek + mass) * light;
    Vs = Wind();
    
    double cs = Theta_S(r, phi);
    heaviside = Heav();
    Br = B_r(heaviside);
    Bp = B_p(heaviside);
    psi = atan(fabs(Bp / Br));

    double B = sqrt(Br * Br + Bp * Bp);

    k_rr = K_rr();
    k_tt = K_tt();
    k_pp = K_pp();

    double r0 = 1.0 * AU;
   
    double gamma = tan(psi);//r * Omega * sin(theta) / Vs;

    double drift = 2 * M_p * V_p * r / (3 * Z * e * light * Bn );//B0 * r0 * r0

    Vdr_gc = -1. * drift / pow(1 + gamma * gamma, 2.) * (-1. * gamma );/// tan(theta)) * polarity;
    Vdt_gc = -1. * drift / pow(1 + gamma * gamma, 2.) *  (2. + gamma * gamma) * gamma * heaviside;
    Vdp_gc = drift / pow(1 + gamma * gamma, 2.) * gamma * gamma ;/// tan(theta);

    double beta = atan(Omega * r * sqrt(fabs(sin(angle) * sin(angle) - cos(cs) * cos(cs))) / (Vs * sin(psi) * sin(cs)));

    double Rg = M_p / (B * Z * e * light);
std::cout << "begin" << r/AU << "  " << theta << "  " << pi/2. << std::endl;
    double d_HCS = fabs(get_HCS_distance());
std::cout << "begin" << std::endl;
    double Vns = 0.;
    if(d_HCS<2.*Rg) Vns = (0.457 - 0.412 * d_HCS / Rg + 0.0915 * d_HCS * d_HCS / Rg / Rg) * V_p;//
    
    Vdr_HCS = -1. * Vns * cos(beta) * sin(psi) * Z * polarity;
    Vdp_HCS = -1. * Vns * cos(beta) * cos(psi) * Z * polarity;
    Vdt_HCS = -1. * Vns * sin(beta) * Z * polarity;

    r0 = r;
    double theta0 = theta;
    double phi0 = phi;
    
    double dwr = dist(gen);
    double dwp = dist(gen);
    double dwt = dist(gen);
    if(3<fabs(dwr)) dwr = dist(gen);
    if(3<fabs(dwp)) dwp = dist(gen);
    if(3<fabs(dwt)) dwt = dist(gen);

    r += (-1. * Vs - Vdr_gc - Vdr_HCS + 2. * k_rr / r) * dt +
         sqrt(2. * k_rr * dt) * dwr;

    theta += (-1. * Vdt_gc / r - Vdt_HCS / r +
              1. / (r * r * sin(theta)) * cos(theta) * k_tt) *
                 dt +
             1. / r * sqrt(2. * k_tt * dt) * dwt;

    phi += (-1. * Vdp_gc - Vdp_HCS) / (r * sin(theta)) * dt +
           sqrt(2. * k_pp * dt) * dwp / (r * sin(theta));

    M_p += 2. * Vs / (3. * r) * M_p * M_p / (Ek + mass)  * dt;
    Ek = sqrt(M_p * M_p + mass * mass) - mass;
    if (r < 0.) {
      r = 0.;
      break;
    }

    if (theta < 0.) {
      theta = fabs(theta);
      phi += pi;
    } else if (pi < theta) {
      theta = 2. * pi - theta;
      phi += pi;
    }

    if (phi < 0.) phi = 2. * pi - phi;
    else if (2. * pi < phi) phi -= 2. * pi;

    f2 << r/AU << "   " << theta << "   " << phi << std::endl;
    // if(60*60*24*365*1.5<record_T) break;
    std::cout << r/AU << "  " << theta << std::endl;
  }
  std::cout << "getOne." << std::endl;
  f2.close();
  getchar();
}