#include "particle.h"
#include "nlopt.hpp"

using namespace std;
using namespace Unit;

particle::particle() {}

particle::particle(const map<string, docopt::value>& args) {
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

  Ek = mass;
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

const double particle::Theta_S() {
  double value;
  value = pi / 2. +
          asin(sin(angle) * sin(phi + Omega * r / (400 * km / sec)) / 0.8354);
  // cout << "theta_s :  " << value << "   " << sin(Omega*r/(400*km/sec))) <<
  // endl; getchar();
  return value;
}

const double particle::Heav() {
  double theta_s = Theta_S();
  double value;
  if (theta < theta_s)
    value = 1.;
  else if (theta_s < theta)
    value = -1.;

  return value;
}

const double particle::B_r(const double& heaviside) {
  const double r0 = 1 * AU;
  // std::cout << "D:   " << heaviside << "   " << polarity << "   " << r << "
  // " << r0 << std::endl;
  //     getchar();
  return B0 * heaviside * polarity / pow(r0 / r, 2.);
}

const double particle::B_p(const double& heaviside) {
  // std::cout << "D:   " << Omega << "   " << sin(theta) << "   " << Vs << " "
  // << polarity<< std::endl; getchar();
  return -1. * B0 * r * Omega * sin(theta) * heaviside * polarity / Vs;
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

  return ky;
}

double particle::HCS_rphi_z(const double &r, const double &phi) const {
  double phi0 = phi + r * Omega / Vs_eq - Omega * (t - t0);
  double alpha = angle;
  double tan_theta_cs = tan(alpha) * sin(phi0);

  return r * tan_theta_cs;
}

double particle::HCS_xy_z(const double &x, const double &y) const {
  double r = sqrt(x * x + y * y);
  double phi = atan2(y, x);
  return HCS_rphi_z(r, phi);
}

struct nlopt_info {
  double x, y, z;
  const particle *p;
};
double distance_to_point(const std::vector<double> &x, std::vector<double> &grad, void *voidp) {
  nlopt_info *info = reinterpret_cast<nlopt_info*>(voidp);
  double z = info->p->HCS_xy_z(x[0], x[1]);

  double dx = x[0] - info->x;
  double dy = x[1] - info->y;
  double dz = z - info->z;
  return sqrt(dx * dx + dy * dy + dz * dz);
}

double particle::get_HCS_distance() const {
  nlopt_info info = { r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta), this };

  nlopt::opt opt(nlopt::LN_SBPLX, 2);
  opt.set_min_objective(distance_to_point, &info);

  double dr = pi / Omega * Vs_eq;
  vector<double> position_hcs_low = { info.x - 2 * dr * cos(phi), info.y - 2 * dr * sin(phi) };
  vector<double> position_hcs_up = { info.x + 2 * dr * cos(phi), info.y + 2 * dr * sin(phi) };
  opt.set_lower_bounds(position_hcs_low);
  opt.set_upper_bounds(position_hcs_up);

  double dlow, dup;
  vector<double> position_hcs = { info.x - dr * cos(phi) , info.y - dr * sin(phi) };
  opt.optimize(position_hcs, dlow);

  position_hcs = { info.x + dr * cos(phi), info.y + dr * sin(phi) };
  opt.optimize(position_hcs, dup);

  return fmin(dlow, dup);
}

void particle::step() {
  random_device rd;
  mt19937 gen(rd());
  double mean = 0.0;
  double dev = 1.0;
  normal_distribution<double> dist(mean, dev);

  double record_T = 0.;

  while (r < boundary || record_T < pow(10., 10.)) {
    record_T += dt;

    M_p = sqrt(Ek * (Ek + 2. * mass));
    rigidity = A / Z * M_p;
    V_p = M_p / (Ek + mass) * light;
    Vs = Wind();
    heaviside = Heav();
    Br = B_r(heaviside);
    Bp = B_p(heaviside);
    psi = atan(fabs(Bp / Br));
    // std::cout << "D:   " << M_p/GeV << "   " << V_p/AU << "   " << Vs/1000.
    // << "   " << Ek/GeV << "   " << mass/GeV << std::endl; getchar();

    double B = sqrt(Br * Br + Bp * Bp);

    //  std::cout << "D:   " << heaviside << "   " << Br << "   " << Bp << "   "
    //  << B << std::endl;
    // getchar();

    k_rr = K_rr();
    k_tt = K_tt();
    k_pp = K_pp();
    //          std::cout << "D:   " << k_rr << "   " << k_pp << "   " << k_tt
    //          << std::endl;
    // getchar();

    double r0 = 1.0 * AU;
    double gamma = r * Omega * sin(theta) / Vs;
    double drift = 2 * M_p * V_p * r / (3 * Z * e * B0 * r0 * r0 * light);

    Vdr_gc = drift / pow(1 + gamma * gamma, 2.) * heaviside *
             (-1. * gamma / tan(theta));
    Vdp_gc = drift / pow(1 + gamma * gamma, 2.) * heaviside *
             (2. + gamma * gamma) * gamma;
    Vdt_gc = drift / pow(1 + gamma * gamma, 2.) * heaviside * gamma * gamma /
             tan(theta);
    // cout << "Z e B0: " << Z << " " <<  e << " " << B0 << " " << r / AU << " "
    // << r0 / AU << " " << light << " " << k_rr / (cm * cm / sec) << " " <<
    // (-1.*Vs - Vdr_gc - Vdr_HCS) * dt / AU << " " << sqrt(k_rr * dt) / (AU) <<
    // endl; cout << "momentum : " << M_p / GeV << " " << V_p / (km/sec) << " "
    // << Z << " " << Z << " " << A << endl;
    //  std::cout << "VD :  " << gamma << "  " << Vdr_gc << "  " << Vdp_gc << "
    //  " << Vdt_gc << "   " << drift << std::endl; std::cout << "VD :  " <<
    //  Vdp_gc << "  " << drift << "  " << pow(1+gamma*gamma, 2.) << "  " <<
    //  heaviside << "   " << (2. + gamma*gamma)
    //  << "   "  << std::endl;
    //  getchar();
    double Vd = drift / (1 + gamma * gamma);

    double L0, Rg;

    double dw = dist(gen);
    double d_HCS = get_HCS_distance();

    // r += (-1.*Vs - Vdr_gc - Vdr_HCS + 2./r) * dt
    //  cout << "rbefore: " << r / AU << endl;
    double ro = r;
    double theta0 = theta;
    double phi0 = phi;
    r += (-1. * Vs - Vdr_gc - Vdr_HCS + 2. * k_rr / r) * dt +
         sqrt(2. * k_rr * dt) * dw;
    // std::cout << "r:   " << Vs << "   " << Vdr_gc<< "   " << 2.*k_rr/r
    // << "   " << sqrt(2. * k_rr * dt) << "   " << k_rr << std::endl;
    // getchar();
    // cout << "rafter: " << r / AU << " " << (-1.*Vs - Vdr_gc - Vdr_HCS) * dt /
    // AU << " " << sqrt(2. * k_rr * dt) / AU << " " << dw << endl;

    theta += (-1. * Vdt_gc / r - Vdt_HCS / r +
              1. / (r * r * sin(theta)) * cos(theta) * k_tt) *
                 dt +
             1. / r * sqrt(2. * k_tt * dt) * dw;
    // std::cout << "theta:   " << Vdt_gc << "   "
    // << 1./(r*r*sin(theta))*cos(theta)*k_tt*r
    // << "   " << 1./r * sqrt(2.*k_tt*dt) << "   " << k_tt << std::endl;
    // getchar();

    phi += (-1. * Vdp_gc - Vdp_HCS) / (r * sin(theta)) * dt +
           sqrt(2. * k_pp * dt) * dw / (r * sin(theta));
    // std::cout << "phi:   " << (-1.*Vdp_gc - Vdp_HCS) / (r * sin(theta)) * dt
    // << "   " << sqrt(2.*k_pp*dt) * dw / (r * sin(theta)) << "   " << k_pp <<
    // std::endl; getchar();

    // std::cout << "D:   " << r-r0 << "   " << (theta-theta0)*r<< "   " <<
    // (phi-phi0)*r*sin(theta) << std::endl; getchar();

    Ek += 2. * Vs / (3. * r) * (Ek * Ek + 2. * Ek * mass) / (Ek + mass) * dt;

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

    if (phi < 0.)
      phi = 2. * pi - phi;
    else if (2. * pi < phi)
      phi -= 2. * pi;
  }
}