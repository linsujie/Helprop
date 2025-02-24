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

const double particle::mp = 0.93827 * GeV;
particle::particle() :
  A(1), Z(1),
  polarity(-1), B0(5 * nT), indexA(2), D0(5 * 1e22 * cm * cm / sec), rigidity0(1 * GeV / e),
  Bn(B0 * AU * AU / 1.35883),
  r(AU), theta(90*deg + 1e-10), phi(1e-10), hcs(Wind())
  {}

particle::particle(const map<string, docopt::value>& args) :
  A(args.at("--A").asLong()), Z(args.at("--Z").asLong()),
  polarity(args.at("--polarity").asLong()), B0(stod(args.at("--B0").asString()) * nT), indexA(stod(args.at("--indexA").asString())), D0(stod(args.at("--D0").asString()) * 1e22 * cm * cm / sec), rigidity0(stod(args.at("--R0").asString()) * GeV / e),
  Bn(B0 * AU * AU / 1.35883),
  r(AU), theta(90*deg + 1e-6), phi(1e-10), hcs(Wind())
{}

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
  return Wind(r, theta, phi, HCS::angle);
}

double particle::Heav() {
  theta_s = hcs.Theta_S(r, phi);
  double value;
  if (theta <= theta_s)
    value = 1.;
  else if (theta_s < theta)
    value = -1.;

  return value;
}

double particle::B_r(double r, double theta, double heaviside) const {
  return Bn / r / r * heaviside * polarity;
}

double particle::B_p(double r, double theta, double heaviside) const {
  return - Bn / r * HCS::Omega * sin(theta) * heaviside * polarity / Vs;
}

double particle::Kpara0() const { // parallel diffusion coefficient at the earth
  return D0 * V_p * pow(rigidity / rigidity0, indexA);
}

double f_perp_t(double theta) {
  if (theta > pi / 2) theta = pi - theta;
  return 1.0 - 0.5 * tanh(8 * (theta - pi / 2 + 35 * deg));
}

void particle::K(double r, double theta, double heaviside, double kpara, double& B, double& krr, double& ktt, double& kpp, double& krp) const {
  double Br = B_r(r, theta, heaviside);
  double Bp = B_p(r, theta, heaviside);
  B = sqrt(Br * Br + Bp * Bp);

  kpara *= B0 / B;
  double k_perp_r = 0.02 * kpara;
  double k_perp_t = k_perp_r * f_perp_t(theta);

  double sin_psi = Bp / B,
         cos_psi = Br / B;

  krr = cos_psi * cos_psi * kpara + sin_psi * sin_psi * k_perp_r;
  kpp = cos_psi * cos_psi * k_perp_r + sin_psi * sin_psi * kpara;
  krp = sin_psi * cos_psi * (k_perp_r - kpara);

  ktt = k_perp_t;
}

void particle::K(double r, double theta, double heaviside, double kpara, double& krr, double& ktt, double& kpp, double& krp) const {
  double B;
  K(r, theta, heaviside, kpara, B, krr, ktt, kpp, krp);
}

void particle::coord_trans(double krr, double ktt, double kpp, double krp, double& dwr, double& dwt, double& dwp) const {
  double mrr = sqrt(2 * krr - 2 * krp * krp / kpp);
  double mrp = sqrt(2) * krp / sqrt(kpp);
  double mtt = sqrt(2 * ktt) / r;
  double mpp = sqrt(2 * kpp) / r / sin(theta);

  dwr = mrr * dwr + mrp * dwp;
  dwt = mtt * dwt;
  dwp = mpp * dwp;
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

  double record_T = 0.;

  // theta = 1e-3;
  const double dh = 1e-3;
  const double h = 1 + dh;
  double Dt = 0;
  double drift = 0;
  double Vdr_gc = 0, Vdt_gc = 0, Vdp_gc = 0;
  double Vkrr = 0, Vkrp = 0;
  double Vdr_HCS = 0, Vdt_HCS = 0, Vdp_HCS = 0;
  double Rg, Vns, d_HCS;
  double dr,dtheta,dphi,dEk,dwr,dwt,dwp;
  double Vs_dr, dr2V_dr;
  double krr, ktt, kpp, krp;
  double krr_dr, ktt_dr, kpp_dr, krp_dr;
  double krr_dt, ktt_dt, kpp_dt, krp_dt;

  ostringstream osname;
  if (!logname.empty()) osname << "s" << seed << "_" << logname;
  std::ofstream logfile(osname.str());

  if (logfile.is_open())
     logfile << "t[month],r[AU],theta[rad],phi[rad],Ek[GeV],dEk[GeV],Vs[km/s],drift[km/s],Vdr_gc[km/s],Vdt_gc[km/s],Vdp_gc[km/s],Vdr_HCS[km/s],Vdt_HCS[km/s],Vdp_HCS[km/s],d_HCS[AU],Rg[AU],dwr[AU],dwt[rad],dwp[rad],dr[AU],dtheta[rad],dphi[rad]" << endl;
  auto write_log = [&]() {
    if (logfile.is_open())
      logfile << Dt/day/30. << "," << r/AU << "," << theta << "," << phi << "," << Ek/GeV << "," << dEk/GeV
        << "," << Vs / (km/sec)
        << "," << drift/(km/sec)
        << "," << Vdr_gc/(km/sec) << "," << Vdt_gc/(km/sec) << "," << Vdp_gc/(km/sec)
        << "," << Vdr_HCS/(km/sec) << "," << Vdt_HCS/(km/sec) << "," << Vdp_HCS/(km/sec)
        << "," << d_HCS/AU << "," << Rg/AU
        << "," << dwr/AU << "," << dwt << "," << dwp
        << "," << dr/AU << "," << dtheta << "," << dphi
        << endl;
  };

  while (r<boundary) {//theta<pi/2
    double E = Ek + mp;
    double p2 = E * E - mp * mp;
    M_p = sqrt(p2);
    rigidity = A / (Z * e) * M_p;
    V_p = M_p / E * c_speed;
    Vs = Wind();

    // std::cout << "Mp:  " << r << "  " << M_p/GeV << "  " << Ek/GeV << "  " << mp/GeV << std::endl;
    // getchar();

    heaviside = Heav();
    double kpara = Kpara0();

    //dt = kpara / (c_speed * c_speed);
    Dt += dt;

    double B = 0;
    K(r, theta, heaviside, kpara, B, krr, ktt, kpp, krp);
    K(r * h, theta, heaviside, kpara, krr_dr, ktt_dr, kpp_dr, krp_dr);
    K(r, theta * h, heaviside, kpara, krr_dt, ktt_dt, kpp_dt, krp_dt);
    Vs_dr = Wind(r * h, theta, phi, HCS::angle);

    double dr2krr_dr = (r * r * h * h * krr_dr - r * r * krr) / r / dh;
    double dstktt_dt = (sin(theta * h) * ktt_dt - sin(theta) * ktt) / theta / dh;
    double dkrp_dp = 0,
           dkpp_dp = 0;
    double drkrp_dr = (r * h * krp_dr - r * krp) / r / dh;

    double gamma = r * HCS::Omega * sin(theta) / Vs;

    //cout << M_p << " " << V_p << " " << r / AU << " " << B0 << " " << Bn << " " << Z << " " << e << endl;
    drift = 2 * A * M_p * V_p * r / (3 * Z * e * Bn) * heaviside * polarity;
    Vdr_gc = drift / pow(1 + gamma * gamma, 2.) * (- gamma ) / tan(theta);
    Vdt_gc = drift / pow(1 + gamma * gamma, 2.) *  (2 + gamma * gamma) * gamma;
    Vdp_gc = drift / pow(1 + gamma * gamma, 2.) * gamma * gamma / tan(theta);

    Vns = 0.;
    Rg = fabs(A * M_p / (B * Z * e * c_speed));
    d_HCS = hcs.get_raw_distance(r, theta);
    if (d_HCS < 2 * Rg) {
      d_HCS = fabs(hcs.get_distance(r, theta, phi));
      if(d_HCS < 2*Rg)
        Vns = (0.457 - 0.412 * d_HCS / Rg + 0.0915 * d_HCS * d_HCS / Rg / Rg) * V_p * polarity * A_drift;//
  
      double Vrx = r * sin(theta_s) * HCS::Omega;
      double Vtx = - tan(HCS::angle) * sin(theta_s) * cos(hcs.phi0(r, phi)) * (hcs.Vs_eq * hcs.Vs_eq + Vrx * Vrx) / hcs.Vs_eq;
      double Vpx = hcs.Vs_eq;
      double Vtot = sqrt(Vrx * Vrx + Vtx * Vtx + Vpx * Vpx);
  
      Vdr_HCS = Vrx / Vtot * Vns;
      Vdt_HCS = Vtx / Vtot * Vns;
      Vdp_HCS = Vpx / Vtot * Vns;
    } else
      Vdr_HCS = Vdt_HCS = Vdp_HCS = 0;

    double Vdr = Vdr_gc + Vdr_HCS;
    double Vdt = Vdt_gc + Vdt_HCS;
    double Vdp = Vdp_gc + Vdp_HCS;

    dwr = dist(gen) * sqrt(dt);
    dwt = dist(gen) * sqrt(dt);
    dwp = dist(gen) * sqrt(dt);

    coord_trans(krr, ktt, kpp, krp, dwr, dwt, dwp);
    //if (r < 4 * AU) dwr = fabs(dwr);

    dr = - (Vs + Vdr
            - 1 / r / r * dr2krr_dr
            - 1 / r / sin(theta) * dkrp_dp
            ) * dt
            + dwr;

    dtheta = - (Vdt
                - 1 / r / sin(theta) * dstktt_dt
                ) / r * dt
                + dwt;

    dphi = - (Vdp
              - 1 / r / sin(theta) * dkpp_dp
              - 1 / r * drkrp_dr
              ) / (r * sin(theta)) * dt
              + dwp;

    dEk = 0;
    if (r >= 1 * AU && r + dr >= 1 * AU) {
      double r_next = r + dr;
      Vs_dr = Wind(r_next, theta, phi, HCS::angle);
      dr2V_dr = (r_next * r_next * Vs_dr - r * r * Vs) / dr;
      dEk = dr2V_dr / 3 / fmax(r, fabs(r_next)) / fmax(r, fabs(r_next)) * p2 / E * dt;
    }

    write_log();

    r += dr;
    if (r < 0.5 * AU) r = 1 * AU - r; // reflect the particle from the near center region to avoid them captured

    theta += dtheta;
    phi += dphi;
    Ek += dEk;

    Ek = fabs(Ek);

    if (r < 0) {
      r = -r;
      theta = pi - theta;
      phi = - phi;
    }

    if (theta < - pi || pi < theta)
      theta -= floor(theta / (2 * pi) + 0.5) * 2 * pi;
    if (theta < 0.) {
      theta = fabs(theta);
      phi += pi;
    }
    theta = fmin(fmax(theta, 1e-3), pi - 1e-3);

    if (phi < 0 || 2 * pi < phi)
      phi -= floor(phi / (2 * pi)) * 2 * pi;
  }

  if (logfile.is_open()) { // Write the final state
    write_log();
    logfile.close();
  }
}
