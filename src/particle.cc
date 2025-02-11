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
  polarity(-1), B0(5 * nT), indexA(2), D(5 * 1e22 * cm * cm / sec),
  Bn(B0 * AU * AU / 1.35883),
  r(AU), theta(90*deg + 1e-10), phi(1e-10), hcs(Wind())
  {}

particle::particle(const map<string, docopt::value>& args) :
  A(args.at("--A").asLong()), Z(args.at("--Z").asLong()),
  polarity(args.at("--polarity").asLong()), B0(stod(args.at("--B0").asString()) * nT), indexA(stod(args.at("--indexA").asString())), D(stod(args.at("--D").asString()) * 1e22 * cm * cm / sec),
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

const double particle::Heav() {
  double theta_s = hcs.Theta_S(r, phi);
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
  return -1. * Bn / r * HCS::Omega * sin(theta) * heaviside * polarity / Vs;
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

  auto r2V_r = [&](double r_) {
    return r_ * r_ * Wind(r_, theta, phi, HCS::angle);
  };

  auto dvdx = [&](double x, const function<double(double)>& v) {
    return (v(x * (1 + 1e-3)) - v(x)) / (1e-3 * x);
  };

  // theta = 1e-3;
  double Dt = 0;
  double M_p0 = sqrt(Ek * (Ek + 2. * mp));
  double drift = 0;
  double Vdr_gc = 0, Vdt_gc = 0, Vdp_gc = 0;
  double Vdr_HCS = 0, Vdt_HCS = 0, Vdp_HCS = 0;
  double k_rr1, k_tt1, k_pp1;
  double Rg, Vns;

  ostringstream osname;
  if (!logname.empty()) osname << "s" << seed << "_" << logname;
  std::ofstream logfile(osname.str());

  if (logfile.is_open())
     logfile << "t[month],r[AU],theta[rad],phi[rad],Ek[GeV],drift[km/s],Vdr_gc[km/s],Vdr_HCS[km/s]" << endl;
  auto write_log = [&]() {
    if (logfile.is_open())
      logfile << Dt/day/30. << "," << r/AU << "," << theta << "," << phi << "," << Ek/GeV
        << "," << drift/(km/sec)
        << "," << Vdr_gc/(km/sec) << "," << Vdr_HCS/(km/sec)
        << endl;
  };

  while (r<boundary) {//theta<pi/2
    Dt += dt;
    M_p = sqrt(Ek * (Ek + 2. * mp));
    rigidity = A / (Z * e) * M_p;
    V_p = M_p / (Ek + mp) * c_speed;
    Vs = Wind();

    // std::cout << "Mp:  " << r << "  " << M_p/GeV << "  " << Ek/GeV << "  " << mp/GeV << std::endl;
    // getchar();

    double cs = hcs.Theta_S(r, phi);
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
        
    k_rr1 = K_rr(r*(1+1e-3), theta*(1+1e-3), phi, psi1, B1, B0, M_p, M_p0, V_p);
    k_tt1 = K_tt(r*(1+1e-3), theta*(1+1e-3), phi, psi1, B1, B0, M_p, M_p0, V_p);
    k_pp1 = K_pp(r*(1+1e-3), theta*(1+1e-3), phi, psi1, B1, B0, M_p, M_p0, V_p);

    double dr2V_dr = dvdx(r, r2V_r);

    double gamma = tan(psi);
    drift = 2 * A * M_p * V_p * r / (3 * Z * e * Bn ) * heaviside * polarity;
    Vdr_gc = drift / pow(1 + gamma * gamma, 2.) * (-1. * gamma ) / fabs(tan(theta));
    Vdt_gc = drift / pow(1 + gamma * gamma, 2.) *  (2. + gamma * gamma) * gamma;
    Vdp_gc = drift / pow(1 + gamma * gamma, 2.) * gamma * gamma / fabs(tan(theta));//;

    double delta = (hcs.Theta_S(r+0.1, phi) - hcs.Theta_S(r, phi)) / 0.1;
    if(delta<0) delta = -1.;
    else if(0<(delta)) delta = 1.;
    double beta = atan(HCS::Omega * r * sqrt(fabs(sin(HCS::angle) * sin(HCS::angle) - cos(cs) * cos(cs))) / (Vs * sin(psi) * sin(cs))) * delta;
    if(Z*polarity<0) beta = pi + beta;
    else if(0<Z*polarity) beta = beta;

    Rg = A * M_p / (B * Z * e * c_speed);
    double d_HCS = fabs(hcs.get_distance(r, theta, phi));
    Vns = 0.;
    if(d_HCS<2.*Rg) Vns = (0.457 - 0.412 * d_HCS / Rg + 0.0915 * d_HCS * d_HCS / Rg / Rg) * V_p * A_drift;//

    double zonal = dis(gen);
    Vdr_HCS = Vns * cos(beta) * sin(zonal);
    Vdt_HCS = Vns * sin(beta);
    Vdp_HCS = Vns * cos(beta) * cos(zonal);

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

    write_log();

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
    double V1 = Wind(r1, theta0, phi0, HCS::angle);
    double E = Ek + mp;
    double p2 = E * E - mp * mp;

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

  if (logfile.is_open()) { // Write the final state
    write_log();
    logfile.close();
  }
}
