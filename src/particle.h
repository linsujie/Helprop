#ifndef particle_h_
#define particle_h_

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include "Vec.hh"
#include "fcache.h"
#include "docopt.h"

namespace Unit {
  const double m = 1;
  const double cm = 1e-2 * m;
  const double km = 1e3 * m;
  const double AU = 1.496*pow(10.,11.) * m;

  const double sec = 1;
  const double min = 60 * sec;
  const double hr = 3600 * sec;
  const double day = 86400 * sec;

  const double c_speed = 2.99792e8 * m / sec;

  const double kg = 1;
  const double g = 1e-3 * kg;
  const double ton = 1e3 * kg;

  const double J = 1 * kg * m * m / sec / sec;

  const double pi = acos(-1);
  const double deg = pi / 180;

  const double T = 1;
  const double nT = 1e-9 * T;
  const double Gauss = 1e-4 * T;

  const double GeV = 1.602177e-10 * J;
  const double MeV = 1e-3 * GeV;
  const double TeV = 1e3 * GeV;

  const double hbar = 1.05457e-34 * J * sec;

  const double C = 1 * J * sec / T / m / m;
  const double e = 1.602e-19 * C;
  const double V = J / C;
  const double epsilon_0 = 8.854e-12 * C / (V * m); //permittivity of free space in F/m
};
class particle {
    public:

    long seed = 0;
    long fix_seed = false;
    double A_drift = 0.5;

    double boundary = 100 * Unit::AU;                                     //boundary condition in AU
    double dt = 100 * Unit::sec;                                             //time interval per step
    double Vs;                                                  //solar wind velocity
    double Vs_eq;
    double Omega = 2*Unit::pi/27.5/Unit::day;                        //angular velocity corresponding to 27.5 day
    double t0 = 0.0;
    double t = 0.0;

    double polarity;                                            //field direction
    static double angle;                                               //tilt angle of HCS
    double B0;                                                  //magnetic strength in the Earth in T

    double k_xx;                                                //field along diffusion coefficient in local fram
    double k_yy;                                                //
    double k_zz;                                                //azimuth diffusion coefficient in local frame
    double indexA;                                              //power index related to rigidity of particle
    double D;                                                   //diffusion factor
    double k_rr;                                                //radial diffusion coefficient
    double k_tt;                                                //pole angle diffusion coefficient
    double k_pp;                                                //azimuthal diffusion coefficient

    double Br;                                                  //magnetic field components in three direction
    double Bp;                                                  //
    double Bn;
    double Bt = 0.0;                                            //
    double psi;                                                 //pitch angle between field and radial direction
    double theta_s;                                             //tile angle of HCS at particle point
    double heaviside;                                           //field direction in particle point

    double rigidity;                                            //rigidity of particle related to kinetic and rest energy
    double A = 1.;                                              //nucleon number, proton by default
    double Z = 1.;                                              //charge number, proton by default
    double mass;                                                  //rest mass

    double Vdr_gc;                                              //drift to radial direction
    double Vdp_gc;                                              //drift to azimuthal direction
    double Vdt_gc;                                              //drift to pole direction
    double Vdr_HCS = 0;                                             //drift at the HCS for thress direction
    double Vdp_HCS = 0;
    double Vdt_HCS = 0;

    inline double phi0(double r, double phi) const {
      return phi + r * Omega / Vs_eq - Omega * (t - t0);
    }
    void r_bound(double r, double phi, double phi0, double& rlow, double& rup) const;

    void spiral_iterate(const Vec& target_point, Vec& p_cs, double& diter) const;
    void wave_iterate(const Vec& target_point, Vec& p_cs, double& diter) const;
    void point_iterate(const Vec& target_point, Vec& p_cs, Vec& dh, double& diter) const;

    Vec norm_vec(const Vec& p_cs) const;

    double Phi0_S_Jokipii_Thomas(double) const;
    double Phi0_S_Kota_Jokipii(double) const;

    double Theta_S_Jokipii_Thomas(double) const;
    double Theta_S_Kota_Jokipii(double) const;

    inline double Theta_S_Jokipii_Thomas(double r, double phi) const {
      return Theta_S_Jokipii_Thomas(phi0(r, phi));
    }
    inline double Theta_S_Kota_Jokipii(double r, double phi) const {
      return Theta_S_Kota_Jokipii(phi0(r, phi));
    }

    public:
    enum HCSFORM { Jokipii_Thomas, Kota_Jokipii };
    static HCSFORM hcsform;

    particle(const std::map<std::string, docopt::value>& args);
    particle();
    ~particle();

    void step(const std::string& logname = "");                                       //simulate trajectory of particle

    const double Wind();                                        //solar wind velocity function
    double Theta_S(double, double) const;                                     //theat_s function
    double Phi0_S(double) const;
    const double Heav();                                        //get heaviside function
    const double B_r(const double& r, const double &heaviside);                  //radial magnetic field function
    const double B_p(const double& r, const double& theta, const double &heaviside);                  //azimuthal magnetic field function
    double get_HCS_distance_old() const;
    double get_HCS_distance() const;

    const double K_rr();
    const double K_tt(double);
    const double K_pp();

    void HCS_rphi(const double &r, const double &phi, double& x, double& y, double& z) const;
    //double HCS_xy_z(const double &x, const double &y) const;

    double Ek;                                                  //kinetic energy
    double M_p;                                                 // momentum of particle
    double V_p;                                                 // velocity of praticle
    double r = (1.)*Unit::AU;                              // radial distance
    double theta = Unit::deg*(90)+pow(10.,-10.);
    double phi = pow(10.,-10.);
};


#endif