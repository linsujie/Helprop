#ifndef particle_h_
#define particle_h_

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include "HCS.h"
#include "docopt.h"
#include "Unit.h"
class particle {
    public:

    long seed = 0;
    long fix_seed = false;
    double A_drift = 1;

    double boundary = 100 * Unit::AU;                                     //boundary condition in AU
    double dt = 500. * Unit::sec;                                             //time interval per step
    double Vs;                                                  //solar wind velocity

    double polarity;                                            //field direction
    double B0;                                                  //magnetic strength in the Earth in T

    double k_xx;                                                //field along diffusion coefficient in local fram
    double k_yy;                                                //
    double k_zz;                                                //azimuth diffusion coefficient in local frame
    double k_rr;                                                //radial diffusion coefficient
    double k_tt;                                                //pole angle diffusion coefficient
    double k_pp;                                                //azimuthal diffusion coefficient
    double k_rr10 = 0.;
    double k_tt10 = 0.;
    double D;                                                   //diffusion factor
    double indexA;                                              //power index related to rigidity of particle

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
    static const double mp;                           //rest mass of proton

    double r, theta, phi;                              // particle position
    double r10 = 0;
    double Ek;                                                  //kinetic energy / nucleon
    double M_p;                                                 // momentum / nucleon
    double V_p;                                                 // velocity of praticle

    HCS hcs;

    particle(const std::map<std::string, docopt::value>& args);
    particle();
    ~particle();

    void step(const std::string& logname = "");                                   //simulate trajectory of particle

    public:
    double Wind() const;                                        //solar wind velocity function
    double Wind(double r, double theta, double phi, double angle) const;                                        //solar wind velocity function
    const double Heav();                                        //get heaviside function
    const double B_r(const double& r, const double& theta, const double& phi, const double& heaviside);                  //radial magnetic field function
    const double B_p(const double& r, const double& theta, const double& phi, const double& heaviside);                  //azimuthal magnetic field function

    const double K_rr(const double& r, const double& theta, const double& phi, const double& psi, const double& B, const double& B0, const double& M_p, const double& M_p0, const double& V_p);
    const double K_tt(const double& r, const double& theta, const double& phi, const double& psi, const double& B, const double& B0, const double& M_p, const double& M_p0, const double& V_p);
    const double K_pp(const double& r, const double& theta, const double& phi, const double& psi, const double& B, const double& B0, const double& M_p, const double& M_p0, const double& V_p);
};


#endif