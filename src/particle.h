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

    double D0;                                                   //diffusion factor
    double rigidity0;                                                  // reference rigidity of D0
    double indexA;                                              //power index related to rigidity of particle

    double Bn;
    double Bt = 0.0;                                            //
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
    double Heav();                                        //get heaviside function
    double B_r(double r, double theta, double heaviside) const;                  //radial magnetic field function
    double B_p(double r, double theta, double heaviside) const;                  //azimuthal magnetic field function

    double Kpara0() const;
    void K(double r, double theta, double heaviside, double kpara, double& krr, double& ktt, double& kpp, double& krp) const;
    void K(double r, double theta, double heaviside, double kpara, double& B, double& krr, double& ktt, double& kpp, double& krp) const;
    void coord_trans(double krr, double ktt, double kpp, double krp, double& dwr, double& dwt, double& dwp) const;
};

#endif