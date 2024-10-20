#ifndef particle_h_
#define particle_h_

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include "docopt.h"
class particle {
    private:

    double AU = 1.496*pow(10.,8.);                              //the Sun-earth distance in km
    double light = 2.9979 * pow(10., 5.)/AU;                    //the speed of light in AU/s
    double pi = 3.1415926;
    double k_n = 1.*pow(10,12.)/(AU*AU);                        //normalization diffusion coefficient in AU^2/s corresponging to 10^22 cm^2/s
    double boundary = 110.;                                     //boundary condition in AU
    double dt = 2.;                                             //time interval per step
    double Vs;                                                  //solar wind velocity in unit AU/s
    double Omega = 2.644 * pow(10.,-6.);                        //angular velocity corresponding to 27.5 day

    double polarity;                                            //field direction
    double angle;                                               //tilt angle of HCS
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
    double Bt = 0. ;                                            //
    double psy;                                                 //pitch angle between field and radial direction
    double theta_s;                                             //tile angle of HCS at particle point
    double heaviside;                                           //field direction in particle point

    double rigidity;                                            //rigidity of particle related to kinetic and rest energy
    double A = 1.;                                              //nucleon number, proton by default
    double Z = 1.;                                              //charge number, proton by default
    double E0;                                                  //rest mass in MeV


    double Vdr_gc;                                              //drift to radial direction
    double Vdp_gc;                                              //drift to azimuthal direction
    double Vdt_gc;                                              //drift to pole direction
    double Vdr_HCS;                                             //drift at the HCS for thress direction
    double Vdp_HCS;
    double Vdt_HCS;



    public:
    particle(const std::map<std::string, docopt::value>& args);
    particle();
    ~particle();

    void step();                                                //simulate trajectory of particle

    const double Wind();                                        //solar wind velocity function
    const double Theta_S();                                     //theat_s function
    const double Heav();                                        //get heaviside function
    const double B_r(const double &heaviside);                  //radial magnetic field function
    const double B_p(const double &heaviside);                  //azimuthal magnetic field function

    const double K_rr();
    const double K_tt();
    const double K_pp();



    double Ek;                                                  //kinetic energy in MeV
    double M_p;                                                 // momentum of particle in MeV
    double V_p;                                                 // velocity of praticle in AU/s
    double r = 1. + pow(10.,-10.);                              // radial distance in AU
    double theta = pi*90./180.+pow(10.,-10.);
    double phy = pow(10.,-10.);



};


#endif