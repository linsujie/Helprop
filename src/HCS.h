#ifndef HCS_H
#define HCS_H
#include "Vec.hh"
#include "fcache.h"
#include "Unit.h"

enum Polygon { Dodecahedron, Icosahedron, test };

class HCS {
  public:
    enum HCSFORM { Jokipii_Thomas, Kota_Jokipii };
    static HCSFORM hcsform;
    static const double Omega;                        //angular velocity corresponding to 27.5 day
    static double angle;                        //tilt angle of HCS
    double t0 = 0.0;
    double t = 0.0;
    double Vs_eq;

    HCS(double Vs_eq_);
    ~HCS();

    double get_distance_old(double r, double theta, double phi, double ftol_abs) const;
    double get_distance(double r, double theta, double phi) const;
    double get_distance_polygon(double r, double theta, double phi, double Rg2, Polygon polygon = Polygon::Dodecahedron) const;

    double get_raw_distance(double r, double theta) const;

    double Theta_S(double, double) const;                                     //theat_s function
    static double Phi0_S(double);

    Vec norm_vec(const Vec& p_cs) const;

    inline double phi0(double r, double phi) const {
      return phi + r * Omega / Vs_eq - Omega * (t - t0);
    }

    void rphi(const double &r, const double &phi, double& x, double& y, double& z) const;
    //double HCS_xy_z(const double &x, const double &y) const;

private:

    void r_bound(double r, double phi, double phi0, double& rlow, double& rup) const;

    bool spiral_iterate(const Vec& target_point, Vec& p_cs, double& diter, double r, double theta, double phi) const;
    bool wave_iterate(const Vec& target_point, Vec& p_cs, double& diter, double r, double theta, double phi) const;
    bool point_iterate(const Vec& target_point, Vec& p_cs, Vec& dh, double& diter) const;

    static double Phi0_S_Jokipii_Thomas(double);
    static double Phi0_S_Kota_Jokipii(double);

    static double Theta_S_Jokipii_Thomas(double);
    static double Theta_S_Kota_Jokipii(double);

    inline double Theta_S_Jokipii_Thomas(double r, double phi) const {
      return Theta_S_Jokipii_Thomas(phi0(r, phi));
    }
    inline double Theta_S_Kota_Jokipii(double r, double phi) const {
      return Theta_S_Kota_Jokipii(phi0(r, phi));
    }
};

#endif /* HCS_H */
