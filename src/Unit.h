#ifndef UNIT_H
#define UNIT_H

#include <cmath>

namespace Unit {
  const double m = 1;
  const double cm = 1e-2 * m;
  const double km = 1e3 * m;
  const double AU = 1.496*std::pow(10.,11.) * m;

  const double sec = 2.99792e8;
  const double min = 60 * sec;
  const double hr = 3600 * sec;
  const double day = 86400 * sec;

  const double c_speed = 2.99792e8 * m / sec;

  const double kg = 1;
  const double g = 1e-3 * kg;
  const double ton = 1e3 * kg;

  const double J = 1 * kg * m * m / sec / sec;

  const double pi = std::acos(-1);
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

  const double GV = GeV / e;
  const double MV = MeV / e;
};

#endif /* UNIT_H */
