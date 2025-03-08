#include <iostream>
#include <iomanip>
#include "newton_solver.h"
#include "root_finding.h"
#include "particle.h"

using namespace std;
using namespace Unit;

void brents_test() {
  int iter = 0;
  double val = brents_method([&](double x) {
    double f = (x*x*x) / 3 - x;
    cout << iter++ << " x: " << x << " " << f << endl; 
    return f;
  },
  0, 4.0, 1e-5, 1e-6);
  //double val = ridders_method([&](double x) -> double {
  //  double f = (x*x*x) / 3 - x;
  //  cout << iter++ << "  x: " << x << " " << f << endl; 
  //  return f;
  //},
  //0.1, 4.0, 1e-6);

  cout << "val: " << val << endl;
}

void newton_test() {
  NewtonSolver s1;
  s1.setFD([&](double x) {
	  double c = (x*x*x) / 3 - x;
	  double d = x*x - 1;
    cout << "x1: " << x << endl;
	  return (c / d);
  });

  NewtonSolver s2;
  s2.setF([&](double x) {
     cout << "x2: " << x << endl;
     return (x*x*x / 3 - x);
  }, 1e-10);

  cout << setprecision(10) << "s1: " << s1.solve(4, 1e-6) << endl;
  cout << setprecision(10) << "s2: " << s2.solve(4, 1e-7) << endl;
}

void phi0_r_test(HCS& hcs) {
  HCS::hcsform = HCS::Kota_Jokipii;
//  double r = 14 * AU;
//  double theta = 70 * deg;
//  double phi = 0;
//  double phi0 = hcs.Phi0_S(theta);
//  double rlow, rup;
//  hcs.r_bound(r, phi, phi0, rlow, rup);
//
//  cout << "r rlow rup" << r / AU << " " << rlow / AU << " " << rup / AU << endl;
//
//  cout << hcs.Theta_S(rlow, phi) / deg << endl;
//  cout << "diff theta: " << setprecision(15) << (theta - hcs.Theta_S(rlow, phi)) / deg
//    << " " << (theta - hcs.Theta_S(rup, phi)) / deg << endl;
}

int main() {
  particle p;
  HCS hcs(p.Wind());
  HCS::hcsform = HCS::Kota_Jokipii;
  hcs.angle = 15 * deg;
  hcs.Vs_eq = p.Wind();

  double x = -1.36394 * AU,
         y =  -4.95836 * AU,
         z =  1.33867 * AU;
  double r = 7138255256618;
  double theta = 1.824591131887;
  double phi = 5.075265712212;
//  double r = 8.8125444e+10;
//  double theta =   1.3001124;
//  double phi = 6.1798528;
  hcs.Theta_S(r, phi);

  cout << r / AU << " " << theta / deg << " " << phi / deg << endl;
  cout << hcs.get_distance(r, theta, phi) / AU << endl;

  return 0;
}
