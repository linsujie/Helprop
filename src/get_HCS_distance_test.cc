#include <iostream>
#include <iomanip>
#include "newton_solver.h"
#include "root_finding.h"
#include "particle.h"
#include <cassert>

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
  assert(1.234 == binaryToDouble(doubleToBinaryString(1.234)));
  double r = binaryToDouble("0100001010101000101110011011000010100100100101000101101101011100");
  double theta = binaryToDouble("0011111111111101010011100001101101101001001011110100011000100111");
  double phi = binaryToDouble("0100000000000110000100100100101101010110110000010111110011101000");
//  double r = 8.8125444e+10;
//  double theta =   1.3001124;
//  double phi = 6.1798528;
  hcs.Theta_S(r, phi);

  hcs.resolution = 1e-6 * AU;
  for (unsigned seed = 0;;seed++) {
    srand(seed);
    for (int i = 0; i < 10000; i++) {
      r = double(rand()) / RAND_MAX * 120 * AU;
      theta = double(rand()) / RAND_MAX * pi;
      phi = double(rand()) / RAND_MAX * 2 * pi;

      cout << setprecision(16) << seed << " " << r << " " << theta << " " << phi
           << flush;
      cout << " " << hcs.get_distance(r, theta, phi) / AU << " "
           << hcs.get_distance_old(r, theta, phi, 1e-4) / AU << endl;
    }
  }

  double par[3] = { 386513072760.0833, 3.13931196722899, 1.404561636981124 };
  cout << hcs.get_distance(par[0], par[1], par[2]) / AU << endl;

  return 0;
}
