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

void phi0_r_test(particle& p) {
  p.hcsform = particle::Kota_Jokipii;
  p.r = 14 * AU;
  p.theta = 70 * deg;
  double phi0 = p.Phi0_S(p.theta);
  double rlow, rup;
  p.r_bound(p.r, p.phi, phi0, rlow, rup);

  cout << "r rlow rup" << p.r / AU << " " << rlow / AU << " " << rup / AU << endl;

  cout << p.Theta_S(rlow, p.phi) / deg << endl;
  cout << "diff theta: " << setprecision(15) << (p.theta - p.Theta_S(rlow, p.phi)) / deg
    << " " << (p.theta - p.Theta_S(rup, p.phi)) / deg << endl;
}

int main() {
  particle p1;
  p1.hcsform = particle::Jokipii_Thomas;
  p1.mass = 0.5 * GeV;
  p1.B0 = 5 * nT;
  p1.polarity = 1;
  p1.angle = 45 * deg;
  p1.D = 5 * 1e22 * cm * cm / sec;
  p1.indexA = 2;
  p1.theta = pi / 2.0 + 1e-6;
  p1.Vs_eq = p1.Wind();
  p1.Bn = p1.B0 * AU * AU / 1.35883;

  p1.r = 8.0 * AU;
  p1.phi = 0.0;
  p1.theta = 30 * deg;
  p1.Theta_S(p1.r, p1.phi);

  vector<double> rs;
  vector<double> phis;
  vector<double> thetas;
  for (double i = 0; i < 50000; i += 1) {
    rs.push_back(double(rand()) / RAND_MAX * 20 * AU);
    phis.push_back(double(rand()) / RAND_MAX * 2 * pi);
    thetas.push_back(double(rand()) / RAND_MAX * pi);
  }

  vector<double> d1, d2;
  d1.reserve(phis.size());
  d2.reserve(phis.size());

  ////newton_test();
  ////phi0_r_test(p1);

  //clock_t t1 = clock();
  for (int i = 0; i < phis.size(); i++) {
    p1.r = rs[i];
    p1.phi = phis[i];
    p1.theta = thetas[i];
    d1.push_back(p1.get_HCS_distance());
  }
  //clock_t t2 = clock();
  for (int i = 0; i < phis.size(); i++) {
    p1.r = rs[i];
    p1.phi = phis[i];
    p1.theta = thetas[i];
    d2.push_back(p1.get_HCS_distance_old());
  }
  //clock_t t3 = clock();

  //cout << "time cost: " << (double)(t2 - t1) / CLOCKS_PER_SEC
  //  << " -> " << (double)(t3 - t2) / CLOCKS_PER_SEC
  //  << "  | " << (double)(t3 - t2) / (t2 - t1) << endl;
  
  double diff = 0;
  for (int i = 0; i < phis.size(); i++) {
    if (2 * (d1[i] - d2[i]) / (d1[i] + d2[i]) > diff)
      cout << "+++ " << rs[i] / AU << " " << thetas[i] / deg << " " << phis[i] / deg << " | " << d1[i] / AU << "  " << d2[i] / AU  << " " << 2 * (d1[i] - d2[i]) / (d1[i] + d2[i]) << endl;
    diff = fmax(diff, 2 * (d1[i] - d2[i]) / (d1[i] + d2[i]));
  }
  cout << "max difference: " << diff << endl;

  p1.r = rs[2197];// 1.983 * AU;
  p1.theta = thetas[2197];// 148.529 * deg;
  p1.phi = phis[2197];// 12.9353 * deg;
  double dold = p1.get_HCS_distance_old() / AU;
  double d = p1.get_HCS_distance() / AU;
  cout << dold << " " << d << endl;
  return 0;
}