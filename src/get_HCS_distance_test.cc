#include <iostream>
#include <iomanip>
#include <unordered_map>
#include "newton_solver.h"
#include "root_finding.h"
#include "particle.h"
#include <cassert>

#include <functional>

#include "rfl.hpp"
#include "rfl/json.hpp"
#include "rfl/bson.hpp"

#include "ThreeDLookup.h"

#include "KDInterp.h"
#include "hcs_interp.h"

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

double func(const vector<double>& x) {
  double sum = 0;
  for (auto& v : x) sum += v*v;
  return sum;
}

int main() {
  particle p;
  HCS::angle = 15 * deg;
  HCS::hcsform = HCS::Kota_Jokipii;
  HCS hcs_real(p.Wind(), false);
  hcs_real.resolution = 1e-7 * AU;
  cout << "generating..." << endl;
  clock_t t1 = clock();
  HCS hcs_intp(p.Wind(), true);
  clock_t t2 = clock();
  cout << "generated..." << endl;
  cout << "time: " << (double(t2) - t1) / CLOCKS_PER_SEC << endl;

  //cout << "generating..." << endl;
  //KDInterp *kd = hcs_interp(hcs, true);
  //kd->store_table("dmap.bson");
  //KDInterp *kd = new KDInterp("dmap.bson");
  //clock_t t2 = clock();
  //KDInterp kd(dist, {9.3035, 85.658, 92.46094}, {9.3036, 85.659, 92.46095}, 0.01, 1, "", tab_corr);

  KDInterp *kd = hcs_intp.kd_tab;
  double r_low = kd->xmin[0], r_up = kd->xmax[0];
  double theta_low = kd->xmin[1], theta_up = kd->xmax[1];
  double phi_low = kd->xmin[2], phi_up = kd->xmax[2];

 long iter = 0;
  cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
  const KDValue* pt;
  vector<double> vpoint;
  double maxval;
  ofstream hist("errhist.dat");

  int N = 1000000;
  vector<double> rs(N), thetas(N), phis(N), dreal(N), dint(N), err(N);
  srand(1);
  for (int i = 0; i < N; i++) {
    rs[i] = double(rand()) / RAND_MAX * (r_up - r_low) + r_low;
    thetas[i] = double(rand()) / RAND_MAX * (theta_up - theta_low) + theta_low;
    phis[i] = double(rand()) / RAND_MAX * (phi_up - phi_low) + phi_low;
  }

  clock_t t0 = clock();
  for (int i = 0; i < N; i++)
    dreal[i] = hcs_real.get_distance(rs[i], thetas[i], phis[i]) / AU;
  clock_t treal = clock();
  for (int i = 0; i < N; i++)
    dint[i] = hcs_intp.get_distance(rs[i], thetas[i], phis[i]) / AU;
  clock_t tint = clock();

  double err_min = -1, err_max = -1, err_med = 0;
  for (int i = 0; i < N; i++) {
    err[i] = fabs(dreal[i] - dint[i]);
    if (err_min == -1 || err[i] < err_min) err_min = err[i];
    if (err_max == -1 || err[i] > err_max) err_max = err[i];
    err_med += err[i] * err[i];
  }
  err_med =  sqrt(err_med / N);
  cout << "err: [" << err_min << ", " << err_max << "] -> " << err_med / N << endl;
  cout << "real time: " << double(treal - t0) / CLOCKS_PER_SEC << "s " <<
    "intp time: " << double(tint - treal) / CLOCKS_PER_SEC << "s" << endl;

  cout << "Omega / Vs_eq = " << hcs_real.Omega / hcs_real.Vs_eq * AU << " (AU^-1)" << endl;

//  double par[3] = { 386513072760.0833, 3.13931196722899, 1.404561636981124 };
//  cout << hcs.get_distance(par[0], par[1], par[2]) / AU << endl;
  return 0;
}
