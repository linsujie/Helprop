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

#include "KDInterpSide.h"
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
  hcs_real.resolution = 5e-4 * AU;
  cout << "generating..." << endl;
  clock_t t1 = clock();
  //HCS hcs_intp(p.Wind(), true);
  //KDInterpSide *kd = hcs_interp(hcs_real, true);
  KDInterpSide *kd = new KDInterpSide("dmap.bson");
  clock_t t2 = clock();
  cout << "generated..." << endl;
  cout << "time: " << (double(t2) - t1) / CLOCKS_PER_SEC << endl;
  cout << "storing..." << endl;
  kd->store_table("dmap.bson");
  cout << "stored..." << endl;

  //cout << "generating..." << endl;
  //KDInterp *kd = hcs_interp(hcs, true);
  //clock_t t2 = clock();
  //KDInterp kd(dist, {9.3035, 85.658, 92.46094}, {9.3036, 85.659, 92.46095}, 0.01, 1, "", tab_corr);

  //KDInterp *kd = hcs_intp.kd_tab;
  //double ang = kd->xmid[0], angw = kd->width[0];
  double r = kd->xmid[0], rw = kd->width[0];
  double theta = kd->xmid[1], thetaw = kd->width[1];
  double phi = kd->xmid[2], phiw = kd->width[2];

  cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
  const KDValueSide* pt;
  vector<double> vpoint;
  double maxval;
  ofstream hist("errhist.dat");

  int N = 100000;
  vector<double> angs(N),
    rs(N), thetas(N), phis(N), dreal(N), dint(N), err(N);
  srand(1);
  for (int i = 0; i < N; i++) {
    rs[i] = double(rand()) / RAND_MAX * 2 * rw + r - rw;
    thetas[i] = double(rand()) / RAND_MAX * 2 * thetaw + theta - thetaw;
    phis[i] = double(rand()) / RAND_MAX * 2 * phiw + phi - phiw;
  }

  clock_t t0 = clock();
  for (int i = 0; i < N; i++) {
    //HCS::angle = angs[i];
    dreal[i] = hcs_real.get_distance(rs[i], thetas[i], phis[i]) / AU;
  }
  clock_t treal = clock();
  for (int i = 0; i < N; i++)
    dint[i] = hcs_interp_eval(rs[i], thetas[i], phis[i], kd, hcs_real) / AU;
  clock_t tint = clock();

  //for (int i = 0; i < N; i++)
  //  cout << rs[i] / AU << " " << thetas[i] / deg << " " << phis[i] / deg << " " << dreal[i] << " " << dint[i] << endl;

  double err_min = -1, err_max = -1, err_med = 0;
  long Nerr = 0;
  for (int i = 0; i < N; i++) {
    if (fabs(dreal[i]) > fabs(dint[i])) continue;
    err[i] = fabs(dreal[i] - dint[i]);
    if (err_min == -1 || err[i] < err_min) err_min = err[i];
    if (err_max == -1 || err[i] > err_max) err_max = err[i];
    err_med += err[i] * err[i];
    //if (err[i] > 1) cout << i << endl;
    hist << err[i] << endl;
    Nerr++;
  }
  //int ix = 998;
  //double vint = hcs_interp_eval(rs[ix], thetas[ix], phis[ix], kd, hcs_real, true);
  //double vreal = hcs_real.get_distance(rs[ix], thetas[ix], phis[ix]);
  //cout << "vint: " << vint / AU << " vreal: " << vreal / AU << endl;

  err_med = sqrt(err_med / Nerr);
  cout << "err: [" << err_min << ", " << err_max << "] -> " << err_med << endl;
  cout << "real time: " << double(treal - t0) / CLOCKS_PER_SEC << "s " <<
    "intp time: " << double(tint - treal) / CLOCKS_PER_SEC << "s" << endl;

  cout << "Omega / Vs_eq = " << hcs_real.Omega / hcs_real.Vs_eq * AU << " (AU^-1)" << endl;

//  double par[3] = { 386513072760.0833, 3.13931196722899, 1.404561636981124 };
//  cout << hcs.get_distance(par[0], par[1], par[2]) / AU << endl;
  return 0;
}
