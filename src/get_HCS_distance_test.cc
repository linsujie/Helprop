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

void hcs_interp_test() {
  particle p;
  HCS::angle = 15 * deg;
  HCS::hcsform = HCS::Kota_Jokipii;
  HCS hcs_real(p.Wind(), "");
  hcs_real.resolution = 5e-4 * AU;

  auto gen_kd = [&](double ang_low, double ang_up, int ix, bool pflag) {
    cout << "generating..." << endl;
    clock_t t1 = clock();
    KDInterp *kd = hcs_interp(ang_low, ang_up, ix, pflag);
    clock_t t2 = clock();

    cout << "generated..." << endl;
    cout << "time: " << (double(t2) - t1) / CLOCKS_PER_SEC << endl;
    cout << "storing..." << endl;
    ostringstream fname;
    fname << "dmap" << ang_low << "_" << ang_up << ".bson";

    kd->store_table(fname.str());
    cout << "stored..." << endl;
    return kd;
  };

  auto read_kd = [&](const string& fname) {
    cout << "generating..." << endl;
    clock_t t1 = clock();
    KDInterp *kd = new KDInterp(fname);
    clock_t t2 = clock();

    cout << "generated..." << endl;
    cout << "time: " << (double(t2) - t1) / CLOCKS_PER_SEC << endl;
    return kd;
  };


  //KDInterp *kd1517 = read_kd("dmap15_20.bson");// gen_kd(15, 17, 0, true);
  KDInterp *kd1517 = gen_kd(15, 20, 0, true);
  //KDInterp *kd1516 = gen_kd(15, 16, 1, true);
  kd1517->show();
  //kd1516->show();
  //kd1617->show();

  //kd1517->kd->children[0]->show();
  //kd1517->kd->children[1]->show();

  //compare(*(kd1517->kd->children[0]), *(kd1516->kd));
  //compare(*(kd1517->kd->children[1]), *(kd1617->kd));

  //vec_t xmid_tot = {0.5, -0.25, -0.90625, 0.25};
  //vec_t xmid = {0, -0.25, -0.90625, 0.25};
  //print_block_d4(kd1517->kd->getkd(xmid_tot));
  //print_block_d4(kd1617->kd->getkd(xmid)->parent);

  KDInterp *kd = kd1517;
  double ang = kd->xmid[0], angw = kd->width[0];
  double r = kd->xmid[1], rw = kd->width[1];
  double theta = kd->xmid[2], thetaw = kd->width[2];
  double phi = kd->xmid[3], phiw = kd->width[3];

  cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
  const KDValueSide* pt;
  vector<double> vpoint;
  double maxval;
  ofstream hist("errhist.dat");

  int N = 10000;
  vector<double> angs(N),
    rs(N), thetas(N), phis(N), dreal(N), dint(N), err(N);
  srand(1);
  for (int i = 0; i < N; i++) {
    angs[i] = double(rand()) / RAND_MAX * 2 * angw + ang - angw;
    rs[i] = double(rand()) / RAND_MAX * 2 * rw + r - rw;
    thetas[i] = double(rand()) / RAND_MAX * 2 * thetaw + theta - thetaw;
    phis[i] = double(rand()) / RAND_MAX * 2 * phiw + phi - phiw;

    thetas[i] = pi / 2 + thetas[i] * angs[i];
  }

  clock_t t0 = clock();
  for (int i = 0; i < N; i++) {
    //cout << "----- " << angs[i] / deg  << " " << rs[i] / AU << " " << thetas[i] / deg << " " << phis[i] / deg << endl;

    HCS::angle = angs[i];
    dreal[i] = hcs_real.get_distance(rs[i], thetas[i], phis[i]) / AU;
  }
  clock_t treal = clock();
  for (int i = 0; i < N; i++)
    dint[i] = hcs_interp_eval(angs[i], rs[i], thetas[i], phis[i], kd, hcs_real) / AU;
  clock_t tint = clock();

  //for (int i = 0; i < N; i++)
  //  cout << rs[i] / AU << " " << thetas[i] / deg << " " << phis[i] / deg << " " << dreal[i] << " " << dint[i] << endl;

  double err_min = -1, err_max = -1, err_med = 0;
  long Nerr = 0;
  int ix = -1;
  for (int i = 0; i < N; i++) {
    if (fabs(dreal[i]) > fabs(dint[i])) continue;
    err[i] = fabs(dreal[i] - dint[i]);
    if (err_min == -1 || err[i] < err_min) err_min = err[i];
    if (err_max == -1 || err[i] > err_max) err_max = err[i];
    err_med += err[i] * err[i];
    if (fabs(err[i]) > 0.04) ix = i;
    hist << err[i] << endl;
    Nerr++;
  }

  if (ix != -1) {
    double vint = hcs_interp_eval(angs[ix], rs[ix], thetas[ix], phis[ix], kd, hcs_real, true);
    HCS::angle = angs[ix];
    double vreal = hcs_real.get_distance(rs[ix], thetas[ix], phis[ix]);
    cout << "vint: " << vint / AU << " vreal: " << vreal / AU << endl;
  }

  cout << "n_min_tol:";
  for (auto& v : kd->n_min_tol) cout << " " << v; cout << endl;
  err_med = sqrt(err_med / Nerr);
  cout << "err: [" << err_min << ", " << err_max << "] -> " << err_med << endl;
  cout << "angle: " << ang / deg << " +- " << angw / deg << endl;
  cout << "real time: " << double(treal - t0) / CLOCKS_PER_SEC << "s " <<
    "intp time: " << double(tint - treal) / CLOCKS_PER_SEC << "s" << endl;

  cout << "Omega / Vs_eq = " << hcs_real.Omega / hcs_real.Vs_eq * AU << " (AU^-1)" << endl;

}

void HCS_get_distance_interp_test() {
  particle p;
  HCS::angle = 15 * deg;
  HCS::hcsform = HCS::Kota_Jokipii;
  HCS hcs_real(p.Wind(), "");
  hcs_real.resolution = 5e-4 * AU;

  HCS hcs_intp(p.Wind(), "dmap");
  HCS::angle = 19 * deg;
  hcs_intp.get_distance(1 * AU, 0, 0);
  HCS::angle = 21 * deg;
  hcs_intp.get_distance(1 * AU, 0, 0);

  double ang = 20 * deg, angw = 4 * deg;
  double r = 60 * AU, rw = 55 * AU;
  double theta = 0, thetaw = 1.1;
  double phi = 180 * deg, phiw = 180 * deg;

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
    angs[i] = double(rand()) / RAND_MAX * 2 * angw + ang - angw;
    rs[i] = double(rand()) / RAND_MAX * 2 * rw + r - rw;
    thetas[i] = double(rand()) / RAND_MAX * 2 * thetaw + theta - thetaw;
    phis[i] = double(rand()) / RAND_MAX * 2 * phiw + phi - phiw;

    thetas[i] = pi / 2 + thetas[i] * angs[i];
  }

  clock_t t0 = clock();
  for (int i = 0; i < N; i++) {
    //cout << "----- " << angs[i] / deg  << " " << rs[i] / AU << " " << thetas[i] / deg << " " << phis[i] / deg << endl;

    HCS::angle = angs[i];
    dreal[i] = hcs_real.get_distance(rs[i], thetas[i], phis[i]) / AU;
  }
  clock_t treal = clock();
  for (int i = 0; i < N; i++) {
    HCS::angle = angs[i];
    dint[i] = hcs_intp.get_distance(rs[i], thetas[i], phis[i]) / AU;
  }
  clock_t tint = clock();

  double err_min = -1, err_max = -1, err_med = 0;
  long Nerr = 0;
  int ix = -1;
  for (int i = 0; i < N; i++) {
    if (fabs(dreal[i]) > fabs(dint[i])) continue;
    err[i] = fabs(dreal[i] - dint[i]);
    if (err_min == -1 || err[i] < err_min) err_min = err[i];
    if (err_max == -1 || err[i] > err_max) err_max = err[i];
    err_med += err[i] * err[i];
    if (fabs(err[i]) > 0.04) ix = i;
    hist << err[i] << endl;
    Nerr++;
  }

  err_med = sqrt(err_med / Nerr);
  cout << "err: [" << err_min << ", " << err_max << "] -> " << err_med << endl;
  cout << "angle: " << ang / deg << " +- " << angw / deg << endl;
  cout << "real time: " << double(treal - t0) / CLOCKS_PER_SEC << "s " <<
    "intp time: " << double(tint - treal) / CLOCKS_PER_SEC << "s" << endl;

  cout << "Omega / Vs_eq = " << hcs_real.Omega / hcs_real.Vs_eq * AU << " (AU^-1)" << endl;


}

int main() {
  HCS_get_distance_interp_test();
  return 0;
}
