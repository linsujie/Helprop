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
  HCS hcs(p.Wind());
  HCS::hcsform = HCS::Kota_Jokipii;

  double x = -1.36394 * AU,
         y =  -4.95836 * AU,
         z =  1.33867 * AU;
  double r = binaryToDouble("0100001010101000101110011011000010100100100101000101101101011100");
  double theta = binaryToDouble("0011111111111101010011100001101101101001001011110100011000100111");
  double phi = binaryToDouble("0100000000000110000100100100101101010110110000010111110011101000");
  hcs.Theta_S(r, phi);
  hcs.resolution = 1e-7 * AU;

  clock_t t1 = clock();
  cout << "generating..." << endl;
  //KDInterp *kd = hcs_interp(hcs, true);
  //kd->store_table("dmap.bson");
  KDInterp *kd = new KDInterp("dmap.bson");
  clock_t t2 = clock();
  cout << "generated..." << endl;
  cout << "time: " << (double(t2) - t1) / CLOCKS_PER_SEC << endl;
  //KDInterp kd(dist, {9.3035, 85.658, 92.46094}, {9.3036, 85.659, 92.46095}, 0.01, 1, "", tab_corr);

  double r_low = kd->xmin[0], r_up = kd->xmax[0];
  double theta_low = kd->xmin[1], theta_up = kd->xmax[1];
  double phi_low = kd->xmin[2], phi_up = kd->xmax[2];

  double err_min = -1,
         err_med = -1,
         err_max = -1;
  long iter = 0;
  cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
  const KDValue* pt;
  vector<double> vpoint;
  double maxval;
  ofstream hist("errhist.dat");
  for (unsigned seed = 1;;seed++) {
    srand(seed);
    for (int i = 0; i < 100000; i++) {
      r = double(rand()) / RAND_MAX * (r_up - r_low) + r_low;
      theta = double(rand()) / RAND_MAX * (theta_up - theta_low) + theta_low;
      phi = double(rand()) / RAND_MAX * (phi_up - phi_low) + phi_low;

      double dreal = hcs.get_distance(r, theta, phi) / AU;
      double dint = hcs_interp_eval(r, theta, phi, kd, hcs) / AU;
      if (fabs(dint) > fabs(dreal)) continue;
      double err = fabs(dreal - dint);
      if (err_min == -1 || err < err_min) err_min = err;
      if (err_max == -1 || err > err_max) {
        vpoint = vector<double>({r / AU, theta / deg, phi / deg});
        maxval = dreal;
        err_max = err;
      }
      hist << err << endl;

      err_med += err;
      iter++;
      //cout << setprecision(10) << seed << " " << r / AU << " " << theta / deg << " " << phi / deg << " " << dreal << " " << dint << " " << err << endl;
    }
    break;
  }

  cout << "err: [" << err_min << ", " << err_max << "] -> " << err_med / iter << endl;

  cout << "Omega / Vs_eq = " << hcs.Omega / hcs.Vs_eq * AU << " (AU^-1)" << endl;

//  double par[3] = { 386513072760.0833, 3.13931196722899, 1.404561636981124 };
//  cout << hcs.get_distance(par[0], par[1], par[2]) / AU << endl;
  return 0;
}
