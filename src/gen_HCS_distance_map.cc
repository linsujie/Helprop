#include <iostream>
#include <sstream>
#include <fstream>
#include <filesystem>

#include <sys/stat.h>
#include <sys/types.h>

#include "particle.h"
#include "Unit.h"
#include "docopt.h"
#include "hcs_interp.h"
#include "HCS.h"

using namespace std;
using namespace Unit;

void map_test(HCS& hcs_real, KDInterpSide* kd, int N = 10000) {
  vector<double> angs(N),
    rs(N), thetas(N), phis(N), dreal(N), dint(N), err(N);

  srand(1);

  double ang = kd->xmid[0], angw = kd->width[0];
  double r = kd->xmid[1], rw = kd->width[1];
  double theta = kd->xmid[2], thetaw = kd->width[2];
  double phi = kd->xmid[3], phiw = kd->width[3];
  for (int i = 0; i < N; i++) {
    angs[i] = double(rand()) / RAND_MAX * 2 * angw + ang - angw;
    rs[i] = double(rand()) / RAND_MAX * 2 * rw + r - rw;
    thetas[i] = double(rand()) / RAND_MAX * 2 * thetaw + theta - thetaw;
    phis[i] = double(rand()) / RAND_MAX * 2 * phiw + phi - phiw;

    thetas[i] = pi / 2 + thetas[i] * angs[i];
  }

  clock_t t0 = clock();
  for (int i = 0; i < N; i++) {
    HCS::angle = angs[i];
    dreal[i] = hcs_real.get_distance(rs[i], thetas[i], phis[i]) / AU;
  }
  clock_t treal = clock();
  for (int i = 0; i < N; i++)
    dint[i] = hcs_interp_eval(angs[i], rs[i], thetas[i], phis[i], kd, hcs_real) / AU;
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
    Nerr++;
  }

  if (ix != -1) {
    cout << ">>>>>>>>>>>>>>>>>>>>>>" << endl;
    double vint = hcs_interp_eval(angs[ix], rs[ix], thetas[ix], phis[ix], kd, hcs_real, true);
    HCS::angle = angs[ix];
    double vreal = hcs_real.get_distance(rs[ix], thetas[ix], phis[ix]);
    cout << "vint: " << vint / AU << " vreal: " << vreal / AU << endl;
    cout << "<<<<<<<<<<<<<<<<<<<<<<" << endl;
  }

  err_med = sqrt(err_med / Nerr);
  cout << "----------------------------------------------------" << endl;
  cout << "angle: " << ang / deg << " +- " << angw / deg << endl;
  cout << "err: [" << err_min << ", " << err_max << "] -> " << err_med << endl;
  cout << "real time: " << double(treal - t0) / CLOCKS_PER_SEC << "s " <<
    "intp time: " << double(tint - treal) / CLOCKS_PER_SEC << "s" << endl;
}

static char USAGE[] = R"(
This Routine is used to generate the table of HCS distance.

    Usage:
      ./gen_HCS_distance_map [options] <dir>

    Options:
      -h --help                         Show this help.
)";
int main(int argc, char *argv[]) {
  map<string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);

  string dir = args["<dir>"].asString();
  mkdir(dir.c_str(), 0755);
  ofstream namelist((dir + "/namelist.txt").c_str());

  particle p;
  HCS::angle = 15 * deg;
  HCS::hcsform = HCS::Kota_Jokipii;
  HCS hcs_real(p.Wind(), false);
  hcs_real.resolution = 5e-6 * AU;

  auto filename = [&](double ang_low, double ang_up) -> string {
    ostringstream fname;
    fname << dir << "/" << ang_low << "_" << ang_up << ".bson";
    return fname.str();
  };

  auto gen_kd = [&](double ang_low, double ang_up, int ix, bool pflag) -> bool {
    filesystem::path fname = filesystem::absolute(filename(ang_low, ang_up));
    namelist << fname.c_str() << endl;
    if (filesystem::exists(fname)) {
      cout << "-- file " << fname.c_str() << " already exists, skip --" << endl;
      return false;
    }

    cout << "generating..." << endl;
    clock_t t1 = clock();
    KDInterpSide *kd = hcs_interp(ang_low, ang_up, ix, pflag);
    clock_t t2 = clock();

    cout << "generated..." << endl;
    cout << "time: " << (double(t2) - t1) / CLOCKS_PER_SEC << endl;
    cout << "storing..." << endl;
    kd->store_table(fname);
    cout << "stored..." << endl;
    delete kd;
    return true;
  };

  double ALOW = 15, AUP = 75, DA = 5;
  vector<double> angs = {15, 20, 25, 30, 35, 40, 45, 50 };

  for (int i = 0; i < angs.size() - 1; i++)
    gen_kd(angs[i], angs[i + 1], 0, true);

  for (int i = 0; i < angs.size() - 1; i++) {
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    clock_t t1 = clock();
    KDInterpSide *kd = new KDInterpSide(filename(angs[i], angs[i + 1]));
    clock_t t2 = clock();
    cout << "reading cost: " << (double(t2) - t1) / CLOCKS_PER_SEC << endl;
    map_test(hcs_real, kd, 100000);
    delete kd;
  }
  
  return 0;
}