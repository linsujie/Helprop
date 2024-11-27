#include <cmath>
#include <iostream>
#include <iomanip>

#include "root_finding.h"

using namespace std;

double ridders_method(const function<double(const double)>& func, double x0, double x1, double ftol) {
  if (x0 == x1) return x0;
  double low = fmin(x0, x1),
         up = fmax(x0, x1),
         point = low,
         flow = func(low),
         fup = func(up),
         fpoint = flow;

  if (flow * fup > 0) {
    cerr << "Ridder's method: Could not find a root for " << x0 << " " << x1 << " in the same side" << endl;
    exit(1);
  }
  if (x1 == 0 && x0 == 0) return 0;

  int iter = 0;
  while (true) {
    double mid = (low + up) / 2,
      fmid = func(mid);
    if (low == low + (mid - low) / 1e1) break;

    point = mid + (mid - low) * (flow < 0 ? -1 : 1) * fmid / sqrt(fmid * fmid - flow * fup);
    fpoint = func(point);
    if (fabs(fpoint) < ftol || fpoint == flow || fpoint == fup) break;

    if (fpoint * flow < 0) {
      up = point;
      fup = fpoint;
    } else {
      low = point;
      flow = fpoint;
    }
  }

  return point;
}

double brents_method(const function<double(const double)>& func, double x0, double x1, double ftol, double xtol) {
  double f0 = func(x0),
         f1 = func(x1);
  if (f0 * f1 > 0) {
    cerr << "Brent's method: Could not find a root for " << x0 << " " << x1 << " in the same side" << endl;
    exit(1);
  }

  if (fabs(f0) < fabs(f1)) {
    double xtmp = x0; x0 = x1; x1 = xtmp;
    double ftmp = f0; f0 = f1; f1 = ftmp;
  }

  double xmid = x0;
  bool mflag = true;

  double s = 0, fs = 1, xval = 0;

  while (f0 != 0 && fs != 0 && fabs(x1 - x0) > ftol) {
    double fmid = func(xmid);
    if (f0 != fmid && f1 != fmid)
      s = x0 * f1 * fmid / ((f0 - f1) * (f0 - fmid)) + x1 * f0 * fmid / ((f1 - f0) * (f1 - fmid)) + xmid * f0 * f1 / ((fmid - f0) * (fmid - f1));
    else
      s = x1 - f1 * (x1 - x0) / (f1 - f0);

    if ((s < (3*x0 + x1) / 4 || x1 < s)
        || (mflag && fabs(s - x1) >= fabs(x1 - xmid) / 2)
        || (!mflag && fabs(s - x1) >= fabs(xmid - xval) / 2)
        || (mflag && fabs(x1 - xmid) < xtol)
        || (!mflag && fabs(xmid - xval) < xtol)) {
        s = (x0 + x1) / 2;
        mflag = true;
    } else mflag = false;

    fs = func(s);
    xval = xmid;
    xmid = x1;

    if (f0 * fs < 0) {
      x1 = s;
      f1 = fs;
    } else {
      x0 = s;
      f0 = fs;
    }

    if (fabs(f0) < fabs(f1)) {
      double xtmp = x0; x0 = x1; x1 = xtmp;
      double ftmp = f0; f0 = f1; f1 = ftmp;
    }
  }

  return s;
}