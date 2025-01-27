#include <algorithm>
#include <cassert>
#include <cmath>

#include "loginterp.h"

using namespace std;

LogInterp::LogInterp(const double *x, const double *y, int n) {
  Init(x, y, n);
}

LogInterp::LogInterp(const std::vector<double>& x, const std::vector<double>& y) {
  assert(x.size() == y.size() && "LogInterp: x and y must have the same size");
  Init(&(x[0]), &(y[0]), x.size());
}

void LogInterp::Init(const double *x, const double *y, int n) {
  xlog.clear();
  ylog.clear();
  xlog.reserve(n);
  ylog.reserve(n);
  for (int i = 0; i < n; i++) {
    xlog.push_back(log(x[i]));
    ylog.push_back(log(y[i]));
  }
}

LogInterp::~LogInterp() {}

double LogInterp::operator()(double x) const {
  double logx = log(x);
  auto uiter = upper_bound(xlog.begin(), xlog.end(), logx);
  if (uiter == xlog.end()) uiter--;
  if (uiter == xlog.begin()) uiter++;

  int iup = uiter - xlog.begin();
  int ilow = iup - 1;

  return exp(ylog[ilow] + (logx - xlog[ilow]) * (ylog[iup] - ylog[ilow]) / (xlog[iup] - xlog[ilow]));
}