#include <cmath>
#include <iostream>
#include "fcache.h"

using namespace std;

fcache::fcache(const std::function<double(double)>& f, int npix_, double x_min_, double x_max_): npix(npix_), x_min(x_min_), x_max(x_max_)
{
  triangularQ = (x_min == 0 && x_max == pi_2);
  x_range = x_max - x_min;
  dx = x_range / (npix - 1);
  vals.reserve(npix);
  for (int i = 0; i < npix; i++)
    vals.push_back(f(x_min + dx * i));
}

fcache::~fcache() {}

double fcache::operator()(double x) const {
  if (triangularQ) {
    x -= floor(x / pi_2) * pi_2;
  }

  int ilow = x / dx;
  if (ilow == npix - 1) ilow = 0;
  int iup = ilow + 1;
  return vals[ilow] + (vals[iup] - vals[ilow]) * (x - ilow * dx) / dx;
}