#include <cmath>
#include <iostream>
#include "fcache.h"

using namespace std;

fcache::fcache(const std::function<double(double)>& f, int npix_): npix(npix_)
{
  dphi = pi_2 / (npix - 1);
  vals.reserve(npix);
  for (int i = 0; i < npix; i++)
    vals.push_back(f(dphi * i));
}

fcache::~fcache() {}

double fcache::operator()(double x) const {
  x -= floor(x / pi_2) * pi_2;

  int ilow = x / dphi;
  if (ilow == npix - 1) ilow = 0;
  int iup = ilow + 1;
  return vals[ilow] + (vals[iup] - vals[ilow]) * (x - ilow * dphi) / dphi;
}