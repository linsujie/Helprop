#ifndef FCACHE_H
#define FCACHE_H

#include <functional>
#include <cmath>
#include <vector>

const double pi_2 = std::acos(-1) * 2;

class fcache {
  public:
    fcache(const std::function<double(double)>& f, int npix_, double x_min = 0, double x_max = pi_2);
    ~fcache();

    double operator()(double x) const;

  public:
    // const double pi_2 = std::acos(-1) * 2;
    double x_min, x_max, x_range;
    bool triangularQ = false;
    int npix;
    double dx;
    std::vector<double> vals;
};

#endif /* FCACHE_H */