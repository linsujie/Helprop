#ifndef FCACHE_H
#define FCACHE_H

#include <functional>
#include <cmath>
#include <vector>

class fcache {
  public:
    fcache(const std::function<double(double)>& f, int npix_);
    ~fcache();

    double operator()(double x) const;

  private:
    const double pi_2 = std::acos(-1) * 2;
    int npix;
    double dphi;
    std::vector<double> vals;
};

#endif /* FCACHE_H */