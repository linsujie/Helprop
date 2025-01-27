#ifndef LOGINTERP_H
#define LOGINTERP_H

#include <vector>

class LogInterp {
public:
    LogInterp(const std::vector<double> &x, const std::vector<double> &y);
    LogInterp(const double *x, const double *y, int n);
    void Init(const double *x, const double *y, int n);
    ~LogInterp();
    double operator()(double x) const;
private:
    std::vector<double> xlog, ylog;
};

#endif /* LOGINTERP_H */
