#ifndef NEWTON_SOLVER_H
#define NEWTON_SOLVER_H
#include <functional>

class NewtonSolver {
  public:
    NewtonSolver() {}
    ~NewtonSolver() {}
    double solve(double x0, double err);
    void setF(const std::function<double(double)>& f_, double dx);
    void setFD(const std::function<double(double)>& fd_);

  private:
    std::function<double(double)> fd;
};

#endif /* NEWTON_SOLVER_H */