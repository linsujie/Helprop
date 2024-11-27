#include "newton_solver.h"
#include <iostream>
#include <cassert>

using namespace std;
double NewtonSolver::solve(double x0, double err) {
  double a = x0;
  double x = x0 - fd(x0);
  int i = 0;
  while (abs(x - a) > err) {
    a = x;
    i++;
    x = x - fd(x);
    if (i > 100) assert(false && "NewtonSolver::solve::too much iterations");
  }
  //cout << "Newton::iter" << i << endl;
  return a;
}

void NewtonSolver::setF(const function<double(double)> &f, double dx) {
  fd = [=](double x) { return f(x) / (f(x + dx) - f(x)) * dx; };
}

void NewtonSolver::setFD(const function<double(double)> &fd_) {
  fd = fd_;
}