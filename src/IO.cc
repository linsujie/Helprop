#include "IO.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

bool IO::readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F) const {
  assert(false && "IO::readspec not implemented for selected type");
}
bool IO::writespec(const std::string& filename, const std::vector<double>& E, const std::vector<double>& F) const {
  assert(false && "IO::writespec not implemented for selected type");
}
bool IO::readmatrix(const std::string& filename, const std::vector<double>& E, std::vector< std::vector<double> >& M) const {
  assert(false && "IO::readmatrix not implemented for selected type");
}
bool IO::writematrix(const std::string& filename, const std::vector<double>& E, const std::vector< std::vector<double> >& M) const {
  assert(false && "IO::writematrix not implemented for selected type");
}

bool IO_TXT::readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F) const {
  E.clear();
  F.clear();

  ifstream spectrumFile(filename);
  if (!spectrumFile.is_open()) {
    cerr << "IO_TXT::readspec: could not open file " << filename << endl;
    return false;
  }

  string line;
  while (getline(spectrumFile, line)) {
    double x, y;
    istringstream iss(line);
    iss >> x >> y;
    E.push_back(x);
    F.push_back(y);
  }

  return true;
}

bool IO_TXT::writespec(const std::string& filename, const std::vector<double>& E, const std::vector<double>& F) const {
  ofstream of(filename);

  if (E.size() != F.size()) {
    cerr << "IO_TXT::writespec: E and F have different sizes" << endl;
    return false;
  }

  of << setprecision(8) << setiosflags(ios::scientific);
  for (int i = 0; i < E.size(); i++)
    of << E[i] << " " << F[i] << endl;

  return true;
}

bool IO_TXT::writematrix(const std::string& filename, const std::vector<double>& E, const std::vector< std::vector<double> >& M) const {
  ofstream of(filename);

  if (E.size() != M.size()) {
    cerr << "IO::writematrix: E and M have different sizes" << endl;
    return false;
  }

  of << setprecision(8) << setiosflags(ios::scientific);
  of << "# ";
  for (int i = 0; i < E.size(); i++)
    of << E[i] << " ";
  of << endl;

  for (int irow = 0; irow < M.size(); irow++) {
    for (int icol = 0; icol < M[irow].size(); icol++)
      of << M[irow][icol] << " ";
    of << endl;
  }

  return true;
}
