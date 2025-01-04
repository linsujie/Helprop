#include "IO.h"
#include "rfl.hpp"
#include "rfl/json.hpp"
#include "rfl/bson.hpp"

#include <bson.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

void IO::set_params(const std::map<std::string, docopt::value>& args) {
  auto fargs = [&](const std::string& key) -> double {
    return atof(args.at(key).asString().c_str());
  };

  for (auto& k : { "number", "mass", "B0", "polarity", "angle", "D", "indexA" })
    params[k] = fargs(string("--") + k);
}

bool IO::readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F, int ientry) {
  assert(false && "IO::readspec not implemented for selected type");
}
bool IO::writespec(const std::string& filename, const std::vector<double>& E, const std::vector<double>& F, WRITEMODE mode) const {
  assert(false && "IO::writespec not implemented for selected type");
}
bool IO::readmatrix(const std::string& filename, std::vector<double>& E, std::vector< std::vector<double> >& M, int ientry) {
  assert(false && "IO::readmatrix not implemented for selected type");
}
bool IO::writematrix(const std::string& filename, const std::vector<double>& E, const std::vector< std::vector<double> >& M, WRITEMODE mode) const {
  assert(false && "IO::writematrix not implemented for selected type");
}

bool IO_TXT::readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F, int ientry) {
  assert(ientry > 0 && "IO_TXT::readspec: ientry must be greater than 0");
  E.clear();
  F.clear();

  ifstream spectrumFile(filename);
  if (!spectrumFile.is_open()) {
    cerr << "IO_TXT::readspec: could not open file " << filename << endl;
    return false;
  }

  string line;
  int idata = 0;
  while (getline(spectrumFile, line)) {
    if (line == "# E F") idata++;
    if (idata == ientry) break;
  }

  while (getline(spectrumFile, line)) {
    if (line[0] == '#') break;

    double x, y;
    istringstream iss(line);
    iss >> x >> y;
    E.push_back(x);
    F.push_back(y);
  }

  spectrumFile.close();
  return true;
}

bool IO_TXT::writespec(const std::string& filename, const std::vector<double>& E, const std::vector<double>& F, WRITEMODE mode) const {
  ofstream of(filename, mode == APPEND ? ios::app : ios::trunc);

  if (E.size() != F.size()) {
    cerr << "IO_TXT::writespec: E and F have different sizes" << endl;
    return false;
  }

  of << "# E F" << endl;
  of << setprecision(8) << setiosflags(ios::scientific);
  for (int i = 0; i < E.size(); i++)
    of << E[i] << " " << F[i] << endl;

  of.close();
  return true;
}

bool IO_TXT::readmatrix(const std::string& filename, std::vector<double>& E, std::vector< std::vector<double> >& M, int ientry) {
  assert(ientry > 0 && "IO_TXT::readmatrix: ientry must be greater than 0");
  E.clear();
  M.clear();

  ifstream data(filename);
  if (!data.is_open()) {
    cerr << "IO_TXT::readmatrix: could not open file " << filename << endl;
    return false;
  }

  string line;
  double val;
  int idata = 0;
  while (getline(data, line)) {
    if (line[0] == '#') idata++;
    if (ientry == idata) break;
  }

  istringstream ishead(line);
  ishead >> val;
  while (ishead >> val) E.push_back(val);

  M.reserve(E.size());
  
  while (getline(data, line)) {
    if (line[0] == '#') break;
    M.resize(M.size() + 1);
    vector<double>& row = M.back();
    row.reserve(E.size());

    istringstream is(line);
    while (is >> val) row.push_back(val);
  }

  data.close();
  return true;
}

bool IO_TXT::writematrix(const std::string& filename, const std::vector<double>& E, const std::vector< std::vector<double> >& M, WRITEMODE mode) const {
  if (E.size() != M.size()) {
    cerr << "IO_TXT::writematrix: E and M have different sizes" << endl;
    return false;
  }

  ofstream of(filename, mode == APPEND ? ios::app : ios::trunc);
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

  of.close();
  return true;
}

vector<string> split(const string& str, const string& splitor)
{
  vector<string> result;

  int c_curr = 0,
      c_next = 0;
  while (c_next >= 0) {
    c_next = str.find(splitor, c_curr);
    result.push_back(str.substr(c_curr, c_next - c_curr));
    c_curr = c_next + 1;
  }

  return result;
}

template<typename T>
std::vector<T> split(const std::string& str, const std::string& splitor = " ") {
  auto strvec = split(str, splitor);
  std::vector<T> result;
  result.reserve(strvec.size());
  T tmp;

  for (auto str : strvec) {
    std::istringstream is(str);
    is >> tmp;
    result.push_back(tmp);
  }

  return result;
}

template <typename T>
std::string join(const std::vector<T>& vecs, const std::string& splitor = " ", int istart = 0, int iend = -1) {
  std::ostringstream os;
  unsigned size = vecs.size();
  istart = (int(istart + (std::abs(istart) / size + 1) * size)) % size;
  iend = (int(iend + (std::abs(iend) / size + 1) * size)) % size;

  if (istart > iend) return "";

  for (unsigned i = istart; i < iend; i++) os << vecs[i] << splitor;
  os << vecs[iend];

  return os.str();
}

bool IO_CSV::readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F, int ientry) {
  assert(ientry > 0 && "IO_CSV::readspec: ientry must be positive");
  E.clear();
  F.clear();

  ifstream data(filename);
  if (!data.is_open()) {
    cerr << "IO_CSV::readspec: cannot open file " << filename << endl;
    return false;
  }

  string line;
  int idata = 0;
  while (getline(data, line)) {
    if (line == "#E,F") idata++;
    if (idata == ientry) break;
  }

  while (getline(data, line)) {
    auto vals = split<double>(line, ",");
    E.push_back(vals[0]);
    F.push_back(vals[1]);
  }

  return true;
}


bool IO_CSV::writespec(const std::string& filename, const std::vector<double>& E, const std::vector<double>& F, WRITEMODE mode) const {
  if (E.size() != F.size()) {
    cerr << "IO_CSV::writespec: E and F have different sizes" << endl;
    return false;
  }

  ofstream of(filename, mode == APPEND ? ios::app : ios::trunc);
  of << "#E,F" << endl;
  of << setprecision(8) << setiosflags(ios::scientific);
  for (int i = 0; i < E.size(); i++)
    of << E[i] << "," << F[i] << endl;

  of.close();
  return true;
}

bool IO_CSV::readmatrix(const std::string& filename, std::vector<double>& E, std::vector< std::vector<double> >& M, int ientry) {
  assert(ientry > 0 && "IO_CSV::readmatrix: ientry must be positive");
  E.clear();
  M.clear();

  ifstream data(filename);
  if (!data.is_open()) {
    cerr << "IO_CSV::readmatrix: could not open file " << filename << endl;
    return false;
  }

  string line;
  int idata = 0;
  while (getline(data, line)) {
    if (line.substr(0, 3) == "#E,") idata++;
    if (idata == ientry) break;
  }

  double val;
  line.erase(0, 3);
  E = split<double>(line, ",");

  getline(data, line);

  M.reserve(E.size());
  while (getline(data, line)) {
    if (line[0] == '#') break;
    M.resize(M.size() + 1);
    M[M.size() - 1] = split<double>(line, ",");
  }

  data.close();
  return true;
}

bool IO_CSV::writematrix(const std::string& filename, const std::vector<double>& E, const std::vector< std::vector<double> >& M, WRITEMODE mode) const {
  if (E.size() != M.size()) {
    cerr << "IO_CSV::writematrix: E and M have different sizes" << endl;
    return false;
  }

  ofstream of(filename, mode == APPEND ? ios::app : ios::trunc);
  of << setprecision(8) << setiosflags(ios::scientific);
  of << "#E,";
  for (int i = 0; i < E.size() - 1; i++)
    of << E[i] << ",";
  of << E[E.size() - 1] << endl;

  of << "#Matrix" << endl;
  for (int irow = 0; irow < M.size(); irow++) {
    for (int icol = 0; icol < M[irow].size() - 1; icol++)
      of << M[irow][icol] << ",";
    of << M[irow][M[irow].size() - 1] << endl;
  }

  of.close();
  return true;
}

vector<char> readbson(const string& filename, int ientry) {
  vector<char> buf;

  FILE* datafile = fopen(filename.c_str(), "rb");
  if (!datafile) {
    cerr << "IO_BSON::readspec/readmatrix: cannot open file " << filename << endl;
    return buf;
  }

  for (int i = 0; i < ientry; i++) {
    buf.resize(4);
    int slen = fread(&buf[0], 1, 4, datafile);
    int entry_size = *((int*)(&buf[0]));
    buf.resize(entry_size);
    int sentry = fread(&buf[4], 1, entry_size - 4, datafile);
    if (slen != 4 || sentry != entry_size - 4) {
      cerr << "IO_BSON::readspec/readmatrix: no " << i << "th entry found in file " << filename << endl;
      buf.clear();
      break;
    }
  }
  fclose(datafile);

  return buf;
}
struct SpecBson {
  std::map<std::string, double> params;
  vector<double> E;
  vector<double> F;
};
bool IO_BSON::readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F, int ientry) {
  vector<char> buf = readbson(filename, ientry);

  if (buf.empty()) return false;

  const auto res = rfl::bson::read<SpecBson>(buf).value();
  E = res.E;
  F = res.F;
  params = res.params;

  return true;
}

bool IO_BSON::writespec(const std::string& filename, const std::vector<double>& E, const std::vector<double>& F, WRITEMODE mode) const {
  if (E.size() != F.size()) {
    cerr << "IO_BSON::writespec: E and F have different sizes" << endl;
    return false;
  }

  const auto spec = SpecBson{.params= params, .E = E, .F = F};
  vector<char> bspec = rfl::bson::write(spec);

  FILE *of = fopen(filename.c_str(), mode == APPEND ? "a" : "w");
  fwrite(&bspec[0], 1, bspec.size(), of);
  fclose(of);

  return true;
}

struct MatrixBson {
  map<string, double> params;
  std::vector<double> E;
  std::vector<std::vector<double> > M;
};
bool IO_BSON::readmatrix(const std::string& filename, std::vector<double>& E, std::vector< std::vector<double> >& M, int ientry) {
  vector<char> buf = readbson(filename, ientry);
  if (buf.empty()) return false;

  const auto res = rfl::bson::read<MatrixBson>(buf).value();
  E = res.E;
  M = res.M;
  params = res.params;

  return true;
}

bool IO_BSON::writematrix(const std::string& filename, const std::vector<double>& E, const std::vector< std::vector<double> >& M, WRITEMODE mode) const {
  if (E.size() != M.size()) {
    cerr << "IO_BSON::writematrix: E and M have different sizes" << endl;
    return false;
  }

  const auto matrix = MatrixBson{.params=params, .E = E, .M = M};
  vector<char> bmatrix = rfl::bson::write(matrix);

  FILE *of = fopen(filename.c_str(), mode == APPEND ? "a" : "w");
  fwrite(&bmatrix[0], 1, bmatrix.size(), of);
  fclose(of);
  return true;
}
