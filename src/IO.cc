#include "IO.h"
#include "Unit.h"
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

vector<double> IO::split_unit(const vector<double>& vec, const double unit) {
  vector<double> res = vec;
  for (auto& x : res) x /= unit;
  return res;
}
vector<double> IO::assign_unit(const vector<double>& vec, const double unit) {
  vector<double> res = vec;
  for (auto& x : res) x *= unit;
  return res;
}
vector<vector<double> > IO::split_unit(const vector<vector<double> >& vec, const double unit) {
  vector<vector<double> > res = vec;
  for (auto& r : res)
    for (auto& x : r) x /= unit;
  return res;
}
vector<vector<double> > IO::assign_unit(const vector<vector<double> >& vec, const double unit) {
  vector<vector<double> > res = vec;
  for (auto& r : res)
    for (auto& x : r) x *= unit;
  return res;
}


void IO::set_params(const std::map<std::string, docopt::value>& args) {
  auto fargs = [&](const std::string& key) -> double {
    return atof(args.at(key).asString().c_str());
  };

  for (auto& k : { "number", "A", "Z", "B0", "polarity", "angle", "D0", "R0", "indexA" })
    params[k] = fargs(string("--") + k);
}

bool IO::readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F, int ientry) {
  assert(false && "IO::readspec not implemented for selected type");
}
bool IO::writespec(const std::string& filename, const std::vector<double>& E, const std::vector<double>& F, WRITEMODE mode) const {
  assert(false && "IO::writespec not implemented for selected type");
}
bool IO::readmatrix(const std::string& filename, std::vector<double>& ETOA, std::vector<double>& ELIS, std::vector< std::vector<double> >& M, int ientry) {
  assert(false && "IO::readmatrix not implemented for selected type");
}
bool IO::writematrix(const std::string& filename, const std::vector<double>& ETOA, const std::vector<double>& ELIS, const std::vector< std::vector<double> >& M, WRITEMODE mode) const {
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
  E = assign_unit(E, eunit);
  F = assign_unit(F, 1.0 / eunit);
  return true;
}

bool IO_TXT::writespec(const std::string& filename, const std::vector<double>& E_, const std::vector<double>& F_, WRITEMODE mode) const {
  ofstream of(filename, mode == APPEND ? ios::app : ios::trunc);

  if (E_.size() != F_.size()) {
    cerr << "IO_TXT::writespec: E and F have different sizes" << endl;
    return false;
  }

  auto E = split_unit(E_, eunit);
  auto F = split_unit(F_, 1.0 / eunit);

  of << "# E F" << endl;
  of << setprecision(8) << setiosflags(ios::scientific);
  for (int i = 0; i < E.size(); i++)
    of << E[i] << " " << F[i] << endl;

  of.close();
  return true;
}

bool IO_TXT::readmatrix(const std::string& filename, std::vector<double>& ETOA, std::vector<double>& ELIS, std::vector< std::vector<double> >& M, int ientry) {
  assert(ientry > 0 && "IO_TXT::readmatrix: ientry must be greater than 0");
  ETOA.clear();
  ELIS.clear();
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
  while (ishead >> val) ELIS.push_back(val);

  M.reserve(ELIS.size());
  
  while (getline(data, line)) {
    if (line[0] == '#') break;
    M.resize(M.size() + 1);
    vector<double>& row = M.back();
    row.reserve(ELIS.size());

    istringstream is(line);
    is >> val;
    ETOA.push_back(val);
    while (is >> val) row.push_back(val);
  }

  data.close();
  ETOA = assign_unit(ETOA, eunit);
  ELIS = assign_unit(ELIS, eunit);
  //M = assign_unit(M, 1.0 / eunit / eunit / eunit);
  return true;
}

bool IO_TXT::writematrix(const std::string& filename, const std::vector<double>& ETOA_, const std::vector<double>& ELIS_, const std::vector< std::vector<double> >& M_, WRITEMODE mode) const {
  if (ETOA_.size() != M_.size()) {
    cerr << "IO_TXT::writematrix: ETOA and M have different sizes" << endl;
    return false;
  }

  auto ETOA = split_unit(ETOA_, eunit);
  auto ELIS = split_unit(ELIS_, eunit);
  auto M = M_;// split_unit(M_, 1.0 / eunit / eunit / eunit);

  ofstream of(filename, mode == APPEND ? ios::app : ios::trunc);
  of << setprecision(8) << setiosflags(ios::scientific);
  of << "# ";
  for (int i = 0; i < ELIS.size(); i++)
    of << ELIS[i] << " ";
  of << endl;

  for (int irow = 0; irow < M.size(); irow++) {
    of << ETOA[irow] << " ";
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

  E = assign_unit(E, eunit);
  F = assign_unit(F, 1.0 / eunit);
  return true;
}


bool IO_CSV::writespec(const std::string& filename, const std::vector<double>& E_, const std::vector<double>& F_, WRITEMODE mode) const {
  if (E_.size() != F_.size()) {
    cerr << "IO_CSV::writespec: E and F have different sizes" << endl;
    return false;
  }

  auto E = split_unit(E_, eunit);
  auto F = split_unit(F_, 1.0 / eunit);

  ofstream of(filename, mode == APPEND ? ios::app : ios::trunc);
  of << "#E,F" << endl;
  of << setprecision(8) << setiosflags(ios::scientific);
  for (int i = 0; i < E.size(); i++)
    of << E[i] << "," << F[i] << endl;

  of.close();
  return true;
}

bool IO_CSV::readmatrix(const std::string& filename, std::vector<double>& ETOA, std::vector<double>& ELIS, std::vector< std::vector<double> >& M, int ientry) {
  assert(ientry > 0 && "IO_CSV::readmatrix: ientry must be positive");
  ETOA.clear();
  ELIS.clear();
  M.clear();

  ifstream data(filename);
  if (!data.is_open()) {
    cerr << "IO_CSV::readmatrix: could not open file " << filename << endl;
    return false;
  }

  string line;
  int idata = 0;
  while (getline(data, line)) {
    if (line.substr(0, 6) == "#ELIS,") idata++;
    if (idata == ientry) break;
  }

  line.erase(0, 3);
  ELIS = split<double>(line, ",");

  getline(data, line);
  line.erase(0, 3);
  ETOA = split<double>(line, ",");

  double val;
  getline(data, line);

  M.reserve(ETOA.size());
  while (getline(data, line)) {
    if (line[0] == '#') break;
    M.resize(M.size() + 1);
    M[M.size() - 1] = split<double>(line, ",");
  }

  data.close();

  ETOA = assign_unit(ETOA, eunit);
  ELIS = assign_unit(ELIS, eunit);
  //M = assign_unit(M, 1.0 / eunit / eunit / eunit);
  return true;
}

bool IO_CSV::writematrix(const std::string& filename, const std::vector<double>& ETOA_, const std::vector<double>& ELIS_, const std::vector< std::vector<double> >& M_, WRITEMODE mode) const {
  if (ETOA_.size() != M_.size()) {
    cerr << "IO_CSV::writematrix: ETOA and M have different sizes" << endl;
    return false;
  }

  auto ETOA = split_unit(ETOA_, eunit);
  auto ELIS = split_unit(ELIS_, eunit);
  auto M = M_;// split_unit(M_, 1.0 / eunit / eunit / eunit);

  ofstream of(filename, mode == APPEND ? ios::app : ios::trunc);
  of << setprecision(8) << setiosflags(ios::scientific);
  of << "#ELIS,";
  for (int i = 0; i < ELIS.size() - 1; i++)
    of << ELIS[i] << ",";
  of << ELIS[ELIS.size() - 1] << endl;

  of << "#ETOA,";
  for (int i = 0; i < ELIS.size() - 1; i++)
    of << ETOA[i] << ",";
  of << ETOA[ETOA.size() - 1] << endl;

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
  vector<long> seed;
  vector<double> etoa, elis;
};
bool IO_BSON::readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F, int ientry) {
  vector<char> buf = readbson(filename, ientry);

  if (buf.empty()) return false;

  const auto res = rfl::bson::read<SpecBson>(buf).value();
  E = res.E;
  F = res.F;
  etoa = res.etoa;
  elis = res.elis;
  params = res.params;

  E = assign_unit(E, eunit);
  F = assign_unit(F, 1.0 / eunit);
  return true;
}

bool IO_BSON::writespec(const std::string& filename, const std::vector<double>& E_, const std::vector<double>& F_, WRITEMODE mode) const {
  if (E_.size() != F_.size()) {
    cerr << "IO_BSON::writespec: E and F have different sizes" << endl;
    return false;
  }

  auto E = split_unit(E_, eunit);
  auto F = split_unit(F_, 1.0 / eunit);

  const auto spec = SpecBson{.params= params, .E = E, .F = F, .seed = seed, .etoa = etoa, .elis = elis};
  vector<char> bspec = rfl::bson::write(spec);

  FILE *of = fopen(filename.c_str(), mode == APPEND ? "a" : "w");
  fwrite(&bspec[0], 1, bspec.size(), of);
  fclose(of);

  return true;
}

struct MatrixBson {
  map<string, double> params;
  std::vector<long> seed;
  std::vector<double> ETOA, ELIS, etoa, elis;
  std::vector<std::vector<double> > M;
};
bool IO_BSON::readmatrix(const std::string& filename, std::vector<double>& ETOA, std::vector<double>& ELIS, std::vector< std::vector<double> >& M, int ientry) {
  vector<char> buf = readbson(filename, ientry);
  if (buf.empty()) return false;

  const auto res = rfl::bson::read<MatrixBson>(buf).value();
  ELIS = res.ELIS;
  ETOA = res.ETOA;
  etoa = res.etoa;
  elis = res.elis;
  M = res.M;
  params = res.params;

  ETOA = assign_unit(ETOA, eunit);
  ELIS = assign_unit(ELIS, eunit);
  //M = assign_unit(M, 1.0 / eunit / eunit / eunit);
  return true;
}

bool IO_BSON::writematrix(const std::string& filename, const std::vector<double>& ETOA_, const std::vector<double>& ELIS_, const std::vector< std::vector<double> >& M_, WRITEMODE mode) const {
  if (ETOA_.size() != M_.size()) {
    cerr << "IO_BSON::writematrix: ETOA and M have different sizes" << endl;
    return false;
  }

  auto ETOA = split_unit(ETOA_, eunit);
  auto ELIS = split_unit(ELIS_, eunit);
  auto M = M_;// split_unit(M_, 1.0 / eunit / eunit / eunit);

  const auto matrix = MatrixBson{.params=params, .seed = seed, .ETOA = ETOA, .ELIS = ELIS, .etoa = etoa, .elis = elis, .M = M};
  vector<char> bmatrix = rfl::bson::write(matrix);

  FILE *of = fopen(filename.c_str(), mode == APPEND ? "a" : "w");
  fwrite(&bmatrix[0], 1, bmatrix.size(), of);
  fclose(of);
  return true;
}
