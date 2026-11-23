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
#include <algorithm>
#include <exception>
#include <cstring>

using namespace std;

namespace {

// Removes the trailing CR of CRLF files.
string chomp(string line) {
  while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
    line.pop_back();
  return line;
}

// Parses whitespace separated doubles. Returns false when a token cannot be
// parsed or when the line does not contain any number.
bool parse_ws_doubles(const string& line, vector<double>& out) {
  out.clear();
  istringstream is(line);
  double val;
  while (is >> val) out.push_back(val);
  return is.eof() && !out.empty();
}

// Parses comma separated doubles. Empty fields (e.g. a trailing comma) are
// ignored, malformed fields make the parsing fail.
bool parse_csv_doubles(const string& line, vector<double>& out) {
  out.clear();
  string token;
  istringstream ss(line);
  while (getline(ss, token, ',')) {
    const auto first = token.find_first_not_of(" \t\r\n");
    if (first == string::npos) continue;
    const auto last = token.find_last_not_of(" \t\r\n");
    token = token.substr(first, last - first + 1);

    istringstream ts(token);
    double val;
    char extra;
    if (!(ts >> val) || (ts >> extra)) return false;
    out.push_back(val);
  }
  return !out.empty();
}

bool erase_prefix(string& line, const string& prefix) {
  if (line.rfind(prefix, 0) != 0) return false;
  line.erase(0, prefix.size());
  return true;
}

} // namespace

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

  for (auto& k : { "number", "A", "Z", "B0", "polarity", "angle", "D0", "R0", "indexA", "indexB" })
    params[k] = fargs(string("--") + k);

  params["matrix_refine"] = fargs("--matrix-refine");
  params["spectrum_mode"] = (args.at("--spectrum-mode").asString() == "events") ? 1.0 : 0.0;
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
  bool found = false;
  while (getline(spectrumFile, line)) {
    line = chomp(line);
    if (line == "# E F" && --ientry == 0) {
      found = true;
      break;
    }
  }
  if (!found) {
    cerr << "IO_TXT::readspec: could not find the requested entry in file " << filename << endl;
    return false;
  }

  while (getline(spectrumFile, line)) {
    line = chomp(line);
    if (line.empty()) continue;
    if (line[0] == '#') break;

    istringstream iss(line);
    double x, y;
    if (!(iss >> x >> y)) {
      cerr << "IO_TXT::readspec: malformed line in file " << filename << ": " << line << endl;
      return false;
    }
    E.push_back(x);
    F.push_back(y);
  }

  spectrumFile.close();

  if (E.empty()) {
    cerr << "IO_TXT::readspec: empty spectrum in file " << filename << endl;
    return false;
  }

  E = assign_unit(E, eunit);
  F = assign_unit(F, 1.0 / eunit);
  return true;
}

bool IO_TXT::writespec(const std::string& filename, const std::vector<double>& E_, const std::vector<double>& F_, WRITEMODE mode) const {
  ofstream of(filename, mode == APPEND ? ios::app : ios::trunc);
  if (!of) {
    cerr << "IO_TXT::writespec: cannot open file " << filename << endl;
    return false;
  }

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
  bool found = false;
  while (getline(data, line)) {
    line = chomp(line);
    if (!line.empty() && line[0] == '#' && --ientry == 0) {
      found = true;
      break;
    }
  }
  if (!found) {
    cerr << "IO_TXT::readmatrix: could not find the requested entry in file " << filename << endl;
    return false;
  }

  // The header line is "# <ELIS values>"; drop the leading '#' before parsing.
  string header = line;
  header.erase(0, header.find('#') + 1);
  if (!parse_ws_doubles(header, ELIS)) {
    cerr << "IO_TXT::readmatrix: malformed ELIS header in file " << filename << endl;
    return false;
  }

  while (getline(data, line)) {
    line = chomp(line);
    if (line.empty()) continue;
    if (line[0] == '#') break;

    vector<double> vals;
    if (!parse_ws_doubles(line, vals) || vals.size() != ELIS.size() + 1) {
      cerr << "IO_TXT::readmatrix: expected " << ELIS.size() + 1
           << " columns in file " << filename << ", got: " << line << endl;
      return false;
    }
    ETOA.push_back(vals.front());
    vals.erase(vals.begin());
    M.push_back(std::move(vals));
  }

  data.close();

  if (M.empty() || M.size() != ETOA.size()) {
    cerr << "IO_TXT::readmatrix: no data rows found in file " << filename << endl;
    return false;
  }

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
  if (ELIS_.empty()) {
    cerr << "IO_TXT::writematrix: empty ELIS axis" << endl;
    return false;
  }
  for (const auto& row : M_)
    if (row.size() != ELIS_.size()) {
      cerr << "IO_TXT::writematrix: matrix row has " << row.size()
           << " columns but ELIS has " << ELIS_.size() << endl;
      return false;
    }

  auto ETOA = split_unit(ETOA_, eunit);
  auto ELIS = split_unit(ELIS_, eunit);
  auto M = M_;// split_unit(M_, 1.0 / eunit / eunit / eunit);

  ofstream of(filename, mode == APPEND ? ios::app : ios::trunc);
  if (!of) {
    cerr << "IO_TXT::writematrix: cannot open file " << filename << endl;
    return false;
  }
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
  bool found = false;
  while (getline(data, line)) {
    line = chomp(line);
    if (line == "#E,F" && --ientry == 0) {
      found = true;
      break;
    }
  }
  if (!found) {
    cerr << "IO_CSV::readspec: could not find the requested entry in file " << filename << endl;
    return false;
  }

  while (getline(data, line)) {
    line = chomp(line);
    if (line.empty()) continue;
    if (line[0] == '#') break;

    vector<double> vals;
    if (!parse_csv_doubles(line, vals) || vals.size() < 2) {
      cerr << "IO_CSV::readspec: malformed line in file " << filename << ": " << line << endl;
      return false;
    }
    E.push_back(vals[0]);
    F.push_back(vals[1]);
  }

  if (E.empty()) {
    cerr << "IO_CSV::readspec: empty spectrum in file " << filename << endl;
    return false;
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
  if (!of) {
    cerr << "IO_CSV::writespec: cannot open file " << filename << endl;
    return false;
  }
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
  bool found = false;
  while (getline(data, line)) {
    line = chomp(line);
    if (line.rfind("#ELIS,", 0) == 0 && --ientry == 0) {
      found = true;
      break;
    }
  }
  if (!found) {
    cerr << "IO_CSV::readmatrix: could not find the requested entry in file " << filename << endl;
    return false;
  }

  string header = line;
  erase_prefix(header, "#ELIS,");
  if (!parse_csv_doubles(header, ELIS)) {
    cerr << "IO_CSV::readmatrix: malformed ELIS header in file " << filename << endl;
    return false;
  }

  if (!getline(data, line)) {
    cerr << "IO_CSV::readmatrix: missing ETOA header in file " << filename << endl;
    return false;
  }
  line = chomp(line);
  if (!erase_prefix(line, "#ETOA,") || !parse_csv_doubles(line, ETOA)) {
    cerr << "IO_CSV::readmatrix: malformed ETOA header in file " << filename << endl;
    return false;
  }

  while (getline(data, line)) {
    line = chomp(line);
    if (line.empty()) continue;
    // A new "#ELIS," line marks the beginning of the next entry.
    if (line.rfind("#ELIS,", 0) == 0) break;
    if (line[0] == '#') continue;

    vector<double> row;
    if (!parse_csv_doubles(line, row) || row.size() != ELIS.size()) {
      cerr << "IO_CSV::readmatrix: expected " << ELIS.size()
           << " columns in file " << filename << ", got: " << line << endl;
      return false;
    }
    M.push_back(std::move(row));
  }

  data.close();

  if (M.empty() || M.size() != ETOA.size()) {
    cerr << "IO_CSV::readmatrix: ETOA and M have different sizes in file " << filename << endl;
    return false;
  }

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
  if (ELIS_.empty()) {
    cerr << "IO_CSV::writematrix: empty ELIS axis" << endl;
    return false;
  }
  for (const auto& row : M_)
    if (row.size() != ELIS_.size()) {
      cerr << "IO_CSV::writematrix: matrix row has " << row.size()
           << " columns but ELIS has " << ELIS_.size() << endl;
      return false;
    }

  auto ETOA = split_unit(ETOA_, eunit);
  auto ELIS = split_unit(ELIS_, eunit);
  auto M = M_;// split_unit(M_, 1.0 / eunit / eunit / eunit);

  ofstream of(filename, mode == APPEND ? ios::app : ios::trunc);
  if (!of) {
    cerr << "IO_CSV::writematrix: cannot open file " << filename << endl;
    return false;
  }
  of << setprecision(8) << setiosflags(ios::scientific);
  of << "#ELIS,";
  for (int i = 0; i < ELIS.size(); i++) {
    if (i > 0) of << ",";
    of << ELIS[i];
  }
  of << endl;

  of << "#ETOA,";
  for (int i = 0; i < ETOA.size(); i++) {
    if (i > 0) of << ",";
    of << ETOA[i];
  }
  of << endl;

  of << "#Matrix" << endl;
  for (int irow = 0; irow < M.size(); irow++) {
    for (int icol = 0; icol < M[irow].size(); icol++) {
      if (icol > 0) of << ",";
      of << M[irow][icol];
    }
    of << endl;
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
    if (fread(&buf[0], 1, 4, datafile) != 4) {
      cerr << "IO_BSON::readspec/readmatrix: entry " << i + 1 << " not found in file " << filename << endl;
      buf.clear();
      break;
    }
    int entry_size = 0;
    memcpy(&entry_size, buf.data(), sizeof(entry_size));
    if (entry_size < 5) {
      cerr << "IO_BSON::readspec/readmatrix: invalid size of entry " << i + 1
           << " in file " << filename << endl;
      buf.clear();
      break;
    }
    buf.resize(entry_size);
    if (fread(&buf[4], 1, entry_size - 4, datafile) != (size_t)(entry_size - 4)) {
      cerr << "IO_BSON::readspec/readmatrix: entry " << i + 1 << " not found in file " << filename << endl;
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
  vector<double> etoa, elis, tcost;
};
bool IO_BSON::readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F, int ientry) {
  vector<char> buf = readbson(filename, ientry);

  if (buf.empty()) return false;

  SpecBson res;
  try {
    // DefaultIfMissing allows specs that only contain E and F (or any other
    // subset of the fields); the missing members keep their default values.
    res = rfl::bson::read<SpecBson, rfl::DefaultIfMissing>(buf).value();
  } catch (const std::exception& e) {
    cerr << "IO_BSON::readspec: cannot parse file " << filename << ": " << e.what() << endl;
    return false;
  }

  if (res.E.empty() || res.E.size() != res.F.size()) {
    cerr << "IO_BSON::readspec: E and F must be non-empty and of equal size in file "
         << filename << endl;
    return false;
  }

  E = assign_unit(res.E, eunit);
  F = assign_unit(res.F, 1.0 / eunit);

  // The metadata of the file is only adopted on demand: by default the
  // parameters of the current run are kept.
  if (read_meta) {
    for (const auto& [key, value] : res.params) params[key] = value;
    if (!res.seed.empty()) seed = res.seed;
    if (!res.etoa.empty()) etoa = res.etoa;
    if (!res.elis.empty()) elis = res.elis;
    if (!res.tcost.empty()) tcost = res.tcost;
  }
  return true;
}

bool IO_BSON::writespec(const std::string& filename, const std::vector<double>& E_, const std::vector<double>& F_, WRITEMODE mode) const {
  if (E_.size() != F_.size()) {
    cerr << "IO_BSON::writespec: E and F have different sizes" << endl;
    return false;
  }

  auto E = split_unit(E_, eunit);
  auto F = split_unit(F_, 1.0 / eunit);

  const auto spec = SpecBson{.params= params, .E = E, .F = F, .seed = seed, .etoa = etoa, .elis = elis, .tcost = tcost};
  vector<char> bspec = rfl::bson::write(spec);

  FILE *of = fopen(filename.c_str(), mode == APPEND ? "a" : "w");
  if (!of) {
    cerr << "IO_BSON::writespec: cannot open file " << filename << endl;
    return false;
  }
  const bool written = fwrite(&bspec[0], 1, bspec.size(), of) == bspec.size();
  fclose(of);
  if (!written) {
    cerr << "IO_BSON::writespec: failed to write file " << filename << endl;
    return false;
  }

  return true;
}

struct MatrixBson {
  map<string, double> params;
  std::vector<long> seed;
  std::vector<double> ETOA, ELIS, etoa, elis, tcost;
  std::vector<std::vector<double> > M;
};
bool IO_BSON::readmatrix(const std::string& filename, std::vector<double>& ETOA, std::vector<double>& ELIS, std::vector< std::vector<double> >& M, int ientry) {
  vector<char> buf = readbson(filename, ientry);
  if (buf.empty()) return false;

  MatrixBson res;
  try {
    res = rfl::bson::read<MatrixBson, rfl::DefaultIfMissing>(buf).value();
  } catch (const std::exception& e) {
    cerr << "IO_BSON::readmatrix: cannot parse file " << filename << ": " << e.what() << endl;
    return false;
  }

  if (res.ETOA.empty() || res.ELIS.empty() || res.M.size() != res.ETOA.size()) {
    cerr << "IO_BSON::readmatrix: inconsistent matrix dimensions in file " << filename << endl;
    return false;
  }
  for (const auto& row : res.M)
    if (row.size() != res.ELIS.size()) {
      cerr << "IO_BSON::readmatrix: matrix row has " << row.size()
           << " columns but ELIS has " << res.ELIS.size() << " in file " << filename << endl;
      return false;
    }

  ELIS = assign_unit(res.ELIS, eunit);
  ETOA = assign_unit(res.ETOA, eunit);
  M = res.M;

  if (read_meta) {
    for (const auto& [key, value] : res.params) params[key] = value;
    if (!res.seed.empty()) seed = res.seed;
    if (!res.etoa.empty()) etoa = res.etoa;
    if (!res.elis.empty()) elis = res.elis;
    if (!res.tcost.empty()) tcost = res.tcost;
  }
  //M = assign_unit(M, 1.0 / eunit / eunit / eunit);
  return true;
}

bool IO_BSON::writematrix(const std::string& filename, const std::vector<double>& ETOA_, const std::vector<double>& ELIS_, const std::vector< std::vector<double> >& M_, WRITEMODE mode) const {
  if (ETOA_.size() != M_.size()) {
    cerr << "IO_BSON::writematrix: ETOA and M have different sizes" << endl;
    return false;
  }
  if (ELIS_.empty()) {
    cerr << "IO_BSON::writematrix: empty ELIS axis" << endl;
    return false;
  }
  for (const auto& row : M_)
    if (row.size() != ELIS_.size()) {
      cerr << "IO_BSON::writematrix: matrix row has " << row.size()
           << " columns but ELIS has " << ELIS_.size() << endl;
      return false;
    }

  auto ETOA = split_unit(ETOA_, eunit);
  auto ELIS = split_unit(ELIS_, eunit);
  auto M = M_;// split_unit(M_, 1.0 / eunit / eunit / eunit);

  const auto matrix = MatrixBson{.params=params, .seed = seed, .ETOA = ETOA, .ELIS = ELIS, .etoa = etoa, .elis = elis, .tcost = tcost, .M = M};
  vector<char> bmatrix = rfl::bson::write(matrix);

  FILE *of = fopen(filename.c_str(), mode == APPEND ? "a" : "w");
  if (!of) {
    cerr << "IO_BSON::writematrix: cannot open file " << filename << endl;
    return false;
  }
  const bool written = fwrite(&bmatrix[0], 1, bmatrix.size(), of) == bmatrix.size();
  fclose(of);
  if (!written) {
    cerr << "IO_BSON::writematrix: failed to write file " << filename << endl;
    return false;
  }
  return true;
}
