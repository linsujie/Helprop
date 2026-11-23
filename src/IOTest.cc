#include "IO.h"
#include "Unit.h"
#include "rfl.hpp"
#include "rfl/bson.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace Unit;

namespace {

int failures = 0;

void check(bool ok, const string& what) {
  if (!ok) {
    cerr << "FAIL: " << what << endl;
    ++failures;
  }
}

bool near(double a, double b) {
  const double scale = max({1.0, abs(a), abs(b)});
  return abs(a - b) <= 1e-6 * scale;
}

bool near(const vector<double>& a, const vector<double>& b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); i++)
    if (!near(a[i], b[i])) return false;
  return true;
}

bool near(const vector<vector<double> >& a, const vector<vector<double> >& b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); i++)
    if (!near(a[i], b[i])) return false;
  return true;
}

// A spec containing nothing but the two mandatory fields.
struct MinimalSpec {
  vector<double> E;
  vector<double> F;
};

void test_roundtrip(IO& io, const string& specname, const string& matname) {
  io.eunit = GeV;
  io.params["B0"] = 5;
  io.params["number"] = 1000;

  const vector<double> E = {0.5, 1.0, 2.0, 4.0};
  const vector<double> F1 = {1, 2, 3, 4};
  const vector<double> F2 = {4, 3, 2, 1};
  const vector<vector<double> > M1 = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
  const vector<vector<double> > M2 = {{16, 15, 14, 13}, {12, 11, 10, 9}, {8, 7, 6, 5}, {4, 3, 2, 1}};

  check(io.writespec(specname, E, F1), "writespec entry 1 " + specname);
  check(io.writespec(specname, E, F2, IO::APPEND), "writespec entry 2 " + specname);
  check(io.writematrix(matname, E, E, M1), "writematrix entry 1 " + matname);
  check(io.writematrix(matname, E, E, M2, IO::APPEND), "writematrix entry 2 " + matname);

  vector<double> e, f;
  check(io.readspec(specname, e, f, 1), "readspec entry 1 " + specname);
  check(near(e, E) && near(f, F1), "spectrum entry 1 " + specname);
  check(io.readspec(specname, e, f, 2), "readspec entry 2 " + specname);
  check(near(e, E) && near(f, F2), "spectrum entry 2 " + specname);

  vector<double> etoa, elis;
  vector<vector<double> > m;
  check(io.readmatrix(matname, etoa, elis, m, 1), "readmatrix entry 1 " + matname);
  check(near(etoa, E) && near(elis, E) && near(m, M1), "matrix entry 1 " + matname);
  check(io.readmatrix(matname, etoa, elis, m, 2), "readmatrix entry 2 " + matname);
  check(near(etoa, E) && near(elis, E) && near(m, M2), "matrix entry 2 " + matname);

  check(!io.readspec(specname, e, f, 3), "missing spectrum entry must fail " + specname);
  check(!io.readmatrix(matname, etoa, elis, m, 3), "missing matrix entry must fail " + matname);

  remove(specname.c_str());
  remove(matname.c_str());
}

void test_rectangular_matrix(IO& io, const string& name) {
  io.eunit = GeV;
  const vector<double> etoa = {0.5, 1.0, 2.0, 4.0};
  const vector<double> elis = {0.5, 1.0, 2.0};
  const vector<vector<double> > M = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}};

  check(io.writematrix(name, etoa, elis, M), "write rectangular matrix " + name);

  vector<double> etoa_read, elis_read;
  vector<vector<double> > M_read;
  check(io.readmatrix(name, etoa_read, elis_read, M_read, 1), "read rectangular matrix " + name);
  check(near(etoa_read, etoa) && near(elis_read, elis) && near(M_read, M),
        "rectangular matrix values " + name);

  remove(name.c_str());
}

void test_bson_meta() {
  const string name = "iotest_meta.bson";
  {
    IO_BSON w;
    w.eunit = 1.0;
    w.params["B0"] = 15;
    w.params["number"] = 1;
    w.seed = {11, 22};
    w.etoa = {0.5, 1.0};
    w.elis = {0.4, 0.8};
    w.tcost = {1.0, 2.0};
    check(w.writespec(name, {1, 2}, {3, 4}), "write spec with metadata");
  }

  vector<double> E, F;

  // By default the metadata of the file must not leak into the current run.
  IO_BSON r;
  r.eunit = 1.0;
  r.params["B0"] = 99;
  r.params["A"] = 4;
  check(r.readspec(name, E, F), "read spec without read_meta");
  check(r.params["B0"] == 99, "params must not be overwritten by default");
  check(r.params["A"] == 4, "params of the current run must be kept");
  check(r.seed.empty() && r.etoa.empty() && r.elis.empty() && r.tcost.empty(),
        "metadata must not be loaded by default");

  // With read_meta the file metadata is merged into the IO object.
  IO_BSON r2;
  r2.eunit = 1.0;
  r2.params["B0"] = 99;
  r2.params["A"] = 4;
  r2.read_meta = true;
  check(r2.readspec(name, E, F), "read spec with read_meta");
  check(r2.params["B0"] == 15, "params of the file must win with read_meta");
  check(r2.params["A"] == 4, "params absent from the file must be kept");
  check(r2.seed.size() == 2 && r2.etoa.size() == 2 && r2.elis.size() == 2 && r2.tcost.size() == 2,
        "sample arrays must be loaded with read_meta");

  // A BSON spec containing only E and F must be readable.
  const string minimal = "iotest_min.bson";
  const auto bytes = rfl::bson::write(MinimalSpec{.E = {1, 2, 3}, .F = {4, 5, 6}});
  FILE* fp = fopen(minimal.c_str(), "wb");
  check(fp != nullptr, "open minimal spec file");
  if (fp) {
    fwrite(bytes.data(), 1, bytes.size(), fp);
    fclose(fp);
  }

  IO_BSON r3;
  r3.eunit = 1.0;
  r3.params["B0"] = 7;
  check(r3.readspec(minimal, E, F), "read spec with only E and F");
  check(near(E, {1, 2, 3}) && near(F, {4, 5, 6}), "values of the minimal spec");
  check(r3.params["B0"] == 7, "minimal spec must not touch params");

  remove(name.c_str());
  remove(minimal.c_str());
}

} // namespace

int main() {
  {
    IO_TXT io;
    test_roundtrip(io, "iotest_spec.txt", "iotest_mat.txt");
    test_rectangular_matrix(io, "iotest_rect.txt");
  }
  {
    IO_CSV io;
    test_roundtrip(io, "iotest_spec.csv", "iotest_mat.csv");
    test_rectangular_matrix(io, "iotest_rect.csv");
  }
  {
    IO_BSON io;
    test_roundtrip(io, "iotest_spec.bson", "iotest_mat.bson");
    test_rectangular_matrix(io, "iotest_rect.bson");
  }
  test_bson_meta();

  if (failures == 0)
    cout << "All IO tests passed." << endl;
  else
    cout << failures << " IO test(s) failed." << endl;

  return failures == 0 ? 0 : 1;
}
