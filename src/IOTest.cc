#include "IO.h"
#include "rfl.hpp"
#include "rfl/json.hpp"
#include "rfl/bson.hpp"
//#include "rfl/yaml.hpp"
#include <iostream>

using namespace std;
void write_test(IO& io, const string& specname, const string& matname) {
  vector<double> E = { 0.5, 1, 2, 4 };
  vector<double> F = { 1, 2, 3, 4 };
  vector<double> F2 = { 1, 3, 2, 6 };
  vector<vector<double> > M = { { 1, 2, 3, 4 }, { 5, 6, 7, 8 }, { 9, 10, 11, 12 }, { 13, 14, 15, 16 } };
  vector<vector<double> > M2 = { { 2, 2, 3, 4 }, { 7, 6, 7, 8 }, { 5, 10, 11, 12 }, { 13, 14, 15, 16 } };

  io.params["number"] = 1;
  io.params["B0"] = 15;
  io.writematrix(matname, E, E, M);
  io.writespec(specname, E, F, IO::RECREATE);

  io.params["number"] = 1000;
  io.params["B0"] = 1.5;

  io.writematrix(matname, E, E, M2, IO::APPEND);
  io.writespec(specname, E, F2, IO::APPEND);
}

void read_test(IO& io, const string& specname, const string& matname) {
  vector<double> Eread, Fread;
  vector<vector<double> > Mread;

  for (int ientry = 1; ientry < 3; ientry++) {
    io.readspec(specname, Eread, Fread, ientry);
    for (const auto& p : io.params)
      cout << p.first << ": " << p.second << endl;
    cout << "# E F" << endl;
    for (int i = 0; i < Eread.size(); i++)
      cout << Eread[i] << " " << Fread[i] << endl;

  
    io.readmatrix(matname, Eread, Eread, Mread, ientry);
    for (const auto& p : io.params)
      cout << p.first << ": " << p.second << endl;
    cout << "# M" << endl;
    for (int i = 0; i < Mread.size(); i++) {
      for (int j = 0; j < Mread[i].size(); j++)
        cout << Mread[i][j] << " ";
      cout << endl;
    }

  }
}

int main() {
  IO_BSON io_w, io_r;
  string specname = "spectest.bson";
  string  matname =  "mattest.bson";

  write_test(io_w, specname, matname);
  read_test(io_r, specname, matname);
  return 0;
}