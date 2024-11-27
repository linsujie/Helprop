#include "IO.h"
#include "rfl.hpp"
#include "rfl/json.hpp"
//#include "rfl/bson.hpp"
#include "rfl/yaml.hpp"
#include <iostream>

using namespace std;
struct spec_json {
  vector<double> E;
  vector<double> F;
};
int main() {
  const auto jsspec = spec_json{.E = { 0.5, 1, 2, 4 }, .F = { 1, 2, 3, 4 } };

  const std::string jsstr = rfl::yaml::write(jsspec);

  return 1;
  IO_BSON io;
  string specname = "spectest.txt";
  string  matname =  "mattest.txt";

  vector<double> E = { 0.5, 1, 2, 4 };
  vector<double> F = { 1, 2, 3, 4 };
  vector<double> F2 = { 1, 3, 3, 4 };
  vector<vector<double> > M = { { 1, 2, 3, 4 }, { 5, 6, 7, 8 }, { 9, 10, 11, 12 }, { 13, 14, 15, 16 } };
  vector<vector<double> > M2 = { { 2, 2, 3, 4 }, { 7, 6, 7, 8 }, { 5, 10, 11, 12 }, { 13, 14, 15, 16 } };

  io.writematrix(matname, E, M);
  io.writematrix(matname, E, M2, IO::APPEND);

  io.writespec(specname, E, F, IO::RECREATE);
  io.writespec(specname, E, F2, IO::APPEND);

  vector<double> Eread, Fread;
  vector<vector<double> > Mread;

  for (int ientry = 1; ientry < 3; ientry++) {
    io.readspec(specname, Eread, Fread, ientry);
  
    cout << "# E F" << endl;
    for (int i = 0; i < E.size(); i++)
      cout << Eread[i] << " " << Fread[i] << endl;
  
//    io.readmatrix(matname, Eread, Mread, ientry);
//    cout << "# M" << endl;
//    for (int i = 0; i < Mread.size(); i++) {
//      for (int j = 0; j < Mread[i].size(); j++)
//        cout << Mread[i][j] << " ";
//      cout << endl;
//    }

  }

  return 0;
}