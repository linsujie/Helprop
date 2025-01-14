#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <cassert>
#include "docopt.h"
#include "particle.h"
#include "IO.h"

using namespace std;
mutex mtx;

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

vector<double> get_ekin(const string& ekin_opt) {
    vector<string> eks = split(ekin_opt, ",");
    double ekmin = stod(eks[0]);
    double ekmax = stod(eks[1]);
    int n = stoi(eks[2]);

    vector<double> ekin;
    for (double i = 0; i < n; i++)
        ekin.push_back(ekmin * pow(ekmax / ekmin, i / (n - 1)));

    return ekin;
}

vector<particle> simulating(const particle& template_particle, int number, int th_num, const string& logname) {
  int n_per_thread = ceil(double(number) / th_num);

  vector<particle> Particle;
  Particle.resize(number);
  vector<thread> threads;
  auto thread_run = [template_particle, &Particle, &logname](int iplow, int ipup) mutable {
    for (int i = iplow; i < ipup; i++) {
      Particle[i] = template_particle;
      if (Particle[i].fix_seed)
        Particle[i].seed += i;

      Particle[i].step(logname);
      cerr << ">>particle " << i << " seed " << Particle[i].seed << ": "
        << " Ek " << template_particle.Ek / Unit::GeV
        << "GeV -> " << Particle[i].Ek / Unit::GeV << "GeV" << endl;
    }
  };

  if (th_num > 1) {
    for (int ith = 0; ith < th_num; ith++)
        threads.emplace_back(thread_run, ith * n_per_thread, min((ith + 1) * n_per_thread, number));

    for (int j = 0; j < th_num; j++)
        threads[j].join();
  } else {
    for (int i = 0; i < number; i++) {
      Particle[i] = template_particle;
      if (Particle[i].fix_seed)
        Particle[i].seed += i;

      Particle[i].step(logname);
      cerr << ">>particle " << i << " seed " << Particle[i].seed << ": "
        << " Ek " << template_particle.Ek / Unit::GeV
        << "GeV -> " << Particle[i].Ek / Unit::GeV << "GeV" << endl;
    }
  }
  return Particle;
}

vector<double> count_GreenFunction(const vector<particle>& Particle, const vector<double>& ekin, int A = 1) {
  const double m_proton = 0.938272 * Unit::GeV;
  vector<double> momentum;
  for (const auto& e : ekin)
    momentum.push_back(sqrt(e * (e + 2. * A * m_proton)));

  assert(ekin.size() >= 2 && "At least two energy grids are required in the generation of Green Function matrix.");
  vector<double> ekin_bound;
  ekin_bound.push_back(ekin[0]*sqrt(ekin[0]/ekin[1]));
  for (int i = 0; i < ekin.size()-1; i++)
    ekin_bound.push_back(sqrt(ekin[i] * ekin[i + 1]));
  ekin_bound.push_back(ekin[ekin.size()-1]*sqrt(ekin[ekin.size()-1]/ekin[ekin.size()-2]));

  vector<double> bin;
  bin.resize(ekin.size());

  int number = Particle.size();
  for (const auto& p : Particle) {
    int ibin = upper_bound(ekin_bound.begin(), ekin_bound.end(), p.Ek / Unit::GeV) - ekin_bound.begin();

    if (0 < ibin && ibin < bin.size() + 1)
      bin[ibin - 1] += 1;
  }

  double sum = 0;
  for (int i = 0; i < bin.size(); i++) {
    bin[i] /= momentum[i] * momentum[i];
    sum += bin[i];
  }
  for (auto& v : bin) v /= sum;

  return bin;
}

static char USAGE[] = R"(
This Routine is used to simulate the modulation of particle within heliosphere.

    Usage:
      ./HelProp [options] <inspec> <outspec>
      ./HelProp [options] <outmatrix>

    Options:
      -h --help                         Show this help.
      -s SEED, --seed SEED              The global seed of this routine, it would be automatically given if not assigned.
      -n NTH, --nthread NTH             The number of threads used in this routine [default: 1].
      --number NUMBER                   The simulation particle number in each bin[default: 1000].
      -m MASS, --mass MASS              The particle mass in GeV [default: 0.93827].
      -B B0, --B0 B0                    The magnetic strength around the Earth in nT [default: 5].
      -p POLARITY, --polarity POLARITY  The direction polarity of the magnetic field [default: -1].
      -a ANGLE, --angle ANGLE           Tilt angle of HCS in deg [default: 15].
      -D D, --D D                       Diffusion factor in unit 1e22 cm^2/s [default: 5].
      --indexA INDEXA                   Diffusion index a [default: 2].
      --ekins EKINS                     The ekin assigned in format min,max,nbin in GeV, this option would only act when no inspec is assigned [default: 0.1,10,40].
      --sample                          If given, to store the samples to the outmatrix or not, only available for BSON format.
      --iotype IOTYPE                   The input/output type (TXT, CSV, or BSON) [default: TXT].
      --append                          Append the output to existing file [default: false].
      --logname LOGNAME                 The output logfile name.
)";
int main(int argc, char* argv[]) {
  std::map<std::string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);

  IO *io = NULL;
  if (args.at("--iotype").asString() == "TXT")
    io = new IO_TXT();
  else if (args.at("--iotype").asString() == "CSV")
    io = new IO_CSV();
  else if (args.at("--iotype").asString() == "BSON")
    io = new IO_BSON();

  io->set_params(args);

  vector<double> ekin;   // set spectrum energy bin
  vector<double> flux;   // boundary differential flux

  if (bool(args.at("<inspec>")))
    io->readspec(args.at("<inspec>").asString(), ekin, flux);

  if (ekin.empty())
    ekin = get_ekin(args.at("--ekins").asString());

  cout << "ekin.size() = " << ekin.size() << endl;
  vector<vector<double>> weight;  // possibility matrix

  int number = args.at("--number").asLong();
  int th_num = args.at("--nthread").asLong();
  particle one(args);
  bool fix_seed = bool(args.at("--seed"));
  long seed = fix_seed ? args.at("--seed").asLong() : 0;

  for (int i = 0; i < ekin.size(); i++) {
    one.Ek = ekin[i] * Unit::GeV;
    one.fix_seed = fix_seed;
    one.seed = seed + i * number;

    cout << "simulating Ek = " << one.Ek / Unit::GeV << endl;
    auto Particle = simulating(one, number, th_num, bool(args.at("--logname")) ? args.at("--logname").asString() : "");

    auto bin = count_GreenFunction(Particle, ekin);
    weight.push_back(bin);
    if (args.at("--sample").asBool())
      for (const auto& p : Particle) {
        io->seed.push_back(p.seed);
        io->ETOA.push_back(one.Ek / Unit::GeV);
        io->ELIS.push_back(p.Ek / Unit::GeV);
      }
  }

  if (bool(args.at("<outmatrix>"))) {
    io->writematrix(args.at("<outmatrix>").asString(), ekin, weight);
    return 0;
  }

  vector<double> Ospec;
  for (int i = 0; i < weight.size(); i++) {
    double value = 0;
    for (int j = 0; j < flux.size(); j++) {
      value += flux[j] * weight[i][j] / ekin[j] / ekin[j] * ekin[i] * ekin[i];
    }
    Ospec.push_back(value);
  }

  io->writespec(args.at("<outspec>").asString(), ekin, Ospec);

  return 0;
}
