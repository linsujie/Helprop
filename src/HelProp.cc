#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
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

vector<particle> simulating(const particle& template_particle, int number, int th_num) {
  int n_per_thread = ceil(double(number) / th_num);

  vector<particle> Particle;
  Particle.resize(number);
  vector<thread> threads;
  auto thread_run = [template_particle, &Particle](int iplow, int ipup) mutable {
    for (int i = iplow; i < ipup; i++) {
      Particle[i] = template_particle;
      if (Particle[i].fix_seed)
        Particle[i].seed += i;

      Particle[i].step();
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

      Particle[i].step();
      cerr << ">>particle " << i << " seed " << Particle[i].seed << ": "
        << " Ek " << template_particle.Ek / Unit::GeV
        << "GeV -> " << Particle[i].Ek / Unit::GeV << "GeV" << endl;
    }
  }
  return Particle;
}

vector<double> count_distribution(const vector<particle>& Particle, const vector<double>& ekin) {
  int number = Particle.size();

  vector<double> bin;

  bin.resize(ekin.size());
  for (int j = 0; j < number; j++) {
    double eng = Particle[j].Ek / Unit::GeV;
    for (int k = 0; k < ekin.size(); k++) {
      if (k == 0) {
        double x1 = log(ekin[k + 1]) / 2. - log(ekin[k]) / 2.;
        if (log(ekin[k]) - x1 <= log(eng) && log(eng) < log(ekin[k]) + x1)
          bin[k] += 1. / number;
      } 
      else if (0 < k && k < ekin.size() - 1) {
        double x0 = log(ekin[k - 1]) / 2. + log(ekin[k]) / 2.;
        double x1 = log(ekin[k + 1]) / 2. + log(ekin[k]) / 2.;
        if (x0 <= log(eng) && log(eng) < x1)
          bin[k] += 1. / number;
      } 
      else if (k == ekin.size() - 1) {
        double x1 = log(ekin[k]) / 2. - log(ekin[k - 1]) / 2.;
        if (log(ekin[k]) - x1 <= log(eng) && log(eng) < log(ekin[k]) + x1)
          bin[k] += 1. / number;
      }
    }
  }

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
      --iotype IOTYPE                   The input/output type (TXT, CSV, or BSON) [default: TXT].
      --append                          Append the output to existing file [default: false].
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
  int th_num = 1;//args.at("--nthread").asLong();
  particle one(args);
  bool fix_seed = bool(args.at("--seed"));
  long seed = fix_seed ? args.at("--seed").asLong() : 0;

  for (int i = 0; i < ekin.size(); i++) {
    one.Ek = ekin[i] * Unit::GeV;
    one.fix_seed = fix_seed;
    one.seed = seed + i * number;

    cout << "simulating Ek = " << one.Ek / Unit::GeV << endl;
    auto Particle = simulating(one, number, th_num);

    auto bin = count_distribution(Particle, ekin);
    weight.push_back(bin);
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