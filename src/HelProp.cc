#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <iomanip>
#include <cassert>
#include "loginterp.h"
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

vector<double> get_ekin(const map<string, docopt::value>& args, const string& key) {
  if (bool(args.at(key)))
    return get_ekin(args.at(key).asString());
  
  return vector<double>();
}

vector<particle> simulating(const particle& template_particle, int number, int th_num, const string& logname) {
  int n_per_thread = ceil(double(number) / th_num);

  vector<particle> Particle;
  Particle.reserve(number);
  vector<thread> threads;
  int ip = 0;
  auto thread_run = [template_particle, &Particle, &logname, &ip, number]() mutable {
    while (true) {
      int ip_local;
      {
        std::lock_guard<std::mutex> lock(mtx);
        if (ip >= number) break;
        ip_local = ip;
        ip++;
      }
      auto particle = template_particle;
      if (particle.fix_seed)
        particle.seed += ip_local;

      cerr << ">>particle " << ip_local << " seed " << particle.seed << ": "
        << " Ek " << template_particle.Ek / Unit::GeV
        << "GeV -> ";
      particle.step(logname);
      cerr << particle.Ek / Unit::GeV << "GeV" << endl;
      if (particle.Ek / Unit::GeV > 110) {
        std::cout << "db - " << __FILE__ << ":" << __LINE__ << " " << particle.seed << std::endl;
        assert(false);
      }
      {
        std::lock_guard<std::mutex> lock(mtx);
        Particle.push_back(particle);
      }
    }
  };

  if (th_num > 1) {
    for (int ith = 0; ith < th_num; ith++)
        threads.emplace_back(thread_run);

    for (int j = 0; j < th_num; j++)
        threads[j].join();

    if (Particle.size() != number) {
        std::cout << "Particle.size(): " << Particle.size() << std::endl;
        assert(false);
    }
  } else {
    for (int i = 0; i < number; i++) {
      Particle[i] = template_particle;
      if (Particle[i].fix_seed)
        Particle[i].seed += i;

      cerr << ">>particle " << i << " seed " << Particle[i].seed << ": "
        << " Ek " << template_particle.Ek / Unit::GeV
        << "GeV -> ";
      Particle[i].step(logname);
      cerr << Particle[i].Ek / Unit::GeV << "GeV" << endl;
    }
  }
  return Particle;
}

const double m_proton = 0.938272;
vector<double> count_GreenFunction(const vector<particle>& Particle, const vector<double>& ekin, int A = 1) {
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
      -A A, --A A                       The nucleon number A of particle [default: 1].
      -Z Z, --Z Z                       The charge number Z of particle [default: 1].
      -B B0, --B0 B0                    The magnetic strength around the Earth in nT [default: 5].
      -p POLARITY, --polarity POLARITY  The direction polarity of the magnetic field [default: -1].
      -a ANGLE, --angle ANGLE           Tilt angle of HCS in deg [default: 15].
      -D D0, --D0 D0                    Reference diffusion coefficient in unit 1e22 cm^2/s [default: 5].
      -R R0, --R0 R0                    Reference rigidity for the diffusion coefficient in unit GV [default: 1].
      --indexA INDEXA                   Diffusion index a [default: 1].
      --etoa ETOA                       The ekin/nucleon of TOA spectrum assigned in format min,max,nbin in GeV, it would follow the input spec if not given.
      --elis ELIS                       The ekin/nucleon of LIS spectrum assigned in format min,max,nbin in GeV, it would follow the input spec or etoa if not given.
      --sample                          If given, to store the samples to the outmatrix or not, only available for BSON format.
      --iotype IOTYPE                   The input/output type (TXT, CSV, or BSON) [default: TXT].
      --append                          Append the output to existing file [default: false].
      --logname LOGNAME                 The output logfile name.
)";
int main(int argc, char* argv[]) {
  std::map<std::string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);

  cout << "==================== HelProp ====================" << endl;

  IO *io = NULL;
  if (args.at("--iotype").asString() == "TXT")
    io = new IO_TXT();
  else if (args.at("--iotype").asString() == "CSV")
    io = new IO_CSV();
  else if (args.at("--iotype").asString() == "BSON")
    io = new IO_BSON();

  IO::WRITEMODE write_mode = args.at("--append").asBool() ? IO::APPEND : IO::RECREATE;

  io->set_params(args);

  HCS::angle = stod(args.at("--angle").asString()) * Unit::deg;
  HCS::hcsform = HCS::Kota_Jokipii;

  // set spectrum energy bin
  vector<double> EIN,
   ELIS = get_ekin(args, "--elis"),
   ETOA = get_ekin(args, "--etoa");
  vector<double> flux;   // boundary differential flux

  if (bool(args.at("<inspec>")))
    io->readspec(args.at("<inspec>").asString(), EIN, flux);

  if (ETOA.empty()) ETOA = EIN;
  assert(!ETOA.empty() && "The ekin axis of TOA spectrum should be given.");

  if (ELIS.empty()) ELIS = EIN.empty() ? ETOA : EIN;

  cout << "ETOA.size() = " << ETOA.size() << endl;
  vector<vector<double>> weight;  // Green function matrix

  int number = args.at("--number").asLong();
  int th_num = args.at("--nthread").asLong();
  particle one(args);
  bool fix_seed = bool(args.at("--seed"));
  long seed = fix_seed ? args.at("--seed").asLong() : 0;
  int A = 1;

  for (int i = 0; i < ETOA.size(); i++) {
    one.Ek = ETOA[i] * Unit::GeV;
    one.fix_seed = fix_seed;
    one.seed = seed + i * number;

    cout << "simulating Ek = " << one.Ek / Unit::GeV << endl;
    time_t start = clock();
    auto Particle = simulating(one, number, th_num, bool(args.at("--logname")) ? args.at("--logname").asString() : "");
    cout << "time costed per particle: " << (clock() - start) / (double)CLOCKS_PER_SEC * 1e3 / number << "ms" << endl;

    auto bin = count_GreenFunction(Particle, ELIS, A);
    weight.push_back(bin);
    if (args.at("--sample").asBool())
      for (const auto& p : Particle) {
        io->seed.push_back(p.seed);
        io->etoa.push_back(one.Ek / Unit::GeV);
        io->elis.push_back(p.Ek / Unit::GeV);
      }
  }

  if (bool(args.at("<outmatrix>"))) {
    io->writematrix(args.at("<outmatrix>").asString(), ETOA, ELIS, weight, write_mode);
    return 0;
  }

  LogInterp f_lis(EIN, flux);
  vector<double> FLIS;
  for (const auto& e : ELIS)
    FLIS.push_back(f_lis(e));

  vector<double> mLIS, mTOA;
  for (const auto& e : ETOA)
    mTOA.push_back(sqrt(e * (e + 2. * A * m_proton)));
  for (const auto& e : ELIS)
    mLIS.push_back(sqrt(e * (e + 2. * A * m_proton)));

  vector<double> Ospec;
  for (int itoa = 0; itoa < weight.size(); itoa++) {
    double value = 0;
    for (int ilis = 0; ilis < ELIS.size(); ilis++) {
      value += weight[itoa][ilis] * FLIS[ilis] / mLIS[ilis] / mLIS[ilis] * mTOA[itoa] * mTOA[itoa];
    }
    Ospec.push_back(value);
  }

  io->writespec(args.at("<outspec>").asString(), ETOA, Ospec, write_mode);

  return 0;
}
