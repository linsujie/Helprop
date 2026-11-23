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
#include <algorithm>
#include <memory>
#include "loginterp.h"
#include "docopt.h"
#include "particle.h"
#include "IO.h"

using namespace std;
using namespace Unit;
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
    double ekmin = stod(eks[0]) * GeV;
    double ekmax = stod(eks[1]) * GeV;
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

// log-spaced grid, values in GeV (used for the refined matrix E_LIS axis)
vector<double> get_ekin_log(double emin, double emax, int n) {
  vector<double> ekin;
  for (int i = 0; i < n; i++)
    ekin.push_back(emin * pow(emax / emin, double(i) / (n - 1)));
  return ekin;
}

vector<particle> simulating(const particle& template_particle, int number, int th_num, const string& logname) {
  int n_per_thread = ceil(double(number) / th_num);

  vector<particle> Particle;
  Particle.resize(number);
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
      Particle[ip_local] = template_particle;
      if (Particle[ip_local].fix_seed)
        Particle[ip_local].seed += ip_local;

      if (ip_local % 1000 == 0)
        cerr << ">>particle " << ip_local << " seed " << Particle[ip_local].seed << ": "
          << " Ek " << template_particle.Ek / Unit::GeV
          << "GeV -> ";
      Particle[ip_local].step(logname);
      if (ip_local % 1000 == 0)
        cerr << Particle[ip_local].Ek / Unit::GeV << "GeV"
          << (Particle[ip_local].available ? " success" : " fail")
          << endl;
      assert(Particle[ip_local].available && "The particle is supposed to energetic enough to get out of the solar center.");
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

      if (i % 1000 == 0)
        cerr << ">>particle " << i << " seed " << Particle[i].seed << ": " << " Ek " << template_particle.Ek / Unit::GeV
          << "GeV -> ";
      Particle[i].step(logname);
      if (i % 1000 == 0)
        cerr << Particle[i].Ek / Unit::GeV << "GeV" << endl;
    }
  }
  return Particle;
}

const double m_proton = 0.938272 * GeV;
double ekin2p(double ekin, int A) { // ekin/nuc -> momentum
  return sqrt(ekin * (ekin + 2. * m_proton)) * A;
}
double p2ekin(double p, int A) { // momentum -> ekin/nuc
  return sqrt(p / A * p / A + m_proton * m_proton) - m_proton;
}
inline vector<double> get_bound(const vector<double>& x) {
  vector<double> bound;
  int n = x.size();
  bound.push_back(x[0]*sqrt(x[0]/x[1]));
  for (int i = 0; i < n-1; i++)
    bound.push_back(sqrt(x[i] * x[i + 1]));
  bound.push_back(x[n-1]*sqrt(x[n-1]/x[n-2]));

  return bound;
}

// counting the Green function matrix, the detail could be checked in the file modulation_matrix.pdf
vector<double> count_GreenFunction(const vector<particle>& Particle, const vector<double>& ekin) {
  assert(ekin.size() >= 2 && "At least two energy grids are required in the generation of Green Function matrix.");

  vector<double> bound = get_bound(ekin);

  vector<double> bin;
  bin.resize(ekin.size());

  int number = Particle.size();
  for (const auto& p : Particle) {
    int ibin = upper_bound(bound.begin(), bound.end(), p.Ek / p.A) - bound.begin();

    if (0 < ibin && ibin < bin.size() + 1)
      bin[ibin - 1] += 1;
  }

  double sum = 0;
  for (int i = 0; i < bin.size(); i++) sum += bin[i];
  for (auto& v : bin) v /= sum;

  return bin; // the returned matrix counting the probability ~ \int G(p, p') dp'
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
      --A-drift-HCS DHCS                The normalization of HCS drift [default: 1].
      --A-drift-gc DGC                  The normalization of gc drift [default: 1].
      -D D0, --D0 D0                    Reference diffusion coefficient in unit 1e22 cm^2/s [default: 5].
      -R R0, --R0 R0                    Reference rigidity for the diffusion coefficient in unit GV [default: 1].
      --indexA INDEXA                   Diffusion index a [default: 1].
      --indexB INDEXB                   Diffusion index b [default: 1].
      --etoa ETOA                       The ekin/nucleon of TOA spectrum assigned in format min,max,nbin in GeV, it would follow the input spec if not given.
      --elis ELIS                       The ekin/nucleon of LIS spectrum assigned in format min,max,nbin in GeV, it would follow the input spec or etoa if not given.
      --sample                          If given, to store the samples to the outmatrix or not, only available for BSON format.
      --read-meta                       Read the metadata (params, seed and sample arrays) stored in the input BSON file and pass it to the output [default: false].
      --spectrum-mode MODE              Produce the output spectrum from the Green function matrix or from per-particle integration: matrix or events [default: matrix].
      --matrix-refine NBIN              If >0, use a refined E_LIS grid with NBIN log-spaced bins covering the E_LIS tails for the matrix output and the matrix-mode spectrum [default: 0].
      --iotype IOTYPE                   The input/output type (TXT, CSV, or BSON) [default: TXT].
      --hcs-table HCS_TABLE             If the dir of the table generated by gen_HCS_distance_map given, to interpolate the HCS distance with it [default: ""].
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

  if (!io) {
    cerr << "Unknown --iotype '" << args.at("--iotype").asString()
         << "', expected TXT, CSV or BSON." << endl;
    return 1;
  }

  io->eunit = GeV;
  io->read_meta = args.at("--read-meta").asBool();
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
    if (!io->readspec(args.at("<inspec>").asString(), EIN, flux)) {
      cerr << "Failed to read the input spectrum from " << args.at("<inspec>").asString() << endl;
      return 1;
    }

  if (ETOA.empty()) ETOA = EIN;
  assert(!ETOA.empty() && "The ekin axis of TOA spectrum should be given.");

  if (ELIS.empty()) ELIS = EIN.empty() ? ETOA : EIN;

  // The output spectrum can be produced either from the Green function matrix
  // (fast, but limited by the E_LIS binning) or directly from the final energy
  // of every simulated particle (no binning bias, but no reusable matrix).
  const string spectrum_mode = args.at("--spectrum-mode").asString();
  if (spectrum_mode != "matrix" && spectrum_mode != "events") {
    cerr << "Unknown --spectrum-mode '" << spectrum_mode
         << "', expected matrix or events." << endl;
    return 1;
  }
  const bool per_event = spectrum_mode == "events";

  const bool write_matrix = bool(args.at("<outmatrix>"));
  const bool write_spectrum = bool(args.at("<outspec>"));
  const bool need_matrix = write_matrix || (write_spectrum && !per_event);
  const bool need_events = write_spectrum && per_event;

  // Optional refined E_LIS grid for the matrix output and matrix-mode spectrum.
  vector<double> EMAT = ELIS;
  const long matrix_refine = args.at("--matrix-refine").asLong();
  if (matrix_refine > 0) {
    if (matrix_refine < 2) {
      cerr << "--matrix-refine must be 0 or at least 2" << endl;
      return 1;
    }
    if (need_matrix) {
      // Cover the E_LIS tails: the kernel can extend by O(1 GeV) at low energy
      // and by O(1.5 GeV) at high energy, so do not tie the upper bound only to
      // ETOA.back() * 1.1 (that would truncate narrow ETOA ranges).
      const double refine_lo = max(ETOA.front() * 0.5, ETOA.front() - 2.0 * GeV);
      const double refine_hi = max(ETOA.back() * 1.1, ETOA.back() + 2.0 * GeV);
      EMAT = get_ekin_log(refine_lo, refine_hi, matrix_refine);
      cout << "refined matrix: " << EMAT.size() << " E_LIS bins in ["
           << EMAT.front() / GeV << ", " << EMAT.back() / GeV << "] GeV" << endl;
    }
  }

  unique_ptr<LogInterp> f_lis;
  if (!EIN.empty())
    f_lis.reset(new LogInterp(EIN, flux));
  if (write_spectrum && !f_lis) {
    cerr << "An input spectrum is required to produce an output spectrum" << endl;
    return 1;
  }

  // Samples describe the current run; never mix them with samples inherited
  // from the input file through --read-meta.
  if (args.at("--sample").asBool()) {
    io->seed.clear();
    io->etoa.clear();
    io->elis.clear();
    io->tcost.clear();
  }

  cout << "ETOA.size() = " << ETOA.size() << endl;
  vector<vector<double>> weight;  // Green function matrix
  vector<double> Ospec_event;     // per-particle spectrum

  int number = args.at("--number").asLong();
  int th_num = args.at("--nthread").asLong();
  bool fix_seed = bool(args.at("--seed"));
  long seed = fix_seed ? args.at("--seed").asLong() : 0;
  int A = args.at("--A").asLong();
  int Z = args.at("--Z").asLong();
  particle one(args);

  for (int i = 0; i < ETOA.size(); i++) {
    one.Ek = ETOA[i] * A;
    one.fix_seed = fix_seed;
    one.seed = seed + i * number;

    cout << "simulating Ek = " << one.Ek / Unit::GeV << endl;
    time_t start = clock();
    auto Particle = simulating(one, number, th_num, bool(args.at("--logname")) ? args.at("--logname").asString() : "");
    cout << "time costed per particle: " << (clock() - start) / (double)CLOCKS_PER_SEC * 1e3 / number << "ms" << endl;

    auto bin = vector<double>();
    if (need_matrix) {
      bin = count_GreenFunction(Particle, EMAT);
      weight.push_back(bin);
    }
    if (need_events) {
      const double ptoa = ekin2p(ETOA[i], A);
      double sum = 0.0;
      for (const auto& p : Particle) {
        const double elis = p.Ek / p.A;
        const double plis = ekin2p(elis, A);
        sum += (*f_lis)(elis) / plis / plis;
      }
      Ospec_event.push_back(sum / double(Particle.size()) * ptoa * ptoa);
    }
    if (args.at("--sample").asBool())
      for (const auto& p : Particle) {
        io->seed.push_back(p.seed);
        io->etoa.push_back(one.Ek / GeV);
        io->elis.push_back(p.Ek / GeV);
        io->tcost.push_back(p.Dt / day);
      }
  }

  if (write_matrix) {
    if (!io->writematrix(args.at("<outmatrix>").asString(), ETOA, EMAT, weight, write_mode)) {
      cerr << "Failed to write the output matrix to " << args.at("<outmatrix>").asString() << endl;
      return 1;
    }
    return 0;
  }

  vector<double> Ospec;
  if (per_event) {
    Ospec = Ospec_event;
  } else {
    vector<double> pLIS, pTOA;
    for (auto E : EMAT) pLIS.push_back(ekin2p(E, A));
    for (auto E : ETOA) pTOA.push_back(ekin2p(E, A));

    vector<double> FLIS;
    for (const auto& e : EMAT)
      FLIS.push_back((*f_lis)(e));

    for (int itoa = 0; itoa < weight.size(); itoa++) {
      double value = 0;

      for (int ilis = 0; ilis < EMAT.size(); ilis++)
        value += weight[itoa][ilis] * FLIS[ilis] / pLIS[ilis] / pLIS[ilis] * pTOA[itoa] * pTOA[itoa];

      Ospec.push_back(value);
    }
  }

  if (!io->writespec(args.at("<outspec>").asString(), ETOA, Ospec, write_mode)) {
    cerr << "Failed to write the output spectrum to " << args.at("<outspec>").asString() << endl;
    return 1;
  }

  return 0;
}
