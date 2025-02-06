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

using namespace std;

static char USAGE[] = R"(
This Routine is used to Check the simulation of one particle.

    Usage:
      ./ParticleTest [options] <outfile>

    Options:
      -h --help                         Show this help.
      -s SEED, --seed SEED              The global seed of this routine, it would be automatically given if not assigned.
      -m MASS, --mass MASS              The particle mass in GeV [default: 0.93827].
      -B B0, --B0 B0                    The magnetic strength around the Earth in nT [default: 5].
      -p POLARITY, --polarity POLARITY  The direction polarity of the magnetic field [default: -1].
      -a ANGLE, --angle ANGLE           Tilt angle of HCS in deg [default: 15].
      -D D, --D D                       Diffusion factor in unit 1e22 cm^2/s [default: 5].
      --indexA INDEXA                   Diffusion index a [default: 2].
      --ekin EKINS                     The ekin of particle [default: 0.5].
)";
int main(int argc, char* argv[]) {
  std::map<std::string, docopt::value> args = docopt::docopt(USAGE, {argv + 1, argv + argc}, true);
  auto fargs = [&](const string& k) { return atof(args.at(k).asString().c_str()); };

  HCS::angle = stod(args.at("--angle").asString()) * Unit::deg;
  HCS::hcsform = HCS::Kota_Jokipii;

  particle one(args);
  one.Ek = fargs("--ekin") * Unit::GeV;
  one.fix_seed = bool(args.at("--seed"));
  one.seed = one.fix_seed ? args.at("--seed").asLong() : 0;
  one.step(args.at("<outfile>").asString());

  return 0;
}