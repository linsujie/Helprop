#ifndef IO_H
#define IO_H

#include <map>
#include <vector>
#include <string>
#include "docopt.h"

class IO {
public:
  IO() {}
  ~IO() {}
  enum WRITEMODE { RECREATE, APPEND };
  std::map<std::string, double> params;
  std::vector<long> seed;
  std::vector<double> etoa, elis;
  double eunit;
  void set_params(const std::map<std::string, docopt::value>& args);
  virtual bool readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F, int ientry = 1);
  virtual bool writespec(const std::string& filename, const std::vector<double>& E, const std::vector<double>& F, WRITEMODE mode = RECREATE) const;

  virtual bool readmatrix(const std::string& filename, std::vector<double>& ETOA, std::vector<double>& ELIS, std::vector< std::vector<double> >& M, int ientry = 1);
  virtual bool writematrix(const std::string& filename, const std::vector<double>& ETOA, const std::vector<double>& ELIS, const std::vector< std::vector<double> >& M, WRITEMODE mode = RECREATE) const;

  static std::vector<double> split_unit(const std::vector<double>& vec, const double unit);
  static std::vector<double> assign_unit(const std::vector<double>& vec, const double unit);
  static std::vector<std::vector<double> > split_unit(const std::vector<std::vector<double> >& vec, const double unit);
  static std::vector<std::vector<double> > assign_unit(const std::vector<std::vector<double> >& vec, const double unit);
};

class IO_TXT : public IO {
public:
  IO_TXT() {}
  ~IO_TXT() {}

  bool readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F, int ientry = 1);
  bool writespec(const std::string& filename, const std::vector<double>& E, const std::vector<double>& F, WRITEMODE mode = RECREATE) const;

  bool readmatrix(const std::string& filename, std::vector<double>& ETOA, std::vector<double>& ELIS, std::vector< std::vector<double> >& M, int ientry = 1);
  bool writematrix(const std::string& filename, const std::vector<double>& ETOA, const std::vector<double>& ELIS, const std::vector< std::vector<double> >& M, WRITEMODE mode = RECREATE) const;
};

class IO_CSV : public IO {
public:
  IO_CSV() {}
  ~IO_CSV() {}

  bool readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F, int ientry = 1);
  bool writespec(const std::string& filename, const std::vector<double>& E, const std::vector<double>& F, WRITEMODE mode = RECREATE) const;

  bool readmatrix(const std::string& filename, std::vector<double>& ETOA, std::vector<double>& ELIS, std::vector< std::vector<double> >& M, int ientry = 1);
  bool writematrix(const std::string& filename, const std::vector<double>& ETOA, const std::vector<double>& ELIS, const std::vector< std::vector<double> >& M, WRITEMODE mode = RECREATE) const;
};

class IO_BSON : public IO {
public:
  IO_BSON() {}
  ~IO_BSON() {}

  bool readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F, int ientry = 1);
  bool writespec(const std::string& filename, const std::vector<double>& E, const std::vector<double>& F, WRITEMODE mode = RECREATE) const;

  bool readmatrix(const std::string& filename, std::vector<double>& ETOA, std::vector<double>& ELIS, std::vector< std::vector<double> >& M, int ientry = 1);
  bool writematrix(const std::string& filename, const std::vector<double>& ETOA, const std::vector<double>& ELIS, const std::vector< std::vector<double> >& M, WRITEMODE mode = RECREATE) const;

};

#endif /* IO_H */