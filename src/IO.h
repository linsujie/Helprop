#ifndef IO_H
#define IO_H

#include <map>
#include <vector>
#include <string>

class IO {
public:
  IO() {}
  ~IO() {}
  std::map<std::string, double> params;
  virtual bool readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F) const;
  virtual bool writespec(const std::string& filename, const std::vector<double>& E, const std::vector<double>& F) const;

  virtual bool readmatrix(const std::string& filename, const std::vector<double>& E, std::vector< std::vector<double> >& M) const;
  virtual bool writematrix(const std::string& filename, const std::vector<double>& E, const std::vector< std::vector<double> >& M) const;
};

class IO_TXT : public IO {
public:
  IO_TXT() {}
  ~IO_TXT() {}

  bool readspec(const std::string& filename, std::vector<double>& E, std::vector<double>& F) const;
  bool writespec(const std::string& filename, const std::vector<double>& E, const std::vector<double>& F) const;
};

#endif /* IO_H */