#ifndef KDINTERP_H
#define KDINTERP_H
#include<vector>
#include<map>
#include<set>
#include<functional>
#include<cassert>
#include <string>

#include"NDIndex.h"

typedef std::vector<double> vec_t;
typedef std::function<double(const vec_t&)> func_t;
struct vector_less_than {
  bool operator() (const vec_t& lhs, const vec_t& rhs) const {
    assert(lhs.size() == rhs.size());

    for (int i = 0; i < lhs.size(); i++)
      if (lhs[i] > rhs[i]) return false; // more than
      else if (lhs[i] < rhs[i]) return true; // less than

    return false; // equal
  }
};

class KDValue;
struct KDPoint {
  vec_t x;
  double val;
  int level;
  long index;
  std::vector<KDValue*> blocks;

  KDPoint(const vec_t& x_, double val_, int level_, long index_) : x(x_), val(val_), level(level_), index(index_) {}
};
typedef std::map<vec_t, KDPoint*, vector_less_than> tab_t;
typedef std::function<std::set<KDPoint*>(std::vector<KDPoint*>& func)> cfunc_t;
std::set<KDPoint*> null_func(const std::vector<KDPoint*>& points);

class KDValue {
  private:
    std::vector<vec_t> get_bounds(const vec_t& xmin, const vec_t& xmax) const;
    KDPoint* get_val(const func_t& func, const vec_t& x, int ind);
    double eval(const vec_t& x) const;

    static double dot(const vec_t& a, const vec_t& b);
    static void sort_minmax(vec_t& xmin, vec_t& xmax);

    bool init_exist_points();
    bool init_points(const func_t& func, const std::vector<vec_t>& bound);
    std::set<KDValue*> kd_correction(const cfunc_t& correction);
    bool breed(const func_t& func, const cfunc_t& correction, int level_depth, double tol);
    void refresh_err();

  public:

  KDValue(const func_t& func, const vec_t& xmin_, const vec_t& xmax_, tab_t& tab_, double tol, const NDIndex& index_,
   const cfunc_t& correction = null_func,
   int level_depth = 0, int level_ = 0, int order_ = -1, KDValue* parent_ = NULL);
  ~KDValue();

  void linear_eval();
  double count_err() const;

  double operator()(const vec_t& x) const;
  const KDValue* getkd(const vec_t& x) const;

  vec_t k;
  double c, err;

  NDIndex index;
  int level;
  int order;
  long kd_index, kd_n;
  KDValue* parent;
  std::vector<KDValue*> children;

  std::vector<KDPoint*> points;
  tab_t& tab;
};

struct KDMap{
  vec_t xmin, xmax;
  double tol;
  int level_depth;
  std::vector<vec_t> x;
  std::vector<double> y;
  std::vector<int> level;
  std::vector<long> index;
  std::vector<int> nkd;
};

class KDInterp {
  public:

  KDInterp(const std::string& tabfile);
  KDInterp(const func_t& func, const vec_t& xmin_, const vec_t& xmax_, double tol_ = 0.01, int level_depth_ = 0, const cfunc_t& correction = null_func);
  ~KDInterp();
  double operator()(const vec_t& x) const;
  bool store_table(const std::string& filename, bool pflag = true) const;

  vec_t xmin, xmax;
  double tol;
  int level_depth;
  tab_t tab;
  KDValue *kd;
};

#endif /* KDINTERP_H */