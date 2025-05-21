#ifndef KDINTERPSIDE_H
#define KDINTERPSIDE_H
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

vec_t operator*(const vec_t& a, const vec_t& b);
vec_t operator+(const vec_t& a, const vec_t& b);
vec_t operator-(const vec_t& a, const vec_t& b);
vec_t operator/(const vec_t& a, const vec_t& b);

class KDValueSide;
class KDInterpSide;
struct KDPoint;
typedef std::function<std::set<KDPoint*>(std::vector<KDPoint*>& func)> cfunc_t;
typedef std::map<vec_t, KDPoint*, vector_less_than> tab_t;
std::set<KDPoint*> null_func(const std::vector<KDPoint*>& points);

struct KDPoint {
  vec_t x;
  double val;
  int level;
  int ix_split;
  KDInterpSide* interp;
  std::vector<KDPoint*> neighbors;
  std::vector<KDValueSide*> blocks;

  KDPoint(KDInterpSide* interp_, const vec_t& x_, double val_, int level_) : interp(interp_), x(x_), val(val_), level(level_), ix_split(-1), neighbors(2 * x.size(), NULL) {}
  vec_t real_x() const;

  bool connect(KDPoint* ref, int ix);
  std::set<KDPoint*> correction(const cfunc_t& correct);
  std::set<KDValueSide*> kd_correction(const cfunc_t& correct);
  void show_neighbors() const;
};
std::ostream& operator<<(std::ostream& os, const KDPoint& p);

class KDValueSide {
  private:
    std::vector<vec_t> get_corners(const vec_t& x, const vec_t& width) const;
    std::vector<vec_t> get_sides(const vec_t& x, const vec_t& width) const;
    KDPoint* get_val(const vec_t& x, const std::vector<KDPoint*>& ref_points);
    double eval(const vec_t& x) const;

    static double dot(const vec_t& a, const vec_t& b);

    bool init_exist_corners();
    bool init_corners(const std::vector<vec_t>& vp);
    bool init_exist_sides();
    bool init_sides(const std::vector<vec_t>& vp);
    bool init_points(const vec_t& x, const vec_t& width);

    bool breed();
    void refresh_err();

  public:

  KDValueSide(KDInterpSide* interp_, const vec_t& x_, const vec_t& width_,
   int level_ = 0, int order_ = -1, KDValueSide* parent_ = NULL);
  ~KDValueSide();

  void linear_eval();
  void count_err();
  bool get_ix_split();

  double operator()(const vec_t& x) const;
  const KDValueSide* getkd(const vec_t& x) const;

  vec_t k, err;
  double c, errmax;

  int level;
  int ix_split, order;
  vec_t width;
  KDValueSide* parent;
  KDInterpSide* interp;
  std::vector<KDValueSide*> children;

  std::vector<KDPoint*> corners, sides;
  KDPoint* pmid;
  bool complete;
};

struct KDMapSide {
  vec_t xmid, width;
  double tol;
  std::vector<int> level_depths;
  std::vector<vec_t> x;
  std::vector<double> y;
  std::vector<int> level;
  std::vector<int> ix_split;
  std::vector<int> nkd;
};

class KDInterpSide {
  public:

  KDInterpSide(const std::string& tabfile);
  KDInterpSide(const func_t& func, const vec_t& xmid_, const vec_t& width_, double tol_ = 0.01, const std::vector<int>& level_depths_ = {}, const cfunc_t& correction = null_func);
  ~KDInterpSide();
  double operator()(const vec_t& x) const;
  bool update_tab(KDPoint* p);
  bool store_table(const std::string& filename, bool pflag = true) const;
  vec_t real_x(const vec_t& x) const;
  vec_t rel_x(const vec_t& x) const;

  vec_t xmid, width;
  double tol;
  std::vector<int> level_depths;
  tab_t tab;
  std::vector<tab_t> ref_tab;
  bool read_mode;
  func_t func;
  cfunc_t correction;
  NDIndex index;
  KDValueSide *kd;
};
#endif /* KDINTERPSIDE_H */
