#include <fstream>
#include <iostream>
#include <cmath>

#include "rfl.hpp"
#include "rfl/json.hpp"
#include "rfl/bson.hpp"

#include "KDInterp.h"

using namespace std;
extern std::set<KDPoint*> null_func(const std::vector<KDPoint*>& points) { return std::set<KDPoint*>(); }

std::vector<vec_t> KDValue::get_bounds(const vec_t& xmin, const vec_t& xmax) const {
  assert(xmin.size() == xmax.size());
  assert(xmin.size() < 16 && "This version use int to implement the bounds, the maximum dimension suported is 16.");
  vector<vec_t> bounds(index.n);
  vector<bool> ivec(index.dim);
  for (int i = 0; i < index.n; i++) {
    index.int2vec(i, ivec);
    for (int ix = 0; ix < index.dim; ix++)
        bounds[i].push_back(ivec[ix] ? xmax[ix] : xmin[ix]);
  }

  return bounds;
}

KDPoint* KDValue::get_val(const func_t& func, const vec_t& x, int ind) {
  auto i = tab.find(x);
  if (i != tab.end())
    return i->second;

  auto p = new KDPoint(x, func(x), level, kd_index * (index.n + 1) + ind);
  tab.insert(pair<vec_t,KDPoint*>(x, p));
  return p;
};

double KDValue::eval(const vec_t& x) const {
  return dot(k, x) + c;
}

double KDValue::dot(const vec_t& a, const vec_t& b) {
  double sum = 0;
  for (int i = 0; i < a.size(); i++)
    sum += a[i] * b[i];
  return sum;
}

void KDValue::sort_minmax(vec_t& xmin, vec_t& xmax) {
  for (int ix = 0; ix < xmin.size(); ix++)
    if (xmin[ix] > xmax[ix]) {
      double tmp = xmin[ix];
      xmin[ix] = xmax[ix];
      xmax[ix] = tmp;
    }
}

bool KDValue::init_exist_points() {
  if (parent == NULL) return false;

  points[order] = parent->points[order]; // assign the corner of parent
  points[index.anti[order]] = parent->points[index.n]; // assign the midpoint of parent

  vector<bool> ivec(index.dim);
  index.int2vec(order, ivec);
  for (int ix = 0; ix < index.dim; ix++) {
    auto neighbor_kd = parent->children[index.slice_anti[ix][order]];

    if (neighbor_kd == NULL) continue;

    for (auto ip : index.slice[ix][1 - ivec[ix]]) {
      if (points[ip] != NULL)
       assert(points[ip] == neighbor_kd->points[index.slice_anti[ix][ip]] && "The overlaping point should be the same.");
      points[ip] = neighbor_kd->points[index.slice_anti[ix][ip]];
    }
  }

  return true;
}

bool KDValue::init_points(const func_t& func, const std::vector<vec_t>& bound) {
  points.resize(bound.size() + 1, NULL);
  //cout << "init_val: " << level << " " << order << " " << kd_index * (index.n + 1) + 0 << endl;

  init_exist_points();

  for (int i = 0; i < bound.size(); i++) {
    if (points[i] != NULL) continue;

    points[i] = get_val(func, bound[i], i + 1);
  }

  vec_t xmid;
  for (int ix = 0; ix < index.dim; ix++)
     xmid.push_back((bound[0][ix] + bound[index.n - 1][ix]) / 2);
  points[index.n] = get_val(func, xmid, 0);

  for (auto& p : points)
    p->blocks.push_back(this);
  return true;
}

set<KDValue*> KDValue::kd_correction(const cfunc_t& correction) {
  auto fresh_points = correction(points);

  set<KDValue*> correction_set;
  
  for (auto& p : fresh_points) {
    for (auto& k : p->blocks)
      if (k != this)
        correction_set.insert(k);
  }

  set<KDValue*> total_correction_set = correction_set;
  for (auto& k : correction_set) {
    auto subset = k->kd_correction(correction);
    total_correction_set.insert(subset.begin(), subset.end());
  }

  total_correction_set.insert(this);
  return total_correction_set;
}
  
KDValue::KDValue(const func_t& func, const vec_t& xmin, const vec_t& xmax, tab_t& tab_, double tol, const NDIndex& index_, const cfunc_t& correction, int level_depth, int level_, int order_, KDValue* parent_)
 : index(index_), level(level_), order(order_), parent(parent_), tab(tab_)
{
  if (parent == NULL) {
    kd_n = 1;
    kd_index = 0;
  } else {
    kd_n = parent->kd_n * (index.n + 1);
    kd_index = parent->kd_index + parent->kd_n * (order + 1);
  }

  auto b = get_bounds(xmin, xmax);
  init_points(func, b);
  auto correction_set = kd_correction(correction);

  //if (correction_set.size() > 1 && level > 10) 
  //  cout << "refreshing number: " << correction_set.size() << endl;
  refresh_err();
  for (auto& k : correction_set) {
    if (k == this) continue;
    double errold = k->err;
    k->refresh_err();
    k->breed(func, correction, level_depth, tol);
    //if (parent == NULL || parent->parent == NULL) continue;

    //if (distance_jump(parent->parent->points) || step_shape(parent->parent->points)) {
    //  distance_jump(parent->parent->points, true);
    //  step_shape(parent->parent->points, true);
    //  for (auto& p : parent->parent->points)
    //    cout << " +-+ " << p->x[0] << " " << p->x[1] << " " << p->x[2] << " " << p->val << endl;
    //  exit(0);
    //}
  }

  breed(func, correction, level_depth, tol);
}

KDValue::~KDValue() {
  for (auto& v : children) delete v;
}

bool KDValue::breed(const func_t& func, const cfunc_t& correction, int level_depth, double tol) {
  if (!children.empty()) return false;

  if (level >= level_depth
      && err < tol * (1 << level)
      ) return false;

  children.resize(index.n, NULL);
  for (int i = 0; i < index.n; i++) {
    vector<double> xmin_c = points[i]->x, xmax_c = points[index.n]->x;
    sort_minmax(xmin_c, xmax_c);
    children[i] = new KDValue(func, xmin_c, xmax_c, tab, tol, index, correction, level_depth, level + 1, i, this);
  }
  return true;
}

void KDValue::refresh_err() {
  linear_eval();
  err = count_err();
}

void KDValue::linear_eval() {
  double vsum = 0;
  for (auto& p : points)
    vsum += p->val;

  k.resize(index.dim);

  for (int ix = 0; ix < index.dim; ix++) {
    double vdiff = 0;

    for (auto i : index.slice[ix][1]) vdiff += points[i]->val;
    for (auto i : index.slice[ix][0]) vdiff -= points[i]->val;

    k[ix] = vdiff / (index.n / 2) / (points[index.n - 1]->x[ix] - points[0]->x[ix]);
  }

  c = vsum;
  for (auto& p : points)
    c -= dot(k, p->x);
  c /= points.size();
}

double KDValue::count_err() const {
  double res = 0;
  for (auto& p : points)
      res += (p->val - this->eval(p->x)) * (p->val - this->eval(p->x));

  res = sqrt(res / points.size());
  KDPoint* midp = points[points.size() - 1];
  return fmax(res, fabs(midp->val - this->eval(midp->x)));
}

double KDValue::operator()(const vec_t& x) const {
  return getkd(x)->eval(x);
}

const KDValue* KDValue::getkd(const vec_t& x) const {
  if (children.empty()) return this;

  int i = 0;
  for (int ix = 0; ix < index.dim; ix++)
    if (x[ix] > points[index.n]->x[ix]) i |= (1 << ix);
  return children[i]->getkd(x);
}

KDInterp::KDInterp(const std::string& tabfile) {
  ifstream ifs(tabfile, ios::binary);
  vector<char> bres((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
  KDMap result = rfl::bson::read<KDMap>(bres).value();
  xmin = result.xmin, xmax = result.xmax;
  tol = result.tol;
  level_depth = result.level_depth;

  assert(result.x.size() == xmin.size() && "The dimension of the table should be the same as the dimension of the function.");

  vector<double> xtmp(result.x.size());
  for (int i = 0; i < result.y.size(); i++) {
    for (int ix = 0; ix < result.x.size(); ix++)
      xtmp[ix] = result.x[ix][i];
    tab.insert(pair<vec_t,KDPoint*>(xtmp, new KDPoint(xtmp, result.y[i], result.level[i], result.index[i])));
  }

  NDIndex index(xmin.size());
  auto func = [&](const vec_t& x) -> double {
    assert(false && "The function should not be called when initializing with exist table.");
    return 0;
  };
  kd = new KDValue(func, xmin, xmax, tab, tol, index, null_func, level_depth);
}

KDInterp::KDInterp(const func_t& func, const vec_t& xmin_, const vec_t& xmax_, double tol_, int level_depth_, const cfunc_t& correction) : xmin(xmin_), xmax(xmax_), tol(tol_), level_depth(level_depth_)
{
  NDIndex index(xmin.size());
  kd = new KDValue(func, xmin, xmax, tab, tol, index, correction, level_depth);
}

KDInterp::~KDInterp() {
  delete kd;
  for (auto& v : tab) delete v.second;
}

double KDInterp::operator()(const vec_t& x) const {
  return (*kd)(x);
}

bool KDInterp::store_table(const std::string& filename, bool pflag) const {
  if (pflag)
    cout << ">> Storing " << tab.size() << " points to " << filename << endl;
  vector<vector<double> > x;
  vector<double> y;
  vector<int> level;
  vector<int> nkd;
  vector<long> ind;
  vector<int> nlevels;

  x.resize(tab.begin()->first.size());
  for (auto& i : tab) {
    for (int ix = 0; ix < i.first.size(); ix++)
      x[ix].push_back(i.first[ix]);
    y.push_back(i.second->val);
    level.push_back(i.second->level);
    ind.push_back(i.second->index);
    nkd.push_back(i.second->blocks.size());

    if (i.second->level >= nlevels.size()) nlevels.resize(i.second->level + 1, 0);
    nlevels[i.second->level]++;
  }

  const auto result = KDMap{.xmin=xmin, .xmax=xmax, .tol=tol, .level_depth=level_depth, .x = x, .y = y, .level = level, .index = ind, .nkd = nkd};
  vector<char> bres = rfl::bson::write(result);
  FILE *of = fopen(filename.c_str(), "w");
  fwrite(&bres[0], 1, bres.size(), of);
  fclose(of);

  if (pflag) {
    cout << "<< Stored:";
    for (auto& n : nlevels) cout << " " << n;
    cout << endl;
  }
  return true;
}