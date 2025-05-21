#include <fstream>
#include <iostream>
#include <cmath>

#include "Unit.h"
#include "rfl.hpp"
#include "rfl/json.hpp"
#include "rfl/bson.hpp"

#include "KDInterpSide.h"

using namespace std;
using namespace Unit;

vec_t operator*(const vec_t& a, const vec_t& b) {
  vec_t c(a.size());
  for (int i = 0; i < a.size(); i++)
    c[i] = a[i] * b[i];
  return c;
}

vec_t operator+(const vec_t& a, const vec_t& b) {
  vec_t c(a.size());
  for (int i = 0; i < a.size(); i++)
    c[i] = a[i] + b[i];
  return c;
}

vec_t operator-(const vec_t& a, const vec_t& b) {
  vec_t c(a.size());
  for (int i = 0; i < a.size(); i++)
    c[i] = a[i] - b[i];
  return c;
}

vec_t operator/(const vec_t& a, const vec_t& b) {
  vec_t c(a.size());
  for (int i = 0; i < a.size(); i++)
    c[i] = a[i] / b[i];
  return c;
}

ostream& operator<<(ostream& os, const vec_t& v) {
  for (auto& x : v)
    os << x << " ";
  return os;
}
ostream& operator<<(ostream& os, const KDPoint& p) {
  os << p.x << " | " << p.val;
  return os;
}

void KDPoint::show_neighbors() const {
  for (auto& n : neighbors)
    if (n == NULL) cout << "-----------------------------------" << endl;
    else cout << (*n) << endl;
}

bool KDPoint::connect(KDPoint* target, int ix) {
  for (int i = 0; i < x.size(); i++)
    if (i != ix && x[i] != target->x[i]) return false;

  assert(x[ix] != target->x[ix] && "There should not be two overlaping points");

  bool dir = x[ix] > target->x[ix];
  KDPoint* bound = target->neighbors[2*ix + dir];
  while (bound != NULL && (x[ix] - bound->x[ix]) * (x[ix] - target->x[ix]) > 0) {
    target = bound;
    bound = target->neighbors[2*ix + dir];
  }

  target->neighbors[2*ix + dir] = this;
  neighbors[2*ix + !dir] = target;

  neighbors[2*ix + dir] = bound;
  if (bound) {
    assert(x[ix] != bound->x[ix] && "There should not be two overlaping points");
    bound->neighbors[2*ix + !dir] = this;
  }
  return true;
}

std::set<KDPoint*> KDPoint::correction(const cfunc_t& correct) {
  std::set<KDPoint*> fresh_points;
  for (int ix = 0; ix < x.size(); ix++) {
    if (neighbors[2*ix] == NULL && neighbors[2*ix + 1] == NULL)
      continue;

    vector<KDPoint*> ps;
    if (neighbors[2*ix] != NULL) ps.push_back(neighbors[2*ix]);
    ps.push_back(this);
    if (neighbors[2*ix + 1] != NULL) ps.push_back(neighbors[2*ix + 1]);

    auto freshtmp = correct(ps);
    fresh_points.insert(freshtmp.begin(), freshtmp.end());

    for (auto& p : freshtmp) {
      auto subfresh = p->correction(correct);
      fresh_points.insert(subfresh.begin(), subfresh.end());
    }
  }

  return fresh_points;
}

std::set<KDValueSide*> KDPoint::kd_correction(const cfunc_t& correct) {
  std::set<KDValueSide*> res;
  auto fresh_points = correction(correct);

  for (auto& p : fresh_points)
    for (auto& val : p->blocks)
      res.insert(val);

  return res;
}

extern std::set<KDPoint*> null_func(const std::vector<KDPoint*>& points) { return std::set<KDPoint*>(); }

std::vector<vec_t> KDValueSide::get_sides(const vec_t& x, const vec_t& width) const {
  assert(x.size() == width.size());

  vector<vec_t> vp(x.size() * 2, x);
  for (int ix = 0; ix < x.size(); ix++) {
    vp[ix * 2][ix] = x[ix] - width[ix];
    vp[ix * 2 + 1][ix] = x[ix] + width[ix];
  }

  return vp;
}
std::vector<vec_t> KDValueSide::get_corners(const vec_t& x, const vec_t& width) const {
  assert(x.size() == width.size());
  assert(x.size() < 16 && "This version use int to implement the bounds, the maximum dimension suported is 16.");

  vec_t xmin = x - width,
        xmax = x + width;

  vector<vec_t> vp(interp->index.n);
  vector<bool> ivec(interp->index.dim);
  for (int i = 0; i < interp->index.n; i++) {
    interp->index.int2vec(i, ivec);
    for (int ix = 0; ix < interp->index.dim; ix++)
        vp[i].push_back(ivec[ix] ? xmax[ix] : xmin[ix]);
  }

  return vp;
}

KDPoint* KDValueSide::get_val(const vec_t& vx, const vector<KDPoint*>& ref_points) {
  auto i = interp->tab.find(vx);
  if (i != interp->tab.end())
    return i->second;

  auto p = new KDPoint(interp, vx, interp->func(interp->real_x(vx)), level);
  assert(ref_points.size() == interp->index.dim);
  for (int ix = 0; ix < interp->index.dim; ix++)
    if (ref_points[ix] != NULL)
      p->connect(ref_points[ix], ix);

  interp->update_tab(p);

  auto kds = p->kd_correction(interp->correction);
  for (auto& kd : kds)
    if (kd->complete) {
      kd->refresh_err();
      kd->breed();
    }

  return p;
};

double KDValueSide::eval(const vec_t& x) const {
  return dot(k, x) + c;
}

double KDValueSide::dot(const vec_t& a, const vec_t& b) {
  double sum = 0;
  for (int i = 0; i < a.size(); i++)
    sum += a[i] * b[i];
  return sum;
}

bool KDValueSide::init_exist_corners() {
  if (parent == NULL) return false;

  auto& slice = interp->index.slice[parent->ix_split][order];
  for (int is = 0; is < slice.size(); is++)
    corners[slice[is]] = parent->corners[slice[is]]; // obtain the corners from parent;

  auto brother = parent->children[!order];
  if (brother != NULL) {
    auto& antislice = interp->index.slice[parent->ix_split][!order];
    for (int is = 0; is < slice.size(); is++)
      corners[antislice[is]] = brother->corners[slice[is]]; // obtain the corners from brother;
  }
  return true;
}

bool KDValueSide::init_corners(const std::vector<vec_t>& vp) {
  corners.resize(interp->index.n, NULL);
  init_exist_corners();

  auto& antislice = interp->index.slice[parent->ix_split][!order];
  for (int i = 0; i < interp->index.n; i++) {
    if (corners[i] != NULL) continue;
    vector<KDPoint*> refs(interp->index.dim, NULL);

    if (parent != NULL)
      for (int ix_ref = 0; ix_ref < interp->index.dim; ix_ref++)
        refs[ix_ref] = corners[interp->index.slice_anti[ix_ref][i]];

    corners[i] = get_val(vp[i], refs);
  }

  for (auto& p : corners) p->blocks.push_back(this);

  return true;
}

bool KDValueSide::init_exist_sides() {
  if (parent == NULL) return false;

   // assign the overlaping side
  sides[parent->ix_split * 2 + order] = parent->sides[parent->ix_split * 2 + order];
  sides[parent->ix_split * 2 + !order] = parent->pmid;
  return true;
}

bool KDValueSide::init_sides(const std::vector<vec_t>& vp) {
  sides.resize(vp.size(), NULL);
  init_exist_sides();

  for (int i = 0; i < vp.size(); i++) {
    if (sides[i] != NULL) continue;
    vector<KDPoint*> refs(interp->index.dim, NULL);
    refs[i / 2] = pmid;

    sides[i] = get_val(vp[i], refs);
  }

  for (auto& p : sides) p->blocks.push_back(this);
  return true;
}

bool KDValueSide::init_points(const vec_t& x, const vec_t& width) {
  pmid = get_val(x, vector<KDPoint*>(interp->index.dim, NULL));
  //cout << x << " " << pmid->ix_split << endl;
  pmid->blocks.push_back(this);

  auto vcorner = get_corners(x, width);
  init_corners(vcorner);

  auto vside = get_sides(x, width);
  init_sides(vside);

  complete = true;
  return true;
}

KDValueSide::KDValueSide(KDInterpSide* interp_, const vec_t& x_, const vec_t& width_, int level_, int order_, KDValueSide* parent_)
 : interp(interp_), level(level_), ix_split(-1), order(order_), width(width_), parent(parent_), complete(false)
{
  init_points(x_, width_);

  refresh_err();
  breed();
}

KDValueSide::~KDValueSide() {
  for (auto& v : children) delete v;
}

bool KDValueSide::get_ix_split() {
  if (interp->read_mode) {
    ix_split = pmid->ix_split;
    return ix_split != -1;
  }

  int dim = interp->index.dim;
  int ix0 = parent == NULL ? 0 : parent->ix_split + 1;

  double wmax = 0;
  for (int ix = 0; ix < dim; ix++)
    wmax = fmax(wmax, width[ix]);

  for (int ix = ix0; ix < ix0 + dim; ix++) {
    ix_split = ix % dim;
    if (!interp->level_depths.empty() && (1 << interp->level_depths[ix_split]) * width[ix_split] > 1) return true;

    if (err[ix_split] * wmax > interp->tol) return true;
  }


  if (errmax * wmax > interp->tol) ix_split = ix0 % dim;
  else ix_split = -1;

  return ix_split != -1;
}

extern bool step_shape_p3(const std::vector<KDPoint*>& points, bool pflag);
extern bool distance_jump_side(const std::vector<KDPoint*>& points, bool pflag);
extern bool step_shape_side(const std::vector<KDPoint*>& points, bool pflag);
bool KDValueSide::breed() {
  if (!children.empty()) return false;
  if (!get_ix_split()) return false;
  if (level > 40) {
    cout << "level: " << parent->level << " " << parent->ix_split << endl;
    for (auto& p : parent->sides)
      cout << " -++ " << p->real_x()[0] / AU << " " << p->real_x()[1] / deg << " " << p->real_x()[2] / deg << " " << p->val / AU << endl;
    for (auto& p : parent->corners)
      cout << " -++ " << p->real_x()[0] / AU << " " << p->real_x()[1] / deg << " " << p->real_x()[2] / deg << " " << p->val / AU << endl;
 
    cout << "level: " << level << endl;
    for (auto& p : sides)
      cout << " --+ " << p->real_x()[0] / AU << " " << p->real_x()[1] / deg << " " << p->real_x()[2] / deg << " " << p->val / AU << endl;
    for (auto& p : corners)
      cout << " --+ " << p->real_x()[0] / AU << " " << p->real_x()[1] / deg << " " << p->real_x()[2] / deg << " " << p->val / AU << endl;
    cout << "ix_split: " << ix_split << endl;
    cout << "err: " << err[0] / AU << " " << err[1] / AU << " " << err[2] / AU << " " << errmax / AU << endl;
    cout << "width: " << width[0] << " " << width[1] << " " << width[2] << endl;
    bool dj = distance_jump_side(sides, true);
    bool ss = step_shape_side(sides, true);
    cout << "distance jump: " << dj << " step shape: " << ss << endl;
    exit(0);
  }

  pmid->ix_split = ix_split;
  children.resize(2, NULL);
  auto width_child = width;
  width_child[ix_split] /= 2;
  for (int i = 0; i < 2; i++) {
    auto xtmp = pmid->x;
    xtmp[ix_split] += width_child[ix_split] * (2 * i - 1);
    children[i] = new KDValueSide(interp, xtmp, width_child, level + 1, i, this);
  }
  return true;
}

void KDValueSide::refresh_err() {
  linear_eval();
  count_err();
}

void KDValueSide::linear_eval() {
  double vsum = 0;
  for (auto& p : sides)
    vsum += p->val;

  k.resize(interp->index.dim);

  for (int ix = 0; ix < interp->index.dim; ix++)
    k[ix] = (sides[ix * 2 + 1]->val - sides[ix * 2]->val) / (2 * width[ix]);

  c = vsum;
  for (auto& p : sides)
    c -= dot(k, p->x);
  c /= sides.size();
}

void KDValueSide::count_err() {
  err.resize(interp->index.dim);

  for (int ix = 0; ix < interp->index.dim; ix++)
    err[ix] = fabs(pmid->val - (sides[ix * 2 + 1]->val + sides[ix * 2]->val) / 2);

  errmax = 0;
  for (auto& p : corners)
    errmax = fmax(errmax, fabs(p->val - eval(p->x)));
}

double KDValueSide::operator()(const vec_t& x) const {
  return getkd(x)->eval(x);
}

const KDValueSide* KDValueSide::getkd(const vec_t& x) const {
  if (children.empty()) return this;

  if (x[ix_split] < pmid->x[ix_split]) return children[0]->getkd(x);

  return children[1]->getkd(x);
}

KDInterpSide::KDInterpSide(const std::string& tabfile) : read_mode(true) {
  ifstream ifs(tabfile, ios::binary);
  vector<char> bres((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
  KDMapSide result = rfl::bson::read<KDMapSide>(bres).value();
  xmid = result.xmid, width = result.width;
  tol = result.tol;
  level_depths = result.level_depths;
  index = NDIndex(xmid.size());
  ref_tab.resize(index.dim);

  assert(result.x.size() == xmid.size() && "The dimension of the table should be the same as the dimension of the function.");

  vector<double> xtmp(result.x.size());
  for (int i = 0; i < result.y.size(); i++) {
    for (int ix = 0; ix < result.x.size(); ix++)
      xtmp[ix] = result.x[ix][i];
    KDPoint* p = new KDPoint(this, xtmp, result.y[i], result.level[i]);
    p->ix_split = result.ix_split[i];
    tab.insert(pair<vec_t,KDPoint*>(p->x, p));
  }

  func = [&](const vec_t& x) -> double {
    auto relx = rel_x(x);
    cout << setprecision(16) << relx[0] << " " << relx[1] << " " << relx[2] << endl;
    auto iter = tab.upper_bound(relx);
    cout << setprecision(16) << iter->first[0] << " " << iter->first[1] << " " << iter->first[2] << endl;
    iter--;
    cout << setprecision(16) << iter->first[0] << " " << iter->first[1] << " " << iter->first[2] << endl;


    assert(false && "The function should not be called when initializing with exist table.");
    return 0;
  };

  correction = null_func;
  vec_t x(xmid.size(), 0),
        w(width.size(), 1);
  kd = new KDValueSide(this, x, w);
}

KDInterpSide::KDInterpSide(const func_t& func_, const vec_t& xmid_, const vec_t& width_, double tol_, const vector<int>& level_depths_, const cfunc_t& correction_) : xmid(xmid_), width(width_), tol(tol_), level_depths(level_depths_), read_mode(false), func(func_), correction(correction_), index(xmid.size())
{
  if (level_depths.empty()) level_depths.resize(xmid.size(), 0);
  ref_tab.resize(index.dim);

  vec_t x(xmid.size(), 0),
        w(width.size(), 1);
  kd = new KDValueSide(this, x, w);
}

KDInterpSide::~KDInterpSide() {
  delete kd;
  for (auto& v : tab) delete v.second;
}

vec_t KDPoint::real_x() const { return interp->real_x(x); }
vec_t KDInterpSide::real_x(const vec_t& x) const {
  return xmid + x * width;
}
vec_t KDInterpSide::rel_x(const vec_t& x) const {
  return (x - xmid) / width;
}

double KDInterpSide::operator()(const vec_t& x) const {
  return (*kd)(rel_x(x));
}

bool KDInterpSide::update_tab(KDPoint* p) {
  tab.insert(pair<vec_t,KDPoint*>(p->x, p));

  bool update_net = false;
  for (int ix = 0; ix < index.dim; ix++) {
    if (p->neighbors[2*ix] != NULL || p->neighbors[2*ix+1] != NULL) continue;

    vec_t list_x;
    for (int i = 0; i < index.dim; i++)
      if (i != ix) list_x.push_back(p->x[i]);
    auto ilist = ref_tab[ix].find(list_x);

    if (ilist == ref_tab[ix].end()) {
      ref_tab[ix].insert(pair<vec_t, KDPoint*>(list_x, p));
      update_net = true;
    } else p->connect(ilist->second, ix);
  }
  return update_net;
}

bool KDInterpSide::store_table(const std::string& filename, bool pflag) const {
  if (pflag)
    cout << ">> Storing " << tab.size() << " points to " << filename << endl;

  vector<vec_t> x;
  vector<double> y;
  vector<int> level;
  vector<int> ix_split;
  vector<int> nkd;
  vector<int> nlevels;

  x.resize(tab.begin()->first.size());
  for (auto& i : tab) {
    for (int ix = 0; ix < i.first.size(); ix++)
      x[ix].push_back(i.first[ix]);
    y.push_back(i.second->val);
    level.push_back(i.second->level);
    ix_split.push_back(i.second->ix_split);
    nkd.push_back(i.second->blocks.size());

    if (i.second->level >= nlevels.size()) nlevels.resize(i.second->level + 1, 0);
    nlevels[i.second->level]++;
  }

  const auto result = KDMapSide{.xmid=xmid, .width=width, .tol=tol, .level_depths=level_depths, .x = x, .y = y, .level = level, .ix_split=ix_split, .nkd = nkd};
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