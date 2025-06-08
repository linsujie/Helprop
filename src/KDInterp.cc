#include <fstream>
#include <iostream>
#include <filesystem>
#include <cmath>

#include "rfl.hpp"
#include "rfl/json.hpp"
#include "rfl/bson.hpp"

#include "KDInterp.h"

using namespace std;

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
      if (p == this) continue;
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

bool verify_kdvalue(KDValueSide* kd) {
  KDValueSide* parent = kd->parent;
  while (true) {
    if (parent == NULL && kd->level == 0) return true;

    if (parent->children.empty() || parent->ix_split < 0 || parent->children[kd->order] != kd) {
      cout << "chain break at " << kd << " " << kd->level << " <-> " << parent << " " << parent->level << endl;
      return false;
    }

    kd = parent;
    parent = kd->parent;
  }

  cout << "verify_kdvalue failed" << endl;
  return false;
}

KDValueSide* KDValueSide::front_offspring() {
  KDValueSide *front = this,
              *cur = this;
  while (cur->parent != NULL) {
    if (cur->parent->alive && !cur->alive) front = cur->parent;
    cur = cur->parent;
  }
  return front;
}

KDPoint* KDValueSide::get_val(const vec_t& vx, const vector<KDPoint*>& ref_points) {
  if (interp->read_mode) {
    KDPoint* p = interp->points[interp->orders[interp->icurr]];
    //cout << "reading: " << interp->icurr << " " << interp->orders[interp->icurr] << " " << vx << " " << p->x << " " << p->ix_split << endl;
    assert(vx == p->x);
    interp->icurr++;

    return p;
  }

  auto i = interp->tab.find(vx);
  if (i != interp->tab.end()) return i->second;

  auto p = new KDPoint(interp, vx, interp->func(interp->real_x(vx)), level);
  assert(ref_points.size() == interp->index.dim);
  for (int ix = 0; ix < interp->index.dim; ix++)
    if (ref_points[ix] != NULL)
      p->connect(ref_points[ix], ix);

  interp->update_tab(p);

  auto kds = p->kd_correction(interp->correction);

  for (auto& kd : kds) {
    if (kd == this || !kd->complete || !kd->alive) continue;
    interp->refresh_blocks.insert(kd);
  }
  return p;
};

bool KDValueSide::is_ancestor_of(const KDValueSide* v) const {
  while (v->parent != NULL) {
    if (v->parent == this) return true;
    v = v->parent;
  }

  return false;
}

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
  if (!interp->read_mode) init_exist_corners();

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
  if (!interp->read_mode) init_exist_sides();

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

KDValueSide::KDValueSide(KDInterp* interp_, const vec_t& x_, const vec_t& width_, int level_, int order_, KDValueSide* parent_)
 : interp(interp_), level(level_), ix_split(-1), order(order_), alive(true), width(width_), parent(parent_), complete(false)
{
  init_points(x_, width_);

  refresh_err();
}

KDValueSide::~KDValueSide() {
  for (auto& v : children) delete v;
}

bool KDValueSide::get_ix_split(int& ix_split_new) const {
  if (interp->read_mode) {
    ix_split_new = pmid->ix_split;
    return ix_split_new != ix_split;
  }

  int dim = interp->index.dim;
  int ix0 = parent == NULL ? interp->ix_split0 : parent->ix_split + 1;

  vec_t tol_local = interp->norm_tol;
  double tmax = 0, tmin = 1e300;
  int ix_min = -1;
  for (int ix = 0; ix < dim; ix++) {
    tol_local[ix] /= width[ix]; // smaller grid are allowed to have larger tolerance
    if (tmax < tol_local[ix]) {
      tmax = tol_local[ix];
    }

    if (tmin > tol_local[ix]) {
      tmin = tol_local[ix];
      ix_min = ix;
    }
  }
  interp->n_min_tol[ix_min]++;

  for (int ix = ix0; ix < ix0 + dim; ix++) {
    ix_split_new = ix % dim;
    if (!interp->level_depths.empty() && (1 << interp->level_depths[ix_split_new]) * width[ix_split_new] > 1) return ix_split_new != ix_split;

    if (err[ix_split_new] > tol_local[ix_split_new]) return  ix_split_new != ix_split;
  }

  if (errmax > tmin) ix_split_new = ix0 % dim;
  else ix_split_new = -1;

  return  ix_split_new != ix_split;
}

extern bool step_shape_p3(const std::vector<KDPoint*>& points, bool pflag);
extern bool distance_jump_side(const std::vector<KDPoint*>& points, bool pflag);
extern bool step_shape_side(const std::vector<KDPoint*>& points, bool pflag);

bool kill(set<KDValueSide*>& dead_list, KDValueSide* kd) {
  if (!kd || !kd->alive) return false;

  kd->alive = false;
  kd->pmid->ix_split = -1;
  dead_list.insert(kd);
  for (auto& c : kd->children)
    if (c && c->alive) kill(dead_list, c);

  return true;
}

bool KDValueSide::breed() {
  int ix_split_new;
  if (!get_ix_split(ix_split_new))
    return false; // if ix_split unchanged, no need to breed
  assert(children.empty());

  ix_split = ix_split_new;
  pmid->ix_split = ix_split;
  if (ix_split == -1) return false;

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

bool pass_through(list<KDValueSide*>& active_blocks, const function<void(list<KDValueSide*>&,KDValueSide*)>& func) {
  int level_min = 9999;
  for (const auto& b : active_blocks)
    level_min = fmin(level_min, b->level);

  list<KDValueSide*> new_active_blocks;
  for (auto iter = active_blocks.begin(); iter != active_blocks.end();) {
    KDValueSide *b = *iter;
   if (b->level == level_min) {
      func(new_active_blocks, b);
      iter = active_blocks.erase(iter);
    } else iter++;
  }
  active_blocks.splice(active_blocks.end(), new_active_blocks);
  return !active_blocks.empty();
}

KDInterp::KDInterp(const std::string& tabfile) : read_mode(true) {
  assert(filesystem::exists(tabfile));

  ifstream ifs(tabfile, ios::binary);
  vector<char> bres((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
  KDMapSide result = rfl::bson::read<KDMapSide>(bres).value();
  xmid = result.xmid, width = result.width;
  tol = result.tol;
  level_depths = result.level_depths;
  orders = result.orders;
  index = NDIndex(xmid.size());
  ref_tab.resize(index.dim);
  n_min_tol.resize(index.dim, 0);

  assert(result.x.size() == xmid.size() && "The dimension of the table should be the same as the dimension of the function.");

  points.reserve(result.y.size());
  vector<double> xtmp(result.x.size());
  for (int i = 0; i < result.y.size(); i++) {
    for (int ix = 0; ix < result.x.size(); ix++)
      xtmp[ix] = result.x[ix][i];

    KDPoint* p = new KDPoint(this, xtmp, result.y[i], result.level[i]);
    p->ix_split = result.ix_split[i];
    points.push_back(p);
  }

  func = [&](const vec_t& x) -> double {
    auto relx = rel_x(x);
    cout << setprecision(16) << relx[0] << " " << relx[1] << " " << relx[2] << " " <<  relx[3] << endl;
    auto iter = tab.upper_bound(relx);
    cout << setprecision(16) << iter->first[0] << " " << iter->first[1] << " " << iter->first[2] << endl;
    iter--;
    cout << setprecision(16) << iter->first[0] << " " << iter->first[1] << " " << iter->first[2] << endl;


    assert(false && "The function should not be called when initializing with exist table.");
    return 0;
  };

  norm_tol = tol;
  for (int i = 0; i < norm_tol.size(); i++)
    norm_tol[i] /= width[i];

  correction = null_func;
  vec_t x(xmid.size(), 0),
        w(width.size(), 1);
  icurr = 0;
  kd = new KDValueSide(this, x, w);
  spring();
}

KDInterp::KDInterp(const func_t& func_, const vec_t& xmid_, const vec_t& width_, const vec_t& tol_, const vector<int>& level_depths_, int ix_split0_, const cfunc_t& correction_) : xmid(xmid_), width(width_), tol(tol_), level_depths(level_depths_), ix_split0(ix_split0_), read_mode(false), func(func_), correction(correction_), index(xmid.size())
{
  assert(xmid.size() == width.size() && width.size() == tol.size());
  if (level_depths.empty()) level_depths.resize(xmid.size(), 0);
  ref_tab.resize(index.dim);

  norm_tol = tol;
  for (int i = 0; i < norm_tol.size(); i++)
    norm_tol[i] /= width[i];
  n_min_tol.resize(index.dim, 0);

  vec_t x(xmid.size(), 0),
        w(width.size(), 1);
  kd = new KDValueSide(this, x, w);
  spring();
  cout << ">> KDInterp: has " << dead_children.size() << " dead children." << endl;

  list<KDValueSide*> active_blocks = { kd };
  set<KDValueSide*> end_blocks;

  while (!active_blocks.empty()) {
    pass_through(active_blocks,
    [&](list<KDValueSide*>& new_active_blocks, KDValueSide* b) {
      int ix_split_new;
      if (b->get_ix_split(ix_split_new)) {
        cout << "block: " << b << " at " << b->level << " seems change from " << b->ix_split << " to " << ix_split_new << endl;
        if (ix_split_new == -1) end_blocks.insert(b);
      }

      for (auto& c : b->children)
        new_active_blocks.push_back(c);
    });
  }
  cout << "There are " << end_blocks.size() << " end blocks unfound." << endl;
}

bool KDInterp::spring() {
  int ix_split_new;
  list<KDValueSide*> active_blocks = { kd };

  while (!active_blocks.empty()) {
    pass_through(active_blocks,
    [](list<KDValueSide*>& new_active_blocks, KDValueSide* b) {
      if (b->breed()) // if there are new borned, adding them to the active blocks.
        for (auto& c : b->children)
          new_active_blocks.push_back(c);
    });

    if (read_mode) continue;
  
    for (auto& b : active_blocks)
      refresh_blocks.erase(b); // to avoid the multiple dealing of the active blocks
  
    // To active the refreshed blocks that has changed their status and kill their origin children.
    for (auto& b : refresh_blocks) {
      b->refresh_err();
      if (b->get_ix_split(ix_split_new)) {
        active_blocks.push_back(b);
        for (auto& c : b->children)
          kill(dead_children, c);
        b->children.clear();
      }
    }
  
    refresh_blocks.clear();
  
    // To ensure all the active blocks are alive.
    for (auto iter = active_blocks.begin(); iter != active_blocks.end();)
      if (!(*iter)->alive) iter = active_blocks.erase(iter);
      else iter++;
  }

 return true;     
}

KDInterp::~KDInterp() {
  delete kd;
  for (auto& v : tab) delete v.second;
  for (auto& v : points) delete v;
}

vec_t KDPoint::real_x() const { return interp->real_x(x); }
vec_t KDInterp::real_x(const vec_t& x) const {
  return xmid + x * width;
}
vec_t KDInterp::rel_x(const vec_t& x) const {
  return (x - xmid) / width;
}

double KDInterp::operator()(const vec_t& x) const {
  return (*kd)(rel_x(x));
}

bool KDInterp::update_tab(KDPoint* p) {
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

inline void push_back(vector<int>& orders, vector<KDPoint*>& points, KDPoint* p) {
  if (p->store_order < 0) {
    points.push_back(p);
    p->store_order = points.size() - 1;
  }
  orders.push_back(p->store_order);
}

void add_points(vector<int>& orders, vector<KDPoint*>& points, const KDValueSide* kd) {
  push_back(orders, points, kd->pmid);
  for (auto& p : kd->corners) push_back(orders, points, p);
  for (auto& p : kd->sides) push_back(orders, points, p);
}

bool KDInterp::store_table(const std::string& filename, bool pflag) const {
  if (pflag)
    cout << ">> Storing to " << filename << endl;

  vector<KDPoint*> points;
  vector<int> orders;
  points.reserve(tab.size());
  orders.reserve(tab.size());

  list<KDValueSide*> active_blocks = { kd };
  while (!active_blocks.empty()) {
    pass_through(active_blocks,
    [&](list<KDValueSide*>& new_active_blocks, KDValueSide* b) {
      add_points(orders, points, b);

      for (auto& c : b->children)
        new_active_blocks.push_back(c);
    });
  }

  cout << ">> To store "  << points.size() << " of " << tab.size() << " points." << endl;

  vector<vec_t> x;
  vector<double> y;
  vector<int> level;
  vector<int> ix_split;
  vector<int> nkd;
  vector<int> nlevels;

  x.resize((*points.begin())->x.size());
  for (auto& i : points) {
    for (int ix = 0; ix < i->x.size(); ix++)
      x[ix].push_back(i->x[ix]);
    y.push_back(i->val);
    level.push_back(i->level);
    ix_split.push_back(i->ix_split);
    nkd.push_back(i->blocks.size());

    if (i->level >= nlevels.size()) nlevels.resize(i->level + 1, 0);
    nlevels[i->level]++;
  }

  const auto result = KDMapSide{.xmid=xmid, .width=width, .tol=tol, .orders = orders, .level_depths=level_depths, .x = x, .y = y, .level = level, .ix_split=ix_split, .nkd = nkd};
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

void count_nlevels(const KDValueSide* kd, int l, vector<int>& nlevels) {
  nlevels.resize(fmax(nlevels.size(), l + 1), 0);
  nlevels[l]++;

  for (auto& i : kd->children)
    count_nlevels(i, l + 1, nlevels);
}

void KDValueSide::show() const {
  cout << "KDValues: " << endl;
  vector<int> nlevels;

  count_nlevels(this, 0, nlevels);
  for (auto& n : nlevels) cout << " " << n;
  cout << endl;
}

void KDInterp::show() const {
  kd->show();
}

void compare_one(const KDValueSide& v1, const KDValueSide& v2, const vec_t& x1, const vec_t& x2, const vec_t& w1, const vec_t& w2) {
  if ((v1.pmid->x - x1) / w1 != (v2.pmid->x - x2) / w2) {
    cout << "  Subblock different at: " << (v1.pmid->x - x1) / w1 << " <-> " << (v2.pmid->x - x2) / w2 << endl;
    return;
  }

  if (v1.ix_split != v2.ix_split) {
    cout << "ix_split different at: " << (v1.pmid->x - x1) / w1 << " | " << v1.ix_split << " <-> " << v2.ix_split << endl;
  }

  if (v1.ix_split != -1 && v2.ix_split != -1)
    for (int i = 0; i < 2; i++)
      compare_one(*v1.children[i], *v2.children[i], x1, x2, w1, w2);
}

void compare(const KDValueSide& v1, const KDValueSide& v2) {
  auto& w1 = v1.width;
  auto& x1 = v1.pmid->x;
  auto& w2 = v2.width;
  auto& x2 = v2.pmid->x;

  compare_one(v1, v2, x1, x2, w1, w2);
}