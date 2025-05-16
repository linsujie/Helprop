#include <map>
#include <vector>
#include <iomanip>
#include "hcs_interp.h"
#include "Unit.h"
#include "Vec.hh"

using namespace std;
using namespace Unit;

KDPoint* extrema(const std::vector<KDPoint*>& p, int low, int up, const std::function<bool(double,double)>& cmp) {
  KDPoint* ptmp = p[low];
  if (up < 0) up += p.size();
  for (int i = low + 1; i <= up; i++)
    if (cmp(p[i]->val, ptmp->val))
      ptmp= p[i];

  return ptmp;  
}

bool lessthan(double a, double b) { return a < b; }
KDPoint* vmin(const std::vector<KDPoint*>& p, int low = 0, int up = -1) {
  return extrema(p, low, up, lessthan);
}

bool absless(double a, double b) { return fabs(a) < fabs(b); }
KDPoint* absmin(const std::vector<KDPoint*>& p, int low = 0, int up = -1) {
  return extrema(p, low, up, absless);
}

bool greaterthan(double a, double b) { return a > b; }
KDPoint* vmax(const std::vector<KDPoint*>& p, int low = 0, int up = -1) {
  return extrema(p, low, up, greaterthan);
}

bool absgreater(double a, double b) { return fabs(a) > fabs(b); }
KDPoint* absmax(const std::vector<KDPoint*>& p, int low = 0, int up = -1) {
  return extrema(p, low, up, absgreater);
}

bool sameside(const std::vector<KDPoint*>& points, bool pflag = false) {
  for (int i = 1; i < points.size(); i++)
    if (points[i]->val * points[0]->val < 0) return false;

  return true;
}

double mid(const vector<double>& v) {
  double xmin = v[0],
  xmax = v[0];

  for (auto& x : v) {
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
  }

  return (xmin + xmax) / 2;
}

double width(const vector<double>& v) {
  double xmin = v[0],
  xmax = v[0];

  for (auto& x : v) {
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
  }

  return (xmax - xmin) / 2;
}

bool distance_jump(std::vector<KDPoint*> points, bool pflag = false) {
  auto pmin = absmin(points),
   pmax = absmax(points);
  double vgap = fabs(pmax->val - pmin->val);

  Vec p1, p2;
  p1.set_spherical(pmin->x[0], pmin->x[1] * deg, pmin->x[2] * deg);
  p2.set_spherical(pmax->x[0], pmax->x[1] * deg, pmax->x[2] * deg);

  if (pflag)
    cout << "vgap: " << vgap << " distance: " << (p1 - p2).len() << endl;
  return (p1 - p2).len() < vgap;
}

bool within(double a, double l, double u) { return l < a && a < u; }
bool step_shape(const std::vector<KDPoint*>& points, bool pflag = false) {
  if (!sameside(points, pflag)) return false;

  KDPoint* pmin = vmin(points, 0, -2),
   *pmax = vmax(points, 0, -2);
 
  double vgap = pmax->val - pmin->val;
  double avg = (pmin->val + pmax->val) / 2;

  double midval = points[points.size()-1]->val;

  if (pflag)
    cout << "vgap: " << vgap << " min: " << pmin->val << " max: " << pmax->val << " avg: " << avg << " midval: " << midval << endl;

  if (fabs(midval - avg) > 10 * vgap) return true;

 vector<double> base, up;
 vector<int> ibase, iup;
  for (int ip = 0; ip < points.size() - 1; ip++)
    if (points[ip]->val < avg) {
      ibase.push_back(ip);
      base.push_back(points[ip]->val);
    } else if (points[ip]->val > avg) {
      iup.push_back(ip);
      up.push_back(points[ip]->val);
    }

  double width_base = width(base),
    mid_base = mid(base),
    width_up = width(up),
    mid_up = mid(up);
  if (pflag)
    cout << "nbase: " << base.size() << " nup: " << up.size() << " base: " << mid_base << " +- " << width_base << " up: " << mid_up << " +- " << width_up << endl;

  if (fmax(width_base, width_up) / vgap > 0.1) return false;
  if (fabs(width_base / mid_base) > 0.05 || fabs(width_up / mid_up) > 0.05) return false;

  if (base.size() == points.size() / 2 && up.size() == points.size() / 2) {
    int label = points.size() - 2; // the last index
    vector<int> *ind = &ibase;
    if (*iup.rbegin() == label) ind = &iup;

    for (auto i : *ind) label &= i;
    bool is_slicing = (label != 0); // the label would be nonzero only when the up points are exactly inside a slice

    if (pflag)
      cout << "is_slicing: " << is_slicing << endl;
    if (is_slicing) return false; 
  }

  return true;
}

KDInterp* hcs_interp(const HCS& hcs, bool pflag) {
  map<vec_t, Vec, vector_less_than> p_cs_tab;
  Vec p_cs;
  auto dist = [&](const vector<double>& x) {
    double r = x[0],
           theta = x[1],
           phi0 = x[2];
    double phi = phi0 - r * hcs.Omega / hcs.Vs_eq;

    double res = hcs.get_distance(r, theta, phi, p_cs);
    p_cs_tab.insert(pair<vec_t, Vec>(x, p_cs));

    static int iter = 0;
    if (pflag && iter++ % 5000 == 0)
      cout << "counting: " << setprecision(16) << r / AU << " " << theta / deg << " " << phi0 / deg << " | " << res / AU << endl;
    return res;
  };

  auto dist_corr = [&](const vec_t& x, Vec& point) -> double {
    double r = x[0],
           theta = x[1],
           phi0 = x[2];
    double phi = phi0 - r * hcs.Omega / hcs.Vs_eq;
    Vec target;
    target.set_spherical(r, theta, phi);
    double res = hcs.sign(r, theta, phi) * hcs.get_distance_from_point(target, point);

    static int iter = 0;
    if (pflag && iter++ % 5000 == 0)
      cout << "counting: " << setprecision(16)
        << r / AU << " " << theta / deg << " " << phi0 / deg << " from " << point.len() / AU << " " << point.theta() / deg << " " << point.phi() / deg
        << " | " << res / AU << endl;
    return res;
  };

  auto tab_corr = [&](const vector<KDPoint*>& points) -> set<KDPoint*> {
    set<KDPoint*> res;

    if (!distance_jump(points) && !step_shape(points)) return res;

    KDPoint *pmin = absmin(points);
    KDPoint *pmax = absmax(points);
    double avg = (pmin->val + pmax->val) / 2;
    for (int i = 0; i < points.size(); i++) {
      if (fabs(points[i]->val) < avg) continue;

      p_cs = p_cs_tab.find(pmin->x)->second;
      double vcorr = dist_corr(points[i]->x, p_cs);

      if (fabs(vcorr) < fabs(points[i]->val)) {
        points[i]->val = vcorr;
        p_cs_tab.find(points[i]->x)->second = p_cs;
        res.insert(points[i]);
      }
    }
   return res;
  };

  return new KDInterp(dist,
                      {0.05 * AU, pi / 2 - HCS::angle - 5 * deg, 0},
                      {120 * AU,  pi / 2 + HCS::angle + 5 * deg, 2 * pi + 1e-5 * deg},
                       1e-4 * AU, 1, tab_corr);
}

double hcs_interp_eval(double r, double theta, double phi, KDInterp *intp, const HCS& hcs) {
  double phi0 = phi + r * hcs.Omega / hcs.Vs_eq;
  phi0 = fmod(phi0, 2 * pi);
  vector<double> x = {r, theta, phi0};
  return (*intp)(x);
}