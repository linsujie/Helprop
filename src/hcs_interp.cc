#include <map>
#include <vector>
#include <iomanip>
#include "particle.h"
#include "HCS.h"
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

bool distance_jump_d4(const std::vector<KDPoint*>& points, bool pflag = false) {
  auto pmin = absmin(points),
   pmax = absmax(points);
  double vgap = fabs(pmax->val - pmin->val);

  auto xmin = pmin->real_x();
  auto xmax = pmax->real_x();

  Vec p1, p2;
  p1.set_spherical(xmin[1], pi / 2 + xmin[2] * xmin[0], xmin[3]);
  p2.set_spherical(xmax[1], pi / 2 + xmax[2] * xmax[0], xmax[3]);

  return (p1 - p2).len() < vgap;
}

bool distance_jump(const std::vector<KDPoint*>& points, bool pflag = false) {
  auto pmin = absmin(points),
   pmax = absmax(points);
  double vgap = fabs(pmax->val - pmin->val);

  auto xmin = pmin->real_x();
  auto xmax = pmax->real_x();

  Vec p1, p2;
  p1.set_spherical(xmin[0], xmin[1], xmin[2]);
  p2.set_spherical(xmax[0], xmax[1], xmax[2]);

  if (pflag)
    cout << "vgap: " << vgap / AU << " distance: " << (p1 - p2).len() / AU << endl;
  return (p1 - p2).len() < vgap;
}

bool distance_jump_side(const std::vector<KDPoint*>& points, bool pflag = false) {
  auto pmid = points.back();

  for (int ix = 0; ix < pmid->x.size(); ix++)
    if (distance_jump({points[2*ix], pmid, points[2*ix+1]}, pflag)) return true;

  return false;
}

bool step_val(double v1, double v2, double v3, bool pflag = false) {
  double l1 = v2 - v1,
         l2 = v3 - v2;
  double w = l1 * l2 > 0 ? fabs(l1 + l2) : fmax(fabs(l1), fabs(l2));
  if (pflag)
    cout << "step_val: " << v1 / AU << " " << v2 / AU << " " << v3 / AU << endl;
  return fmin(fabs(l1 + l2), fmin(fabs(l1), fabs(l2))) / w < 0.05;
}

bool step_shape_p3(const std::vector<KDPoint*>& points, bool pflag = false) {
  if (points.size() < 3) return false;
  assert(points.size() == 3);
  if (!sameside(points, pflag)) return false;

  return step_val(points[0]->val, points[1]->val, points[2]->val, pflag);
}

bool step_shape_side(const std::vector<KDPoint*>& points, bool pflag = false) {
  if (!sameside(points, pflag)) return false;

  auto pmid = points.back();
  for (int ix = 0; ix < pmid->x.size(); ix++)
    if (step_val(points[2*ix]->val, pmid->val, points[2*ix+1]->val, pflag)) return true;

  return false;
}

bool within(double a, double l, double u) { return l < a && a < u; }
bool step_shape(const std::vector<KDPoint*>& points, bool pflag = false) {
  if (!sameside(points, pflag)) return false;

  KDPoint* pmin = vmin(points, 0, -2),
   *pmax = vmax(points, 0, -2);
 
  double vgap = pmax->val - pmin->val;
  double avg = (pmin->val + pmax->val) / 2;

  double midval = points.back()->val;

  if (pflag)
    cout << "vgap: " << vgap / AU << " min: " << pmin->val / AU << " max: " << pmax->val / AU << " avg: " << avg / AU << " midval: " << midval / AU << endl;

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
    cout << "nbase: " << base.size() << " nup: " << up.size() << " base: " << mid_base / AU << " +- " << width_base / AU << " up: " << mid_up / AU << " +- " << width_up / AU << endl;

  if (fmax(width_base, width_up) / vgap > 0.1) return false;
  if (fabs(width_base / mid_base) > 0.05 || fabs(width_up / mid_up) > 0.05) return false;

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

  auto dist_corr = [&](const vec_t& x, Vec& point, double ref_val) -> double {
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
        << " | " << res / AU << "->" << ref_val/AU << " " << (fabs(res) < fabs(ref_val)) << endl;
    return res;
  };

  auto tab_corr = [&](const vector<KDPoint*>& points) -> set<KDPoint*> {
    set<KDPoint*> res;

    if (!distance_jump(points) && !step_shape_p3(points))
      return res;

    KDPoint *pmin = absmin(points);
    KDPoint *pmax = absmax(points);
    double avg = (pmin->val + pmax->val) / 2;
    for (int i = 0; i < points.size(); i++) {
      if (fabs(points[i]->val) < avg) continue;

      auto ip_cs = p_cs_tab.find(pmin->real_x());
      assert(ip_cs != p_cs_tab.end());
      p_cs = ip_cs->second;
      double vcorr = dist_corr(points[i]->real_x(), p_cs, points[i]->val);

      if (fabs(vcorr) < fabs(points[i]->val)) {
        points[i]->val = vcorr;
        p_cs_tab.find(points[i]->real_x())->second = p_cs;
        res.insert(points[i]);
      }
    }
    return res;
  };

  return new KDInterp(dist,
                      {60.05 * AU, pi / 2, 180.0001 * deg},
                      {60 * AU,  HCS::angle + 5 * deg, 180.0001 * deg},
                       { 5e-4 * AU, 5e-4 * AU, 5e-4 * AU }, { 2, 2, 2 }, 0, tab_corr);
}

double sum(const vec_t& vec) {
  double res = 0;
  for (auto v : vec) res += v;
  return res;
}
double hcs_interp_eval(double r, double theta, double phi, KDInterp *intp, const HCS& hcs, bool pflag) {
  double phi0 = phi + r * hcs.Omega / hcs.Vs_eq;
  phi0 = fmod(phi0, 2 * pi);
  vector<double> x = {r, theta, phi0};
  return (*intp)(x);
}

bool approx(double a, double b) {
  return fabs(a - b) / fabs(a + b) < 1e-6;
}

void get_r_theta_phi(const vector<double>& x, double& r, double& theta, double& phi, const HCS& hcs) {
  HCS::angle = x[0];
  r = x[1];
  theta = pi / 2 + x[2] * HCS::angle;
  double phi0 = x[3];
  phi = phi0 - r * hcs.Omega / hcs.Vs_eq;
}
void print_block_d4(const KDValueSide* kd) {
  particle p;
  HCS::angle = 15 * deg;
  HCS::hcsform = HCS::Kota_Jokipii;
  HCS hcs(p.Wind(), "");
  hcs.resolution = 1e-7 * AU;

  auto print_point = [&](KDPoint* p) {
    double r, theta, phi;
    get_r_theta_phi(p->real_x(), r, theta, phi, hcs);

    cout << setprecision(16) << p->x[0] << " " << p->x[1] << " " << p->x[2] << " " << p->x[3] << " | " << HCS::angle / deg << " " << r / AU << " " << theta / deg << " " << phi / deg << " | " << p->val / AU <<  endl;
  };

  cout << "--------------sides-----------------" << endl;
  for (auto p : kd->sides) print_point(p);
  cout << "-------------corners----------------" << endl;
  for (auto p : kd->corners) print_point(p);
  cout << "---------------mid------------------" << endl;
  print_point(kd->pmid);
  cout << "errs: " << kd->err[0] / AU << " " << kd->err[1] / AU << " " << kd->err[2] / AU << " " << kd->err[3] / AU << " | " << kd->errmax / AU << endl;
  cout << "wids: " << kd->width[0] << " " << kd->width[1] << " " << kd->width[2] << " " << kd->width[3] << endl;

  auto tol_local = kd->interp->norm_tol;
  for (int ix = 0; ix < kd->interp->index.dim; ix++)
    tol_local[ix] /= kd->width[ix]; // smaller grid are allowed to have larger tolerance
  cout << "tols: " << tol_local[0] / AU << " " << tol_local[1] / AU << " " << tol_local[2] / AU << " " << tol_local[3] / AU << endl;

  cout << "ix_split: " << kd->ix_split << endl;
}

KDInterp* hcs_interp(double angle_low, double angle_up, double resolution, int ix_split, bool pflag) {
  particle p;
  HCS::angle = 15 * deg;
  HCS::hcsform = HCS::Kota_Jokipii;
  HCS hcs(p.Wind(), "");
  hcs.resolution = 1e-7 * AU;

  map<vec_t, Vec, vector_less_than> p_cs_tab;
  Vec p_cs;

  auto dist = [&](const vector<double>& x) -> double {
    double r, theta, phi;
    get_r_theta_phi(x, r, theta, phi, hcs);
    double res = hcs.get_distance(r, theta, phi, p_cs);

    static int iter = 0;
    if (pflag && iter++ % 5000 == 0) {
      Vec target;
      target.set_spherical(r, theta, phi);
      cout << "counting: " << iter << " " << setprecision(16) << HCS::angle / deg << " " << r / AU << " " << theta / deg << " " << phi / deg << " | " << res / AU << " " << (target - p_cs).len() / AU << endl;
    }

    assert(p_cs_tab.find(x) == p_cs_tab.end());
    p_cs_tab.insert(pair<vec_t, Vec>(x, p_cs));

    return res;
  };

  auto dist_corr = [&](KDPoint* p, Vec& point, int iter) -> double {
    const auto x = p->real_x();
    double r, theta, phi;
    get_r_theta_phi(x, r, theta, phi, hcs);

    Vec target;
    target.set_spherical(r, theta, phi);
    double res = hcs.sign(r, theta, phi) * hcs.get_distance_from_point(target, point);

    static int ncorr = 0;
    ncorr++;
    if (pflag && iter++ % 5000 == 0)
      cout << "counting: " << ncorr << " " << setprecision(16) << HCS::angle /deg << " "
        << r / AU << " " << theta / deg << " " << phi / deg << " from " << point.len() / AU << " " << point.theta() / deg << " " << point.phi() / deg
        << " | " << res / AU << " " << p->val / AU << " " << (fabs(res) < fabs(p->val)) << " " << (target - point).len() / AU << endl;
    return res;
  };

  auto tab_corr = [&](const vector<KDPoint*>& points) -> set<KDPoint*> {
    set<KDPoint*> res;

    bool dj = distance_jump_d4(points);
    bool ss = step_shape_p3(points);
    if (!dj && !ss) return res;

    static int iter = 0;
    iter++;
    KDPoint *pmin = absmin(points);
    KDPoint *pmax = absmax(points);
    double avg = fabs(pmin->val + pmax->val) / 2;
    if (pflag && iter % 5000 == 0) {
      cout << "-----------correction-------------" << endl;
      for (auto& p : points)
        cout << p->x[0] << " " << p->x[1] << " " << p->x[2] << " " << p->x[3] << " | " << p->val / AU << endl;
      for (auto& p : points) {
        double r, theta, phi;
        get_r_theta_phi(p->real_x(), r, theta, phi, hcs);
        Vec target;
        target.set_spherical(r, theta, phi);
        auto pcs = p_cs_tab.at(p->real_x());

        cout << setprecision(16) << HCS::angle / deg << " " << r / AU << " " << theta / deg << " " << phi / deg << " | " << p->val / AU << " | " << (target - pcs).len() / AU << endl;
      }

      cout << "min: " << pmin->real_x()[0] / deg << " " << pmin->real_x()[1] / AU << " " << (pi / 2 + pmin->real_x()[0] * pmin->real_x()[2]) / deg << " " << pmin->real_x()[3] / deg << " | " << pmin->val / AU << endl;
      cout << "dj: " << dj << " ss: " << ss << endl;
    }
    
    for (int i = 0; i < points.size(); i++) {
      if (fabs(points[i]->val) < avg) continue;

      auto ip = p_cs_tab.find(pmin->real_x());
      assert(ip != p_cs_tab.end());
      p_cs = ip->second;
      double vcorr = dist_corr(points[i], p_cs, iter);

      if (fabs(points[i]->val) - fabs(vcorr) > 1e-5 * AU) {
        points[i]->val = vcorr;
        auto ipcs = p_cs_tab.find(points[i]->real_x());
        assert(ipcs != p_cs_tab.end());
        ipcs->second = p_cs;
        res.insert(points[i]);
      }
    }
    return res;
  };

  double angle_mid = (angle_low + angle_up) / 2;
  double angle_width = (angle_up - angle_low) / 2;

  double tol0 = resolution * AU;
  vec_t ref_width = { 20 * deg, 50 * AU, 2,  2 * pi };
  vec_t tol;
  // The tolerance is inversely proportional to the width; consequently, the reference tolerance parameter is defined as tol × width.
  for (auto& w : ref_width) tol.push_back(tol0 * w);

  return new KDInterp(dist,
                      {angle_mid * deg, 60.05 * AU, 0, 180.00001 * deg },
                      {angle_width * deg, 60 * AU,  1.1, 180.00001 * deg},
                      tol,
                      {0, 0, 0, 1}, ix_split, tab_corr);
}

double hcs_interp_eval(double angle, double r, double theta, double phi, KDInterp *intp, const HCS& hcs, bool pflag) {
  double phi0 = phi + r * hcs.Omega / hcs.Vs_eq;
  phi0 = fmod(phi0, 2 * pi);

  double theta_rel = (theta - pi / 2) / angle;
  vector<double> x = {angle, r, theta_rel, phi0};
 
  return (*intp)(x);
}