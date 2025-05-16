#ifndef HCS_INTERP_H
#define HCS_INTERP_H
#include "KDInterp.h"

class HCS;
KDInterp* hcs_interp(const HCS& hcs, bool pflag = false);
double hcs_interp_eval(double r, double theta, double phi, KDInterp *intp, const HCS& hcs);
#endif /* HCS_INTERP_H */
