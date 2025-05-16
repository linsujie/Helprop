#include "KDInterp.h"
#include "HCS.h"

KDInterp* hcs_interp(const HCS& hcs, bool pflag = false);
double hcs_interp_eval(double r, double theta, double phi, KDInterp *intp, const HCS& hcs);