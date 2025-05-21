#ifndef HCS_INTERP_H
#define HCS_INTERP_H
#include "KDInterpSide.h"

class HCS;
KDInterpSide* hcs_interp(const HCS& hcs, bool pflag = false);
double hcs_interp_eval(double r, double theta, double phi, KDInterpSide *intp, const HCS& hcs, bool pflag = false);
KDInterpSide* hcs_interp(bool pflag = false);
double hcs_interp_eval(double angle, double r, double theta, double phi, KDInterpSide *intp, const HCS& hcs);
#endif /* HCS_INTERP_H */
