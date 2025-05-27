#ifndef HCS_INTERP_H
#define HCS_INTERP_H
#include "KDInterp.h"

class HCS;
KDInterp* hcs_interp(const HCS& hcs, bool pflag = false);
double hcs_interp_eval(double r, double theta, double phi, KDInterp *intp, const HCS& hcs, bool pflag = false);
KDInterp* hcs_interp(double angle_low, double angle_up, double resolution, int ix_split = 0, bool pflag = false);
double hcs_interp_eval(double angle, double r, double theta, double phi, KDInterp *intp, const HCS& hcs, bool pflag = false);

void print_block_d4(const KDValueSide* kd);
#endif /* HCS_INTERP_H */
