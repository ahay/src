#ifndef __XWAVE_ASG_READ__
#define __XWAVE_ASG_READ__

#include "defaults.h"
#include "sgn_indices.h"
#include "model.h"
#include "gridio.h"

int asg_readgrid(PARARRAY *, FILE *, IMODEL *);
int asg_readtimegrid(PARARRAY *, FILE *, IMODEL *);
int asg_readmedia(PARARRAY *, FILE *, IMODEL *, int);

int init_acoustic_geom_par(grid * _g, PARARRAY par, FILE * fp);
#endif
