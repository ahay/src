#ifndef __ESGN_READ__
#define __ESGN_READ__

#include "defaults.h"
#include "esgn_indices.h"
#include "model.h"
#include "gridio.h"

int esg_readgrid(PARARRAY *, FILE *, IMODEL *);
int esg_readtimegrid(PARARRAY *, FILE *, IMODEL *);
int esgn_readmedia(PARARRAY *, FILE *, IMODEL *, int);
//int esg_readmedia(PARARRAY *, FILE *, IMODEL *, int);

int init_elastic_geom_par(grid * _g, PARARRAY par, FILE * fp);
#endif
