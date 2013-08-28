#include "esgn.h"

#ifdef __ESG_STEPS__

#include "rdomain.h"

int esgn_gts2d_24 ( RDOM *, RDOM *, RDOM *, int, void * );
// int esgn_gts3d_24 ( RDOM *, RDOM *, RDOM *, int, void * );
// int esgn_gts2d_210( RDOM *, RDOM *, RDOM *, int, void * );
// int esgn_gts3d_210( RDOM *, RDOM *, RDOM *, int, void * );

int esg_step_s( RDOM * dom, void * pars );
int esg_step_v( RDOM * dom, void * pars );

#endif
