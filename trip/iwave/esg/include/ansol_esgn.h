/*
  Analytical solution version for esgn.h
*/

#ifndef __ANSOL_ESGN_H_
#define __ANSOL_ESGN_H_


#include "fd.h"
#include "esgn_read.h"
#include "create_sten_2k.h"
#include "model_action.h"
#include "defaults.h"
#include "traceio.h"
#include "pointsrc.h"

// #include "ansol_esgsteps.h"


/* Source info struct for analytical solution.
   Values set in ansol_HI_esg_srcinfo. */
typedef struct {
	float amp;	
	float fpeak;
	RPNT  src_d; //source direction
} RICKER_INFO;

/* Medium info struct for analytical solution.
   Values set in ansol_HI_esg_medinfo. */
typedef struct {
	float mu;
	float lambda;
	float rho;   //density
	float alpha; //p-velocity
	float beta;  //s-velocity
} ESG_HIMED_INFO;

//parameter struct
typedef struct {
	//old parameters
	ireal dt;	//set in ansol_esg_readschemeinfo
  	int   k;	//set in ansol_HI_esg_modelinit
  	int   ndim;	//set in ansol_esg_readschemeinfo
  	RPNT  lam;

	//new parameters
	IMODEL * link_mdl;	//link to fd MODEL needed in time step function computations (i.e., physical dom and time t)
	RPNT  dx; 		//step sizes, set in ansol_esg_readschemeinfo
	RPNT  o_coord;		//origin coordinates, set in ansol_esg_readschemeinfo
	IPNT  o_index;		//origin indexes, assumed in primal grid, set in ansol_esg_readschemeinfo
	RPNT  src_x;		//coordinates of source term, set in ansol_HI_esg_srcinfo
	ireal t_off;		//total time offset

	IPNT gtype[RDOM_MAX_NARR];

	void * srcinfo;	//source information
	void * medinfo;	//medium information
} ANSOL_ESG_PARS;


int ansol_HI_esg_modelinit( PARARRAY *pars, FILE *stream, IMODEL *model );
int ansol_HI_esg_step( RDOM*, int, void* );


//post model construction functions that set srcinfo and medinfo
//NOTE: these functions must be called after initializing IMODEL structs for both fd and sol.
int ansol_HI_esg_medinfo( FILE * stream, IMODEL * model );
int ansol_HI_esg_srcinfo( FILE * stream, PARARRAY * pars, IMODEL * model, tracegeom * tg, POINTSRC * pntsrc );
int ansol_link_mdl( FILE * stream, IMODEL * sol_mdl, IMODEL * fd_mdl );
void ansol_esg_fprintf( FILE * stream, IMODEL * sol_mdl );



#endif /*__ANSOL_ESGN_H_*/
