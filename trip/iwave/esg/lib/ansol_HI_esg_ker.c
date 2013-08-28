//Author: Mario Bencomo (summer 2013)


/* Analytical solution kernels for the 2D and 3D linear elaticity equations for isotropic homogeneous medium with point source f:
	f(x,t) = delta(xs-x)psi(t)
where psi(t) is the Ricker wavelet.

For the analytical formulas from which this kernel is based on, I refer to:
	Aki and Richards, "Quantitative Seismology: Theory and Methods", volume 1, 2002.*/

#include <math.h>
#include <stdio.h>

#include "ansol_esgn.h"
#include "ansol_esgsteps.h"
#include "esgn_indices.h"
#include "gauss.h"

#define SQR(x) ((x)*(x))
#define WARN 131313

#ifndef TOL
#define TOL 1.0e-8
#endif 

// #define MARIO_VERBOSE

//2D FUNCTIONS
int ansol_HI_esg_2dp0 ( RDOM * dom, void * pars ); //updates pz
int ansol_HI_esg_2dp1 ( RDOM * dom, void * pars ); //updates px
int ansol_HI_esg_2dv0 ( RDOM * dom, void * pars ); //updates vz
int ansol_HI_esg_2dv1 ( RDOM * dom, void * pars ); //updates vx
int ansol_HI_esg_2dss0( RDOM * dom, void * pars ); //updates szx

//3D FUNCTIONS
int ansol_HI_esg_3dp0 ( RDOM * dom, void * pars ); //updates pz
int ansol_HI_esg_3dp1 ( RDOM * dom, void * pars ); //updates px
int ansol_HI_esg_3dp2 ( RDOM * dom, void * pars ); //updates py
int ansol_HI_esg_3dv0 ( RDOM * dom, void * pars ); //updates vz
int ansol_HI_esg_3dv1 ( RDOM * dom, void * pars ); //updates vx
int ansol_HI_esg_3dv2 ( RDOM * dom, void * pars ); //updates vy
int ansol_HI_esg_3dss0( RDOM * dom, void * pars ); //updates szx
int ansol_HI_esg_3dss1( RDOM * dom, void * pars ); //updates sxy
int ansol_HI_esg_3dss2( RDOM * dom, void * pars ); //updates syz


//some auxilary function declarations
//============================================================================//

//computes integral_(r/alpha)^(r/beta) Ricker(t-tau) dtau
float Int1( float t, 
            float r, 
            float alpha, 
            float beta, 
            float fpeak );

//computes integral_(r/alpha)^(r/beta) tau Ricker(t-tau) dtau
float Int2( float t, 
	    float r, 
	    float alpha, 
	    float beta, 
	    float fpeak );

//computes integral_(r/alpha)^(r/beta) tau Ricker_prime(t-tau) dtau
float Int3( float t, 
	    float r, 
	    float alpha, 
	    float beta, 
	    float fpeak );

//computes integral_(r/alpha)^(r/beta) tau^2 Ricker_prime(t-tau) dtau
float Int4( float t, 
	    float r, 
	    float alpha, 
	    float beta, 
	    float fpeak );


//computes derivative of Ricker
float Ricker_prime( float t,
	      	    float fpeak  );

//Kronocker delta, delta(i,j) = 1. if i=j, else 0.
float Kdelta( int i, int j );

//computes the kth derivative of the ith component of the displacement field u, do to a point source in the jth direction.
float uij_primexk( int   * indexes, 
		   int     dim, 
		   float * x, 
		   float * xs,
		   float   t,
		   float   alpha, 
		   float   beta, 
		   float   fpeak,
		   FILE  * stream   );

//computes the time derivative of the ith component of the displacement field u, do to a point source in the jth direction.
float uij_primet( int   * indexes, 
		  int     dim, 
		  float * x, 
		  float * xs,
		  float   t,
		  float   alpha, 
		  float   beta, 
		  float   fpeak,
         	  FILE  * stream  );

//============================================================================//


//#############################################################################//
//########################### 2D CODE #########################################//
//#############################################################################//

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_ker2d( RDOM * dom, int iarr, void * pars ){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,">> Inside ansol_HI_esg_ker2d\n");
#endif
 	if ( iarr == D_P0 ) return ansol_HI_esg_2dp0 ( dom, pars );
 	if ( iarr == D_P1 ) return ansol_HI_esg_2dp1 ( dom, pars );
 	if ( iarr == D_V0 ) return ansol_HI_esg_2dv0 ( dom, pars );
 	if ( iarr == D_V1 ) return ansol_HI_esg_2dv1 ( dom, pars );
 	if ( iarr == D_S0 ) return ansol_HI_esg_2dss0( dom, pars );

  	return E_NOTIMESTEP;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_2dp0( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,"> Inside ansol_HI_esg_2dp0\n");
#endif

	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ijk_0, ijk_1;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primexk;
	
	float r2;
	RPNT x_temp;
	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	int iarr = D_P0; 	

	//field to be updated
	register ireal **_p0 = s[iarr]._s2;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting parameters from medium struct */
	register ireal multi0 = himedinfo->lambda + 2.0*himedinfo->mu;
	register ireal multi1 = himedinfo->lambda;
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;
 
	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp     = rickerinfo->amp;
	register ireal fpeak   = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt;
// fprintf(stderr,"time = %e, and ricker(time) = %e\n",t,comprick(t,fpeak));

	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 2D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

// fprintf(stderr,"offset = (%g,%g)\n",offset0,offset1);

	/* loop in 1-direction */
	for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
		/* updating 1-component of x */
		x[1] = offset[1] + index[1]*dx[1];

		/* loop in 0-direction */
		for ( index[0]=0; index[0]<nx[0]; index[0]++) {
			/* updating 0-component of x */
			x[0] = offset[0] + index[0]*dx[0];

			//checking for case where x = src_x
			r2 = 0;
			for (idim=0; idim<ndim; idim++)
				r2 += SQR(x[idim]-src_x[idim]);
			//case where x = src_x
			if (fabs(r2)<1.0e-9){
#ifdef MARIO_VERBOSE
				fprintf(stderr,"WARNING: x = srx_x for [i0,i1] = [%d,%d]\n",index[0],index[1]);
				fprintf(stderr,"x = (%g,%g)\n",x[0],x[1]);
				fprintf(stderr,"src_x = (%g,%g)\n",src_x[0],src_x[1]);	
#endif
				flag++;
				if (flag>1){
				//case where flag was already set off
					fprintf(stderr,"ERROR: x=src_x in more than one point! ABORT\n");
					exit(1);
				}

				_p0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] = 0;

				//computing analytical solution at that src_x using average of neighbors
				for (idim=0; idim<ndim; idim++){
					// computing field at x+dx
					RASN(x_temp,x);
					x_temp[idim] = x[idim]+dx[idim];
 
					/* updating ijk indexes */
					ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 0; //(i,j,k)=(0,0,0)
					ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 0; //(i,j,k)=(0,1,0)
					/* updating field */
					u_primexk[0] = uij_primexk( ijk_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primexk[1] = uij_primexk( ijk_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_p0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi0*( src_d[0]*u_primexk[0] 
                                                                                  + src_d[1]*u_primexk[1] );
					/* updating ijk indexes */
					ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 1; //(i,j,k)=(1,0,1)
					ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 1; //(i,j,k)=(1,1,1)
					/* updating field */
					u_primexk[0] = uij_primexk( ijk_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primexk[1] = uij_primexk( ijk_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					
					_p0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi1*( src_d[0]*u_primexk[0] 
										  + src_d[1]*u_primexk[1] );
					// computing field at x-dx
					RASN(x_temp,x);
					x_temp[idim] = x[idim]-dx[idim];

					/* updating ijk indexes */
					ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 0; //(i,j,k)=(0,0,0)
					ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 0; //(i,j,k)=(0,1,0)
					/* updating field */
					u_primexk[0] = uij_primexk( ijk_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primexk[1] = uij_primexk( ijk_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_p0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi0*( src_d[0]*u_primexk[0] 
                                                                                  + src_d[1]*u_primexk[1] );
					/* updating ijk indexes */
					ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 1; //(i,j,k)=(1,0,1)
					ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 1; //(i,j,k)=(1,1,1)
					/* updating field */
					u_primexk[0] = uij_primexk( ijk_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primexk[1] = uij_primexk( ijk_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					
					_p0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi1*( src_d[0]*u_primexk[0] 
										  + src_d[1]*u_primexk[1] );
				}
				_p0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] /= 4.0;


			}
			//case where x!=src_x
			else {
				/* updating ijk indexes */
				ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 0; //(i,j,k)=(0,0,0)
				ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 0; //(i,j,k)=(0,1,0)
	
				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				_p0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] = multi0*( src_d[0]*u_primexk[0] 
									  + src_d[1]*u_primexk[1] );
	
				/* updating ijk indexes */
				ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 1; //(i,j,k)=(1,0,1)
				ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 1; //(i,j,k)=(1,1,1)
	
				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				
				_p0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi1*( src_d[0]*u_primexk[0] 
									  + src_d[1]*u_primexk[1] );
	
			}
			_p0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;

		}
	}
	return 0;
}


/*----------------------------------------------------------------------------*/
int ansol_HI_esg_2dp1( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,"> Inside ansol_HI_esg_2dp1\n");
#endif

	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ijk_0, ijk_1;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primexk;

	float r2;
	RPNT x_temp;
	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	register int   iarr = D_P1;

	//field to be updated
	register ireal **_p1 = s[iarr]._s2;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting parameters from medium struct */
	register ireal multi0 = himedinfo->lambda + 2.0*himedinfo->mu;
	register ireal multi1 = himedinfo->lambda;
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;

	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp     = rickerinfo->amp;
	register ireal fpeak = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );	
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt;

	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 2D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

// fprintf(stderr,"offset = (%g,%g)\n",offset0,offset1);

	/* loop in 1-direction */
	for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
		/* updating 1-component of x */
		x[1] = offset[1] + index[1]*dx[1];

		/* loop in 0-direction */
		for ( index[0]=0; index[0]<nx[0]; index[0]++) {
			/* updating 0-component of x */
			x[0] = offset[0] + index[0]*dx[0];

			//checking for case where x = src_x
			r2 = 0;
			for (idim=0; idim<ndim; idim++)
				r2 += SQR(x[idim]-src_x[idim]);
			//case where x = src_x
			if (fabs(r2)<1.0e-9){
#ifdef MARIO_VERBOSE
				fprintf(stderr,"WARNING: x = srx_x for [i0,i1] = [%d,%d]\n",index[0],index[1]);
				fprintf(stderr,"x = (%g,%g)\n",x[0],x[1]);
				fprintf(stderr,"src_x = (%g,%g)\n",src_x[0],src_x[1]);	
#endif
				flag++;
				if (flag>1){
				//case where flag was already set off
					fprintf(stderr,"ERROR: x=src_x in more than one point! ABORT\n");
					exit(1);
				}

				_p1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] = 0;
				//computing analytical solution at that src_x using average of neighbors
				for (idim=0; idim<ndim; idim++){
					// computing field at x+dx
					RASN(x_temp,x);
					x_temp[idim] = x[idim]+dx[idim];
 
					/* updating ijk indexes */
					ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 1; //(i,j,k)=(1,0,1)
					ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 1; //(i,j,k)=(1,1,1)
		
					/* updating field */
					u_primexk[0] = uij_primexk( ijk_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primexk[1] = uij_primexk( ijk_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_p1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi0*( src_d[0]*u_primexk[0] 
										  + src_d[1]*u_primexk[1] );
		
					/* updating ijk indexes */
					ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 0; //(i,j,k)=(0,0,0)
					ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 0; //(i,j,k)=(0,1,0)
		
					/* updating field */
					u_primexk[0] = uij_primexk( ijk_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primexk[1] = uij_primexk( ijk_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_p1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi1*( src_d[0]*u_primexk[0] 
										  + src_d[1]*u_primexk[1] );
					// computing field at x-dx
					RASN(x_temp,x);
					x_temp[idim] = x[idim]-dx[idim];

					/* updating ijk indexes */
					ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 1; //(i,j,k)=(1,0,1)
					ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 1; //(i,j,k)=(1,1,1)
		
					/* updating field */
					u_primexk[0] = uij_primexk( ijk_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primexk[1] = uij_primexk( ijk_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_p1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi0*( src_d[0]*u_primexk[0] 
										  + src_d[1]*u_primexk[1] );
		
					/* updating ijk indexes */
					ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 0; //(i,j,k)=(0,0,0)
					ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 0; //(i,j,k)=(0,1,0)
		
					/* updating field */
					u_primexk[0] = uij_primexk( ijk_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primexk[1] = uij_primexk( ijk_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_p1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi1*( src_d[0]*u_primexk[0] 
										  + src_d[1]*u_primexk[1] );
				}
				_p1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] /= 4.0;


			}
			//case where x!=src_x
			else {
				/* updating ijk indexes */
				ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 1; //(i,j,k)=(1,0,1)
				ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 1; //(i,j,k)=(1,1,1)
	
				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				_p1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] = multi0*( src_d[0]*u_primexk[0] 
									  + src_d[1]*u_primexk[1] );
	
				/* updating ijk indexes */
				ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 0; //(i,j,k)=(0,0,0)
				ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 0; //(i,j,k)=(0,1,0)
	
				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				_p1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi1*( src_d[0]*u_primexk[0] 
									  + src_d[1]*u_primexk[1] );	
			}
			_p1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;
		}
	}
	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_2dv0( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,"> Inside ansol_HI_esg_2dv0\n");
#endif

	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ij_0, ij_1;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primet;

	float r2;
	RPNT x_temp;
	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	register int   iarr = D_V0;

	//field to be updated
	register ireal **_v0 = s[iarr]._s2;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting parameters from medium struct */
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;

	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp   = rickerinfo->amp;
	register ireal fpeak = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt; //compensating for substep time offset

	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

#ifdef MARIO_VERBOSE
	fprintf(stderr,"    time = %e\n",get_time( fd_mdl->tsind )-.5*ansolpars->dt);
	fprintf(stderr,"time+off = %e\n",t);
	fprintf(stderr,"     src = %e\n",comprick(t,fpeak));
	fprintf(stderr," src*amp = %e\n",amp*comprick(t,fpeak));
#endif

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 2D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

// fprintf(stderr,"offset = (%g,%g)\n",offset0,offset1);

	/* loop in 1-direction */
	for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
		/* updating 1-component of x */
		x[1] = offset[1] + index[1]*dx[1];

		/* loop in 0-direction */
		for ( index[0]=0; index[0]<nx[0]; index[0]++) {
			/* updating 0-component of x */
			x[0] = offset[0] + index[0]*dx[0];

			//checking for case where x = src_x
			r2 = 0;
			for (idim=0; idim<ndim; idim++)
				r2 += SQR(x[idim]-src_x[idim]);
			//case where x = src_x
			if (fabs(r2)<1.0e-9){
#ifdef MARIO_VERBOSE
				fprintf(stderr,"WARNING: x = srx_x for [i0,i1] = [%d,%d]\n",index[0],index[1]);
				fprintf(stderr,"x = (%g,%g)\n",x[0],x[1]);
				fprintf(stderr,"src_x = (%g,%g)\n",src_x[0],src_x[1]);	
#endif
				flag++;
				if (flag>1){
				//case where flag was already set off
					fprintf(stderr,"ERROR: x=src_x in more than one point! ABORT\n");
					exit(1);
				}

				_v0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] = 0;
				//computing analytical solution at that src_x using average of neighbors
				for (idim=0; idim<ndim; idim++){
					// computing field at x+dx
					RASN(x_temp,x);
					x_temp[idim] = x[idim]+dx[idim];
 
					/* updating ij indexes */
					ij_0[0] = 0; ij_0[1] = 0; //(i,j)=(0,0)
					ij_1[0] = 0; ij_1[1] = 1; //(i,j)=(0,1)
		
					/* updating field */
					u_primet[0] = uij_primet( ij_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primet[1] = uij_primet( ij_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_v0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += src_d[0]*u_primet[0] 
										  + src_d[1]*u_primet[1];

					// computing field at x-dx
					RASN(x_temp,x);
					x_temp[idim] = x[idim]-dx[idim];

					/* updating ij indexes */
					ij_0[0] = 0; ij_0[1] = 0; //(i,j)=(0,0)
					ij_1[0] = 0; ij_1[1] = 1; //(i,j)=(0,1)
		
					/* updating field */
					u_primet[0] = uij_primet( ij_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primet[1] = uij_primet( ij_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_v0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += src_d[0]*u_primet[0] 
										  + src_d[1]*u_primet[1];
				}
				_v0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] /= 4.0;


			}
			//case where x!=src_x
			else {
				/* updating ij indexes */
				ij_0[0] = 0; ij_0[1] = 0; //(i,j)=(0,0)
				ij_1[0] = 0; ij_1[1] = 1; //(i,j)=(0,1)
	
				/* updating field */
				u_primet[0] = uij_primet( ij_0, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				u_primet[1] = uij_primet( ij_1, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				_v0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] = src_d[0]*u_primet[0] 
									  + src_d[1]*u_primet[1];
			}
			_v0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;
		}
	}
	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_2dv1( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,"> Inside ansol_HI_esg_2dv1\n");
#endif

	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ij_0, ij_1;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primet;

	float r2;
	RPNT x_temp;
	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	register int iarr = D_V1; 		//field dependent
	
	//field to be updated
	register ireal ** _v1 = s[iarr]._s2;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting parameters from medium struct */
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;

	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp   = rickerinfo->amp;
	register ireal fpeak = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );	
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt; //compensating for substep time offset

	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 2D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

// fprintf(stderr,"offset = (%g,%g)\n",offset0,offset1);

	/* loop in 1-direction */
	for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
		/* updating 1-component of x */
		x[1] = offset[1] + index[1]*dx[1];

		/* loop in 0-direction */
		for ( index[0]=0; index[0]<nx[0]; index[0]++) {
			/* updating 0-component of x */
			x[0] = offset[0] + index[0]*dx[0];

			//checking for case where x = src_x
			r2 = 0;
			for (idim=0; idim<ndim; idim++)
				r2 += SQR(x[idim]-src_x[idim]);
			//case where x = src_x
			if (fabs(r2)<1.0e-9){
#ifdef MARIO_VERBOSE
				fprintf(stderr,"WARNING: x = srx_x for [i0,i1] = [%d,%d]\n",index[0],index[1]);
				fprintf(stderr,"x = (%g,%g)\n",x[0],x[1]);
				fprintf(stderr,"src_x = (%g,%g)\n",src_x[0],src_x[1]);	
#endif
				flag++;
				if (flag>1){
				//case where flag was already set off
					fprintf(stderr,"ERROR: x=src_x in more than one point! ABORT\n");
					exit(1);
				}

				_v1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] = 0;
				//computing analytical solution at that src_x using average of neighbors
				for (idim=0; idim<ndim; idim++){
					// computing field at x+dx
					RASN(x_temp,x);
					x_temp[idim] = x[idim]+dx[idim];
 
					/* updating ij indexes */
					ij_0[0] = 1; ij_0[1] = 0; //(i,j)=(1,0)
					ij_1[0] = 1; ij_1[1] = 1; //(i,j)=(1,1)
		
					/* updating field */
					u_primet[0] = uij_primet( ij_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primet[1] = uij_primet( ij_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_v1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += src_d[0]*u_primet[0] 
										  + src_d[1]*u_primet[1];

					// computing field at x-dx
					RASN(x_temp,x);
					x_temp[idim] = x[idim]-dx[idim];

					/* updating ij indexes */
					ij_0[0] = 1; ij_0[1] = 0; //(i,j)=(1,0)
					ij_1[0] = 1; ij_1[1] = 1; //(i,j)=(1,1)
		
					/* updating field */
					u_primet[0] = uij_primet( ij_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primet[1] = uij_primet( ij_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_v1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += src_d[0]*u_primet[0] 
										  + src_d[1]*u_primet[1];
				}
				_v1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] /= 4.0;
			}
			//case where x!=src_x
			else {
				/* updating ij indexes */
				ij_0[0] = 1; ij_0[1] = 0; //(i,j)=(1,0)
				ij_1[0] = 1; ij_1[1] = 1; //(i,j)=(1,1)
	
				/* updating field */
				u_primet[0] = uij_primet( ij_0, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				u_primet[1] = uij_primet( ij_1, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				_v1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] = src_d[0]*u_primet[0] 
									  + src_d[1]*u_primet[1];
			}
			_v1 [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;
		}
	}
	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_2dss0( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
#ifdef MARIO_VERBOSE
	fprintf(stderr,"> Inside ansol_HI_esg_2dss0\n");
#endif

	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ijk_0, ijk_1;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primexk;

	float r2;
	RPNT x_temp;
	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	register int   iarr = D_S0;
	
	//field to be updated
	register ireal ** _ss0 = s[iarr]._s2;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting parameters from medium struct */
	register ireal multi     = himedinfo->mu;
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;

	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp   = rickerinfo->amp;
	register ireal fpeak = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt;

	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 2D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

// fprintf(stderr,"offset = (%g,%g)\n",offset0,offset1);

	/* loop in 1-direction */
	for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
		/* updating 1-component of x */
		x[1] = offset[1] + index[1]*dx[1];

		/* loop in 0-direction */
		for ( index[0]=0; index[0]<nx[0]; index[0]++) {
			/* updating 0-component of x */
			x[0] = offset[0] + index[0]*dx[0];
			//checking for case where x = src_x
			r2 = 0;
			for (idim=0; idim<ndim; idim++)
				r2 += SQR(x[idim]-src_x[idim]);
			//case where x = src_x
			if (fabs(r2)<1.0e-9){
#ifdef MARIO_VERBOSE
				fprintf(stderr,"WARNING: x = srx_x for [i0,i1] = [%d,%d]\n",index[0],index[1]);
				fprintf(stderr,"x = (%g,%g)\n",x[0],x[1]);
				fprintf(stderr,"src_x = (%g,%g)\n",src_x[0],src_x[1]);	
#endif
				flag++;
				if (flag>1){
				//case where flag was already set off
					fprintf(stderr,"ERROR: x=src_x in more than one point! ABORT\n");
					exit(1);
				}

				_ss0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] = 0;
				//computing analytical solution at that src_x using average of neighbors
				for (idim=0; idim<ndim; idim++){
					// computing field at x+dx
					RASN(x_temp,x);
					x_temp[idim] = x[idim]+dx[idim];

					/* updating ijk indexes */
					ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 1;//(i,j,k)=(0,0,1)
					ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 1;//(i,j,k)=(0,1,1)
		
					/* updating field */
					u_primexk[0] = uij_primexk( ijk_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primexk[1] = uij_primexk( ijk_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_ss0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi * (src_d[0]*u_primexk[0] 
									 	   + src_d[1]*u_primexk[1]);
					/* updating ijk indexes */
					ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 0;//(i,j,k)=(1,0,0)
					ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 0;//(i,j,k)=(1,1,0)
		
					/* updating field */
					u_primexk[0] = uij_primexk( ijk_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primexk[1] = uij_primexk( ijk_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_ss0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi * (src_d[0]*u_primexk[0] 
										   + src_d[1]*u_primexk[1]);
					// computing field at x-dx
					RASN(x_temp,x);
					x_temp[idim] = x[idim]-dx[idim];

					/* updating ijk indexes */
					ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 1;//(i,j,k)=(0,0,1)
					ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 1;//(i,j,k)=(0,1,1)
		
					/* updating field */
					u_primexk[0] = uij_primexk( ijk_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primexk[1] = uij_primexk( ijk_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_ss0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi * (src_d[0]*u_primexk[0] 
									 	   + src_d[1]*u_primexk[1]);
		
					/* updating ijk indexes */
					ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 0;//(i,j,k)=(1,0,0)
					ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 0;//(i,j,k)=(1,1,0)
		
					/* updating field */
					u_primexk[0] = uij_primexk( ijk_0, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					u_primexk[1] = uij_primexk( ijk_1, 
								ndim, 
								x_temp, 
								src_x,
								t,
								alpha, 
								beta, 
								fpeak,
								stderr   );
					_ss0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi * (src_d[0]*u_primexk[0] 
										   + src_d[1]*u_primexk[1]);
				}
				_ss0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] /= 4.0;
			}
			//case where x!=src_x
			else {
				/* updating ijk indexes */
				ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 1;//(i,j,k)=(0,0,1)
				ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 1;//(i,j,k)=(0,1,1)
	
				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				_ss0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] = multi * (src_d[0]*u_primexk[0] 
									   + src_d[1]*u_primexk[1]);
				/* updating ijk indexes */
				ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 0;//(i,j,k)=(1,0,0)
				ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 0;//(i,j,k)=(1,1,0)
	
				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							ndim, 
							x, 
							src_x,
							t,
							alpha, 
							beta, 
							fpeak,
							stderr   );
				_ss0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi * (src_d[0]*u_primexk[0] 
								  	   + src_d[1]*u_primexk[1]);
			}
			_ss0 [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;
		}
	}
	return 0;
}


//#############################################################################//
//########################### 3D CODE #########################################//
//#############################################################################//

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_ker3d( RDOM * dom, int iarr, void * pars ){
/*----------------------------------------------------------------------------*/
 	if ( iarr == D_P0 ) return ansol_HI_esg_3dp0 ( dom, pars );
 	if ( iarr == D_P1 ) return ansol_HI_esg_3dp1 ( dom, pars );
 	if ( iarr == D_P2 ) return ansol_HI_esg_3dp2 ( dom, pars );
 	if ( iarr == D_V0 ) return ansol_HI_esg_3dv0 ( dom, pars );
 	if ( iarr == D_V1 ) return ansol_HI_esg_3dv1 ( dom, pars );
 	if ( iarr == D_V2 ) return ansol_HI_esg_3dv2 ( dom, pars );
 	if ( iarr == D_S0 ) return ansol_HI_esg_3dss0( dom, pars );
 	if ( iarr == D_S1 ) return ansol_HI_esg_3dss1( dom, pars );
 	if ( iarr == D_S2 ) return ansol_HI_esg_3dss2( dom, pars );

  	return E_NOTIMESTEP;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_3dp0( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ijk_0, ijk_1, ijk_2;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primexk;

	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	register int   iarr = D_P0; 		
	
	//field to be updated
	register ireal *** _p0 = s[iarr]._s3;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting multiples and parameters from medium struct */
	register ireal multi0 = himedinfo->lambda + 2.0*himedinfo->mu; //lambda + 2mu
	register ireal multi1 = himedinfo->lambda; //lambda
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;

	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp     = rickerinfo->amp;
	register ireal fpeak = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );	
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt;

	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 3D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

	/* loop in 2-direction */
	for ( index[2]=0; index[2]<nx[2]; index[2]++ ) {
		/* updating 2-component of x */
		x[2] = offset[2] + index[2]*dx[2];

		/* loop in 1-direction */
		for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
			/* updating 1-component of x */
			x[1] = offset[1] + index[1]*dx[1];

			/* loop in 0-direction */
			for ( index[0]=0; index[0]<nx[0]; index[0]++ ) {
				/* updating 0-component of x */
				x[0] = offset[0] + index[0]*dx[0];

				/* updating ijk indexes */
				ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 0; //(i,j,k)=(0,0,0)
				ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 0; //(i,j,k)=(0,1,0)
				ijk_2[0] = 0; ijk_2[1] = 2; ijk_2[2] = 0; //(i,j,k)=(0,2,0)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				_p0 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] =
									multi0*(  src_d[0]*u_primexk[0]
					                 			+ src_d[1]*u_primexk[1]
					                 	  		+ src_d[2]*u_primexk[2] );
	
				/* updating ijk indexes */
				ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 1; //(i,j,k)=(1,0,1)
				ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 1; //(i,j,k)=(1,1,1)
				ijk_2[0] = 1; ijk_2[1] = 2; ijk_2[2] = 1; //(i,j,k)=(1,2,1)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				_p0 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] +=
							  		multi1*(  src_d[0]*u_primexk[0] 
					                   	 		+ src_d[1]*u_primexk[1]
					                   	  		+ src_d[2]*u_primexk[2] );

				/* updating ijk indexes */
				ijk_0[0] = 2; ijk_0[1] = 0; ijk_0[2] = 2; //(i,j,k)=(2,0,2)
				ijk_1[0] = 2; ijk_1[1] = 1; ijk_1[2] = 2; //(i,j,k)=(2,1,2)
				ijk_2[0] = 2; ijk_2[1] = 2; ijk_2[2] = 2; //(i,j,k)=(2,2,2)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				_p0 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] +=  
									multi1*(  src_d[0]*u_primexk[0] 
					                          		+ src_d[1]*u_primexk[1]
					                          		+ src_d[2]*u_primexk[2] );


				_p0 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;
			}
		}
	}
	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_3dp1( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ijk_0, ijk_1, ijk_2;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primexk;

	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	register int   iarr = D_P1; 		

	//field to be updated
	register ireal *** _p1 = s[iarr]._s3;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting multiples and parameters from medium struct */
	register ireal multi0 = himedinfo->lambda + 2.0*himedinfo->mu; //lambda + 2mu
	register ireal multi1 = himedinfo->lambda; //lambda
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;

	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp     = rickerinfo->amp;
	register ireal fpeak = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );	
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt;

	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 3D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

	/* loop in 2-direction */
	for ( index[2]=0; index[2]<nx[2]; index[2]++ ) {
		/* updating 2-component of x */
		x[2] = offset[2] + index[2]*dx[2];

		/* loop in 1-direction */
		for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
			/* updating 1-component of x */
			x[1] = offset[1] + index[1]*dx[1];

			/* loop in 0-direction */
			for ( index[0]=0; index[0]<nx[0]; index[0]++ ) {
				/* updating 0-component of x */
				x[0] = offset[0] + index[0]*dx[0];

				/* updating ijk indexes */
				ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 1; //(i,j,k)=(1,0,1)
				ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 1; //(i,j,k)=(1,1,1)
				ijk_2[0] = 1; ijk_2[1] = 2; ijk_2[2] = 1; //(i,j,k)=(1,2,1)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				_p1 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] = 
									multi0*(  src_d[0]*u_primexk[0] 
					                   	  		+ src_d[1]*u_primexk[1]
					                 	  		+ src_d[2]*u_primexk[2] );
	
				/* updating ijk indexes */
				ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 0; //(i,j,k)=(0,0,0)
				ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 0; //(i,j,k)=(0,1,0)
				ijk_2[0] = 0; ijk_2[1] = 2; ijk_2[2] = 0; //(i,j,k)=(0,2,0)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				_p1 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] +=  
									multi1*(  src_d[0]*u_primexk[0] 
					                   	   		+ src_d[1]*u_primexk[1]
					                   	   		+ src_d[2]*u_primexk[2] );

				/* updating ijk indexes */
				ijk_0[0] = 2; ijk_0[1] = 0; ijk_0[2] = 2; //(i,j,k)=(2,0,2)
				ijk_1[0] = 2; ijk_1[1] = 1; ijk_1[2] = 2; //(i,j,k)=(2,1,2)
				ijk_2[0] = 2; ijk_2[1] = 2; ijk_2[2] = 2; //(i,j,k)=(2,2,2)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							    ndim, 
							    x, 
							    src_x,
							    t,
							    alpha, 
							    beta, 
							    fpeak,
							    stderr   );
				_p1 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] +=  
									multi1*(  src_d[0]*u_primexk[0] 
					                   	   		+ src_d[1]*u_primexk[1]
					                   	   		+ src_d[2]*u_primexk[2] );


				_p1 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;
			}
		}
	}
	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_3dp2( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ijk_0, ijk_1, ijk_2;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primexk;

	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	register int   iarr = D_P2; 		

	//field to be updated
	register ireal *** _p2 = s[iarr]._s3;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting multiples and parameters from medium struct */
	register ireal multi0 = himedinfo->lambda + 2.0*himedinfo->mu; //lambda + 2mu
	register ireal multi1 = himedinfo->lambda; //lambda
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;

	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp     = rickerinfo->amp;
	register ireal fpeak = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );	
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt;

	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 2D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

	/* loop in 2-direction */
	for ( index[2]=0; index[2]<nx[2]; index[2]++ ) {
		/* updating 2-component of x */
		x[2] = offset[2] + index[2]*dx[2];

		/* loop in 1-direction */
		for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
			/* updating 1-component of x */
			x[1] = offset[1] + index[1]*dx[1];

			/* loop in 0-direction */
			for ( index[0]=0; index[0]<nx[0]; index[0]++ ) {
				/* updating 0-component of x */
				x[0] = offset[0] + index[0]*dx[0];

				/* updating ijk indexes */
				ijk_0[0] = 2; ijk_0[1] = 0; ijk_0[2] = 2; //(i,j,k)=(2,0,2)
				ijk_1[0] = 2; ijk_1[1] = 1; ijk_1[2] = 2; //(i,j,k)=(2,1,2)
				ijk_2[0] = 2; ijk_2[1] = 2; ijk_2[2] = 2; //(i,j,k)=(2,2,2)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				_p2 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] = multi0*( src_d[0]*u_primexk[0] 
								  + src_d[1]*u_primexk[1]
					                 	  + src_d[2]*u_primexk[2] );
	
				/* updating ijk indexes */
				ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 0; //(i,j,k)=(0,0,0)
				ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 0; //(i,j,k)=(0,1,0)
				ijk_2[0] = 0; ijk_2[1] = 2; ijk_2[2] = 0; //(i,j,k)=(0,2,0)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				_p2 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] +=  multi1*( src_d[0]*u_primexk[0] 
					                           + src_d[1]*u_primexk[1]
					                           + src_d[2]*u_primexk[2] );

				/* updating ijk indexes */
				ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 1; //(i,j,k)=(1,0,1)
				ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 1; //(i,j,k)=(1,1,1)
				ijk_2[0] = 1; ijk_2[1] = 2; ijk_2[2] = 1; //(i,j,k)=(1,2,1)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );

				_p2 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] += multi1*( src_d[0]*u_primexk[0] 
					                           + src_d[1]*u_primexk[1]
					                           + src_d[2]*u_primexk[2] );


				_p2 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;
			}
		}
	}
	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_3dv0( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ij_0, ij_1, ij_2;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primet;

	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	register int iarr = D_V0; 		

	//field to be updated
	register ireal ***_v0 = s[iarr]._s3;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting multiples and parameters from medium struct */
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;

	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp     = rickerinfo->amp;
	register ireal fpeak = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );	
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt;
	
	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 2D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

	/* loop in 2-direction */
	for ( index[2]=0; index[2]<nx[2]; index[2]++ ) {
		/* updating 2-component of x */
		x[2] = offset[2] + index[2]*dx[2];

		/* loop in 1-direction */
		for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
			/* updating 1-component of x */
			x[1] = offset[1] + index[1]*dx[1];

			/* loop in 0-direction */
			for ( index[0]=0; index[0]<nx[0]; index[0]++ ) {
				/* updating 0-component of x */
				x[0] = offset[0] + index[0]*dx[0];

				/* updating ij indexes */
				ij_0[0] = 0; ij_0[1] = 0; //(i,j)=(0,0)
				ij_1[0] = 0; ij_1[1] = 1; //(i,j)=(0,1)
				ij_2[0] = 0; ij_2[1] = 2; //(i,j)=(0,2)

				/* updating field */
				u_primet[0] = uij_primet( ij_0, 
							 ndim, 
							 x, 
							 src_x,
							 t,
							 alpha, 
							 beta, 
							 fpeak,
							 stderr   );
				u_primet[1] = uij_primet( ij_1, 
							 ndim, 
							 x, 
							 src_x,
							 t,
							 alpha, 
							 beta, 
							 fpeak,
							 stderr   );
				u_primet[2] = uij_primet( ij_2, 
							 ndim, 
							 x, 
							 src_x,
							 t,
							 alpha, 
							 beta, 
							 fpeak,
							 stderr   );
				_v0 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] = src_d[0]*u_primet[0] 
					        	  	  + src_d[1]*u_primet[1]
					        		  + src_d[2]*u_primet[2];


				_v0 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;
			}
		}
	}
	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_3dv1( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ij_0, ij_1, ij_2;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primet;

	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	register int iarr = D_V1; 		

	//field to be updated
	register ireal ***_v1 = s[iarr]._s3;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting multiples and parameters from medium struct */
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;

	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp     = rickerinfo->amp;
	register ireal fpeak = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt;

	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 2D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

	/* loop in 2-direction */
	for ( index[2]=0; index[2]<nx[2]; index[2]++ ) {
		/* updating 2-component of x */
		x[2] = offset[2] + index[2]*dx[2];

		/* loop in 1-direction */
		for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
			/* updating 1-component of x */
			x[1] = offset[1] + index[1]*dx[1];

			/* loop in 0-direction */
			for ( index[0]=0; index[0]<nx[0]; index[0]++ ) {
				/* updating 0-component of x */
				x[0] = offset[0] + index[0]*dx[0];

				/* updating ij indexes */
				ij_0[0] = 1; ij_0[1] = 0; //(i,j)=(1,0)
				ij_1[0] = 1; ij_1[1] = 1; //(i,j)=(1,1)
				ij_2[0] = 1; ij_2[1] = 2; //(i,j)=(1,2)

				/* updating field */
				u_primet[0] = uij_primet( ij_0, 
							 ndim, 
							 x, 
							 src_x,
							 t,
							 alpha, 
							 beta, 
							 fpeak,
							 stderr   );
				u_primet[1] = uij_primet( ij_1, 
							 ndim, 
							 x, 
							 src_x,
							 t,
							 alpha, 
							 beta, 
							 fpeak,
							 stderr   );
				u_primet[2] = uij_primet( ij_2, 
							 ndim, 
							 x, 
							 src_x,
							 t,
							 alpha, 
							 beta, 
							 fpeak,
							 stderr   );
				_v1 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] = src_d[0]*u_primet[0] 
					        		  + src_d[1]*u_primet[1]
					        		  + src_d[2]*u_primet[2];


				_v1 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;
			}
		}
	}
	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_3dv2( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ij_0, ij_1, ij_2;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primet;

	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	register int   iarr = D_V2; 		

	//field to be updated
	register ireal ***_v2 = s[iarr]._s3;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting multiples and parameters from medium struct */
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;

	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp     = rickerinfo->amp;
	register ireal fpeak = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt;

	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 2D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

	/* loop in 2-direction */
	for ( index[2]=0; index[2]<nx[2]; index[2]++ ) {
		/* updating 2-component of x */
		x[2] = offset[2] + index[2]*dx[2];

		/* loop in 1-direction */
		for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
			/* updating 1-component of x */
			x[1] = offset[1] + index[1]*dx[1];

			/* loop in 0-direction */
			for ( index[0]=0; index[0]<nx[0]; index[0]++ ) {
				/* updating 0-component of x */
				x[0] = offset[0] + index[0]*dx[0];

				/* updating ij indexes */
				ij_0[0] = 2; ij_0[1] = 0; //(i,j)=(2,0)
				ij_1[0] = 2; ij_1[1] = 1; //(i,j)=(2,1)
				ij_2[0] = 2; ij_2[1] = 2; //(i,j)=(2,2)

				/* updating field */
				u_primet[0] = uij_primet( ij_0, 
							 ndim, 
							 x, 
							 src_x,
							 t,
							 alpha, 
							 beta, 
							 fpeak,
							 stderr   );
				u_primet[1] = uij_primet( ij_1, 
							 ndim, 
							 x, 
							 src_x,
							 t,
							 alpha, 
							 beta, 
							 fpeak,
							 stderr   );
				u_primet[2] = uij_primet( ij_2, 
							 ndim, 
							 x, 
							 src_x,
							 t,
							 alpha, 
							 beta, 
							 fpeak,
							 stderr   );
				_v2 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] = src_d[0]*u_primet[0] 
					        		  + src_d[1]*u_primet[1]
					        		  + src_d[2]*u_primet[2];


				_v2 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;
			}
		}
	}
	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_3dss0( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ijk_0, ijk_1, ijk_2;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primexk;

	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	register int   iarr = D_S0; 		

	//field to be updated
	register ireal ***_ss0 = s[iarr]._s3;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting multiples and parameters from medium struct */
	register ireal mu = himedinfo->mu;
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;

	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp     = rickerinfo->amp;
	register ireal fpeak = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );	
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt;

	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 2D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

	/* loop in 2-direction */
	for ( index[2]=0; index[2]<nx[2]; index[2]++ ) {
		/* updating 2-component of x */
		x[2] = offset[2] + index[2]*dx[2];

		/* loop in 1-direction */
		for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
			/* updating 1-component of x */
			x[1] = offset[1] + index[1]*dx[1];

			/* loop in 0-direction */
			for ( index[0]=0; index[0]<nx[0]; index[0]++ ) {
				/* updating 0-component of x */
				x[0] = offset[0] + index[0]*dx[0];

				/* updating ijk indexes */
				ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 1; //(i,j,k)=(0,0,1)
				ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 1; //(i,j,k)=(0,1,1)
				ijk_2[0] = 0; ijk_2[1] = 2; ijk_2[2] = 1; //(i,j,k)=(0,2,1)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				_ss0 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] = src_d[0]*u_primexk[0] 
					          		   + src_d[1]*u_primexk[1]
					          		   + src_d[2]*u_primexk[2];
	
				/* updating ijk indexes */
				ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 0; //(i,j,k)=(1,0,0)
				ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 0; //(i,j,k)=(1,1,0)
				ijk_2[0] = 1; ijk_2[1] = 2; ijk_2[2] = 0; //(i,j,k)=(1,2,0)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				_ss0 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] += src_d[0]*u_primexk[0] 
					           		    + src_d[1]*u_primexk[1]
					           		    + src_d[2]*u_primexk[2];
				_ss0 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= mu;


				_ss0 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;
			}
		}
	}
	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_3dss1( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ijk_0, ijk_1, ijk_2;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primexk;

	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	register int   iarr = D_S1; 		
	
	//field to be updated
	register ireal ***_ss1 = s[iarr]._s3;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting multiples and parameters from medium struct */
	register ireal mu = himedinfo->mu;
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;

	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp     = rickerinfo->amp;
	register ireal fpeak = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );	
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt;

	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 2D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

	/* loop in 2-direction */
	for ( index[2]=0; index[2]<nx[2]; index[2]++ ) {
		/* updating 2-component of x */
		x[2] = offset[2] + index[2]*dx[2];

		/* loop in 1-direction */
		for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
			/* updating 1-component of x */
			x[1] = offset[1] + index[1]*dx[1];

			/* loop in 0-direction */
			for ( index[0]=0; index[0]<nx[0]; index[0]++ ) {
				/* updating 0-component of x */
				x[0] = offset[0] + index[0]*dx[0];

				/* updating ijk indexes */
				ijk_0[0] = 1; ijk_0[1] = 0; ijk_0[2] = 2; //(i,j,k)=(1,0,2)
				ijk_1[0] = 1; ijk_1[1] = 1; ijk_1[2] = 2; //(i,j,k)=(1,1,2)
				ijk_2[0] = 1; ijk_2[1] = 2; ijk_2[2] = 2; //(i,j,k)=(1,2,2)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				_ss1 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] = src_d[0]*u_primexk[0] 
					          		   + src_d[1]*u_primexk[1]
					          		   + src_d[2]*u_primexk[2];
	
				/* updating ijk indexes */
				ijk_0[0] = 2; ijk_0[1] = 0; ijk_0[2] = 1; //(i,j,k)=(2,0,1)
				ijk_1[0] = 2; ijk_1[1] = 1; ijk_1[2] = 1; //(i,j,k)=(2,1,1)
				ijk_2[0] = 2; ijk_2[1] = 2; ijk_2[2] = 1; //(i,j,k)=(2,2,1)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				_ss1 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] += src_d[0]*u_primexk[0] 
					           		    + src_d[1]*u_primexk[1]
					           	     	    + src_d[2]*u_primexk[2];
				_ss1 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= mu;


				_ss1 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;
			}
		}
	}
	return 0;
}

/*----------------------------------------------------------------------------*/
int ansol_HI_esg_3dss2( RDOM * dom, void * pars ){
/*----------------------------------------------------------------------------*/
	ANSOL_ESG_PARS * ansolpars  = (ANSOL_ESG_PARS *)(pars);
	RICKER_INFO    * rickerinfo = (RICKER_INFO *)(ansolpars->srcinfo);
	ESG_HIMED_INFO * himedinfo  = (ESG_HIMED_INFO *)(ansolpars->medinfo);
	IMODEL	       * fd_mdl     = ansolpars->link_mdl;
	RARR           * s          = dom->_s;
	RARR	       * s_phys     = (fd_mdl->ld_p)._s;

	RPNT x; 
	IPNT ijk_0, ijk_1, ijk_2;
	RPNT goff;
	RPNT offset;
	IPNT index;
	RPNT u_primexk;

	int idim;
	int ndim = ansolpars->ndim;
	int flag = 0;

	register int   iarr = D_S2; 		
	
	//field to be updated
	register ireal ***_ss2 = s[iarr]._s3;

	/* Extracting parameters from ansol struct */
	RPNT dx;
	RPNT o_coord;
	IPNT o_index;
	RPNT src_x;
	
	RASN( dx,      ansolpars->dx );
	RASN( o_coord, ansolpars->o_coord );
	IASN( o_index, ansolpars->o_index );
	RASN( src_x,   ansolpars->src_x );

	/* Extracting multiples and parameters from medium struct */
	register ireal mu = himedinfo->mu;
	register ireal rho    = himedinfo->rho;
	register ireal alpha  = himedinfo->alpha;
	register ireal beta   = himedinfo->beta;

	/* Extracting parameters from Ricker wavelet struct */
	register ireal amp     = rickerinfo->amp;
	register ireal fpeak = rickerinfo->fpeak;
	RPNT src_d;
	RASN( src_d, rickerinfo->src_d );
	
	/* Extracting time t */
	register ireal t = get_time( fd_mdl->tsind );	
	t -= ansolpars->t_off;
	t -= 0.5*ansolpars->dt;

	/* scaling factor */
	ireal sc_fac = amp/(4.0*M_PI*rho);

	/* Extracting information about grid geometry from physical domain */
	IPNT gs;	//global starting indexes
	IPNT nx;	//number of points per axis grid
	for (idim=0; idim<ndim; idim++){
		gs[idim] = s_phys[iarr]._dims[idim].gs;
		nx[idim] = s_phys[iarr]._dims[idim].n;
		
	}

	/* computing index offset due to grid type: 
		e.g., if field is of dual grid type in 0-axis then there is an 1/2 offset in 0-axis. */
	for (idim=0; idim<ndim; idim++){
		goff[idim] = 0.5 * ansolpars->gtype[iarr][idim];
	}

// fprintf(stderr,"goff = (%g,%g)\n",goff[0],goff[1]);

	//----------------------------------------------------------------------------//
	//--main spatial 2D loop------------------------------------------------------//
	RASN( x, RPNT_0 );
	for (idim=0; idim<ndim; idim++)
		offset[idim] = o_coord[idim] + ( gs[idim] + goff[idim] - o_index[idim] )*dx[idim];

	/* loop in 2-direction */
	for ( index[2]=0; index[2]<nx[2]; index[2]++ ) {
		/* updating 2-component of x */
		x[2] = offset[2] + index[2]*dx[2];

		/* loop in 1-direction */
		for ( index[1]=0; index[1]<nx[1]; index[1]++ ) {
			/* updating 1-component of x */
			x[1] = offset[1] + index[1]*dx[1];

			/* loop in 0-direction */
			for ( index[0]=0; index[0]<nx[0]; index[0]++ ) {
				/* updating 0-component of x */
				x[0] = offset[0] + index[0]*dx[0];

				/* updating ijk indexes */
				ijk_0[0] = 2; ijk_0[1] = 0; ijk_0[2] = 0; //(i,j,k)=(2,0,0)
				ijk_1[0] = 2; ijk_1[1] = 1; ijk_1[2] = 0; //(i,j,k)=(2,1,0)
				ijk_2[0] = 2; ijk_2[1] = 2; ijk_2[2] = 0; //(i,j,k)=(2,2,0)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				_ss2 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] += src_d[0]*u_primexk[0] 
					            		    + src_d[1]*u_primexk[1]
					            		    + src_d[2]*u_primexk[2];
	
				/* updating ijk indexes */
				ijk_0[0] = 0; ijk_0[1] = 0; ijk_0[2] = 2; //(i,j,k)=(0,0,2)
				ijk_1[0] = 0; ijk_1[1] = 1; ijk_1[2] = 2; //(i,j,k)=(0,1,2)
				ijk_2[0] = 0; ijk_2[1] = 2; ijk_2[2] = 2; //(i,j,k)=(0,2,2)

				/* updating field */
				u_primexk[0] = uij_primexk( ijk_0, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[1] = uij_primexk( ijk_1, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				u_primexk[2] = uij_primexk( ijk_2, 
							   ndim, 
							   x, 
							   src_x,
							   t,
							   alpha, 
							   beta, 
							   fpeak,
							   stderr   );
				_ss2 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] += src_d[0]*u_primexk[0] 
					           		    + src_d[1]*u_primexk[1]
					           	   	    + src_d[2]*u_primexk[2];
				_ss2 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= mu;


				_ss2 [ index[2]+gs[2] ] [ index[1]+gs[1] ] [ index[0]+gs[0] ] *= sc_fac;
			}
		}
	}
	return 0;
}


//#############################################################################//
//############### DEFINITION OF AUXILARY FUNCTIONS ############################//
//#############################################################################//

/*----------------------------------------------------------------------------*/
float uij_primexk( int   * indexes, 
		   int     ndim, 
		   float * x, 
		   float * xs,
		   float   t,
		   float   alpha, 
		   float   beta, 
		   float   fpeak,
		   FILE  * stream   ){
/*----------------------------------------------------------------------------*/
	if (ndim!=2 && ndim!=3){
		fprintf(stream,"Error: in uij_primexk, invalid entry for ndim: %d. ABORT\n",ndim);
		exit(1);
	}

	int idim,i,j,k;
	float r, r2, r3, r4;
	float a2, a3, b2, b3;
	float tra, trb;
	float gamma_i, gamma_j, gamma_k;
	float gamma_i_k, gamma_j_k;
	float output;
	
	i = indexes[0];
	j = indexes[1];
	k = indexes[2];

	//computing r,r3,r2,r4
	r2 = 0;
	for (idim=0; idim<ndim; idim++){
		r2 += SQR( x[idim]-xs[idim] );
	}	
	
	//case where x = xs
	if ( fabs(r2)<TOL ){    
		fprintf(stream,"ERROR: in uij_primexk, r = 0!\n");
		exit(1);
	}	
	r = sqrt(r2);
	r3 = r*r2;
	r4 = SQR(r2);

	//computing a2, a3, b2, b3
	a2 = SQR(alpha); a3 = alpha*a2;
	b2 = SQR(beta);  b3 = beta*b2;
	
	//computing tra, trb
	tra = t - r/alpha;
	trb = t - r/beta;

	//computing gamma_i, gamma_j, gamma_k
	gamma_i = (x[i]-xs[i])/r;
	gamma_j = (x[j]-xs[j])/r;
	gamma_k = (x[k]-xs[k])/r;

	//computing gamma_i_k, gamma_j_k
	gamma_i_k = (Kdelta(i,k)-gamma_i*gamma_k)/r;
	gamma_j_k = (Kdelta(j,k)-gamma_j*gamma_k)/r;

	//first term
	output = ( 3.0*(gamma_i_k*gamma_j + gamma_i*gamma_j_k)/r3 - gamma_k*(3.0*gamma_i*gamma_j-Kdelta(i,j))/r4 )*Int2( t, r, alpha, beta, fpeak );
        output -= gamma_k*(3.0*gamma_i*gamma_j-Kdelta(i,j))*Int4( t, r, alpha, beta, fpeak)/r4;

	//second term
	output += ( -gamma_i*gamma_j*gamma_k/r + gamma_i_k*gamma_j + gamma_i*gamma_j_k )*comprick(tra,fpeak)/(r*a2);
	output -= gamma_i*gamma_j*gamma_k*Ricker_prime(tra,fpeak)/(r*a3);

	//third term
	output += ( gamma_k*(gamma_i*gamma_j-Kdelta(i,j))/r - gamma_i_k*gamma_j - gamma_i*gamma_j_k )*comprick(trb,fpeak)/(r*b2);
	output += gamma_k*(gamma_i*gamma_j-Kdelta(i,j))*Ricker_prime(trb,fpeak)/(r*b3);

	return output;
}

/*----------------------------------------------------------------------------*/
float uij_primet( int   * indexes, 
		  int     ndim, 
		  float * x, 
		  float * xs,
		  float   t,
		  float   alpha, 
		  float   beta, 
		  float   fpeak,
         	  FILE  * stream  ){
/*----------------------------------------------------------------------------*/
	if (ndim!=2 && ndim!=3){
		fprintf(stream,"Error: in uij_primexk, invalid entry for ndim: %d. ABORT\n",ndim);
		exit(1);
	}

	int idim,i,j;
	float r, r2, r3;
	float tra, trb;
	float gamma_i, gamma_j;
	float a2, b2;
	float output;
	
	i = indexes[0];
	j = indexes[1];

	//computing r,r2,r3
	r2 = 0;
	for (idim=0; idim<ndim; idim++){
		r2 += SQR( x[idim]-xs[idim] );
	}	

	//case where x = xs	
	if ( fabs(r2)<TOL ){    
		fprintf(stream,"ERROR: in uij_primexk, r = 0!\n");
		exit(1);
	}	
	r = sqrt(r2);
	r3 = r*r2;

	//computing tra, trb, a2, b2
	tra = t - r/alpha;
	trb = t - r/beta;
	a2 = SQR(alpha);
	b2 = SQR(beta);

	//computing gamma_i, gamma_j
	gamma_i = (x[i]-xs[i])/r;
	gamma_j = (x[j]-xs[j])/r;

	//first term
	output = ( 3.0*gamma_i*gamma_j - Kdelta(i,j) )*Int3( t, r, alpha, beta, fpeak )/r3;

	//second term
	output += gamma_i*gamma_j*Ricker_prime(tra,fpeak)/(a2*r);

	//third term
	output -= ( gamma_i*gamma_j - Kdelta(i,j) )*Ricker_prime(trb,fpeak)/(b2*r);

	return output;
}

/*----------------------------------------------------------------------------*/
float Int1( float t, 
            float r, 
            float alpha, 
            float beta, 
            float fpeak ){
/*----------------------------------------------------------------------------*/
	float tra, trb, output;
	tra = t-r/alpha; 
	trb = t-r/beta;

	output = -compdgauss( trb, fpeak ) + compdgauss( tra, fpeak );
	return output;
}

/*----------------------------------------------------------------------------*/
float Int2( float t, 
	    float r, 
	    float alpha, 
	    float beta, 
	    float fpeak ){
/*----------------------------------------------------------------------------*/
	float tra, trb, output;
	tra = t-r/alpha; 
	trb = t-r/beta;

	output = (trb-t)*compdgauss(trb,fpeak) - (tra-t)*compdgauss(tra,fpeak);
	output += -compgauss(trb,fpeak) + compgauss(tra,fpeak);

	return output;
}

/*----------------------------------------------------------------------------*/
float Int3( float t, 
	    float r, 
	    float alpha, 
	    float beta, 
	    float fpeak ){
/*----------------------------------------------------------------------------*/
	float tra, trb, output;
	tra = t-r/alpha; 
	trb = t-r/beta;

	output = (trb-t)*comprick(trb,fpeak) - (tra-t)*comprick(tra,fpeak);
	output += -compdgauss(trb,fpeak) + compdgauss(tra,fpeak);

	return output;
}

/*----------------------------------------------------------------------------*/
float Int4( float t, 
	    float r, 
	    float alpha, 
	    float beta, 
	    float fpeak ){
/*----------------------------------------------------------------------------*/
	float tra, trb, output;
	tra = t-r/alpha; 
	trb = t-r/beta;

	output = -SQR(trb-t)*comprick(trb,fpeak) + SQR(tra-t)*comprick(tra,fpeak);
	output += 2.0*(trb-t)*compdgauss(trb,fpeak) - 2.0*(tra-t)*compdgauss(tra,fpeak);
	output += -2.0*compgauss(trb,fpeak) + 2.0*compgauss(tra,fpeak);

	return output;
}

/*----------------------------------------------------------------------------*/
float Ricker_prime( float t,
	      	    float fpeak  ){
/*----------------------------------------------------------------------------*/
	float st = SQR(M_PI)*SQR(fpeak);
	float output = (-6.0*t + 4.0*st*SQR(t)*t)*st*exp(-st*SQR(t));

	return output;
}


/*----------------------------------------------------------------------------*/
float Kdelta( int i, int j ){
/*----------------------------------------------------------------------------*/
	if( i==j ) return 1.0;
	return 0.0;
}
