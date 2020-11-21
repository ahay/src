/* Chain of Operators Interface */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <math.h>

#include <rsf.h>
/*^*/

#include "tdpi.h"
#include "t2warp.h"
#include "mig3.h"
#include "azpwd.h"
#ifdef _OPENMP
#include <omp.h>
#endif

static int nt2, nt, nx, ny, n12, nt12, apt;
static float ot, ox, oy, dt, dx, dy, v_1, v_2, v_3, v_4, eps, epst2, passthr, *vel, rho, angle;
static float *model, *datat2, *outputt2, *data, *output, *pwddata;
static char antialias;
static bool domod, sm, dopi, omp;
static float ot2, dt2;

void piintface_init(int n1, int n2, int n3  /* data size */, 
		 float d1, float d2, float d3 /* sampling */,
                 float o1, float o2, float o3 /* origins */,
		 float fpassthr /* pass threshold */, 
		 float v1 /* left cut velocity */,
		 float v2 /* left full pass velocity */,
		 float v3 /* right full pass velocity */,
		 float v4 /* right cut velocity  */,
		 float feps /* damper for pi */,
                 float fepst2 /* damper for t2warp */,
                 int fnt2 /* new axis length for t2warp */,
                 float *fvel /* velocity for Kirchhoff */,
                 float frho /* phase factor for Kirchhoff */,
                 char fantialias /* Kirchoff antialiasing method */,
                 int nw /* [1,2,3] accuracy order for PWD */,
		 float *dip1, float *dip2 /* dip distribution */,
		 float *faz /* azimuth distribution */,
                 int nj1, int nj2 /* antialiasing for PWD */,
		 bool fdomod, bool fsm, bool fdopi /* disable operator flags */,
	         bool fomp /* OpenMP */,
		 int fapt, float fangle)
/*< Initialize Chain of Operators  >*/
{

    /* disable operator flags */
    domod = fdomod;
    sm = fsm;
    dopi = fdopi;

    omp = fomp;

    apt = fapt;
    angle = fangle;

    /* Copying Dimensions */
    nt = n1;
    nt2 = fnt2;
    nx = n2;
    ny = n3;
    ot = o1;
    dt = d1;
    ox = o2;
    dx = d2;
    oy = o3;
    dy = d3;

    /* Copying Path-Summation Integral Parameters */
    v_1 = v1;
    v_2 = v2;
    v_3 = v3;
    v_4 = v4;
    eps = feps;
    passthr = fpassthr;

    epst2 = fepst2;

    /* Copying Kirchhoff Parameters */
    vel = fvel;
    rho = frho;
    antialias = fantialias;

    /* Allocating Memory */
    n12  = nt*nx*ny; 	
    nt12 = nt2*nx*ny;   

    data = sf_floatalloc(n12);
    model = sf_floatalloc(n12);
    datat2 = sf_floatalloc(nt12); 
    outputt2 = sf_floatalloc(nt12);
    output = sf_floatalloc(n12);
    pwddata = sf_floatalloc(n12);

    /* t2warping axis evaluation */
    ot2 = ot*ot;
    dt2 = ot+(nt-1)*dt;
    dt2 = (dt2*dt2 - ot2)/(nt2-1);	
		
    /* take in account different output trace length */
    t2warp_init(nt,nt2,nx*ny,ot,dt,ot2,dt2,epst2);

    /* init pi filter */
    tdpi_init(nt2, nx,ny, dt2, dx,dy, 0.001, v_1, v_2, v_3, v_4, eps);

    /* AzPWD init */
    azpwd_init(nt,nx,ny,dt,dx,dy,ot,ox,oy,
		 nw /* [1,2,3] accuracy order */,
		 dip1,dip2 /* dip distribution */,
		 faz /* azimuth distribution */,
                 nj1,nj2 /* antialiasing */);

}

void piintface_lop (bool adj, bool add, int fnx, int fny, float* x, float* y)
/*< Apply Chain of Operators >*/
{

    int i, ch=0;

    if((nt*nx*ny != n12) || (fnx != n12) || (fny != n12)) sf_error("Wrong dimensions");

    sf_adjnull(adj,add,fnx,fny,x,y);

    /* zeroing intermediate allocations */
    for (i=0; i < n12; i++){

	data[i] = 0.;
	model[i] = 0.;
	datat2[i] = 0.;
	outputt2[i] = 0.;
	output[i] = 0.;
	pwddata[i] = 0.;

    }

    if(n12 != nt12){

    	for (i=n12; i<nt12; i++){

		datat2[i] = 0.;
		outputt2[i] = 0.;		

	}

    }

    #ifdef _OPENMP
    float tstart;
    tstart = omp_get_wtime();
    #endif

    /* addition flag comment:
       we are chaining up 3 operators,
       which when combined can be treated as a single one.
       Three operators in adj=false mode:
       1) Kirchhoff (add=false)
       2) AzPWD (add=false)
       3) Path Summation (in chain3 add=true):
          a) t2warp
          b) Path-Summation filtering (with ffts inside)
	  c) it2warp 
       This flag distribution will guarantee
       that the initial guess residual will be
       added after all operators in the chain 
       have been executed.
       Therefore for adjoint:
       1) Path-Summation filtering (add=false)
       2) AzPWD (add=false)
       3) Kirchhoff (add=true). */
   
    if (adj == false){

        /* if perform Kirchhoff modeling (from x to data) */
	if (domod){

		mig3_lop (false,false, nt,nx,ny, dt,dx,dy, ot,ox,oy,
		data /* zero-offset */,
		x /* image */,
		vel /* velocity */,
		rho /* phase factor */,
		antialias /* antialiasing method */,
		omp, apt, angle);

		for(i=0; i<n12; i++){

			if( ((fabsf(data[i]) > 3.0e+38) || (isinf(data[i])) || (isnan(data[i]))) && !ch )
			{ 
				sf_warning("piintface.c adj=false: unstable after Kirchhoff modeling"); ch=1;
			}

		}/* debugging for */

        } else {/* copy data otherwise */

		for(i=0; i<n12; i++){

			data[i] = x[i];

		}/* for */

	}/* else */

	/* apply AzPWD (data to data in place) */
	if (sm) {

		azpwd_lop (false,false,n12,n12,data,pwddata);
		
		for(i=0; i<n12; i++){

			data[i] = pwddata[i];

			if( ((fabsf(data[i]) > 3.0e+38) || (isinf(data[i])) || (isnan(data[i]))) && !ch )
			{ 
				sf_warning("piintface.c adj=false: unstable after AzPWD"); ch=1;
			}

		}/* for */

	}/* AzPWD */

    }/* adj=false */
	
    if (dopi) {/* if filter by path-summation integration operator (from y to output) */

    	if(adj == true) {

		sf_chain3(t2warp_inv,tdpi_lop,t2warp,adj,false,nt*nx*ny,nt2*nx*ny,nt2*nx*ny,nt*nx*ny,output,y,outputt2,datat2);

		for(i=0; i<n12; i++){

			if( ((fabsf(output[i]) > 3.0e+38) || (isinf(output[i])) || (isnan(output[i]))) && !ch )
			{
				sf_warning("piintface.c adj=true: unstable after path summation"); ch=1;
			}/* if */

		}/* for */

    	} else {/* from data to y */
	
		sf_chain3(t2warp_inv,tdpi_lop,t2warp,false,add,nt*nx*ny,nt2*nx*ny,nt2*nx*ny,nt*nx*ny,data,y,datat2,outputt2);
		
		for(i=0; i<n12; i++){

			if( ((fabsf(y[i]) > 3.0e+38) || (isinf(y[i])) || (isnan(y[i]))) && !ch )
			{
				sf_warning("piintface.c adj=false: unstable after path summation"); ch=1;
			}/* if */

		}/* for */

    	}/* else */

    } else {/* just copy */

	if (adj){

		for(i=0; i<n12; i++){

			output[i] = y[i];

		}/* for */

	} else {

		/* if path-summation filtering is disabled addition flag
		   is taken care of here */
		for(i=0; i<n12; i++){

			y[i] += data[i];

		}/* for */

	}/* adj */

    }/* simply copying */

    if (adj == true) {

	if(sm) {/* adjoint AzPWD (output to output in place) */

		azpwd_lop (true,false,n12,n12,pwddata,output);
		
		for(i=0; i<n12; i++){

			output[i] = pwddata[i];

			if( ((fabsf(output[i]) > 3.0e+38) || (isinf(output[i])) || (isnan(output[i]))) && !ch )
			{
				sf_warning("piintface.c adj=true: unstable after AzPWD"); ch=1;
			}

		}/* for */

	}/* sm */

	if(domod) {/* adjoint Kirchhoff modeling = migration (from output to x)*/

		mig3_lop (true,add, nt,nx,ny, dt,dx,dy, ot,ox,oy,
		output /* zero-offset */,
		x /* image */,
		vel /* velocity */,
		rho /* phase factor */,
		antialias /* antialiasing method */,
		omp, apt, angle);

		for(i=0; i<n12; i++){

			if( ((fabsf(x[i]) > 3.0e+38) || (isinf(x[i])) || (isnan(x[i]))) && !ch )
			{
				sf_warning("piintface.c adj=true: unstable after Kirchhoff migration"); ch=1;
			}

		}

	} else {
	
		for(i=0; i<n12; i++){

			x[i] += output[i];

		}

	}/* else */
    	

    }/* adj */

    #ifdef _OPENMP
    float tend, elapsed;
    tend = omp_get_wtime();
    elapsed = tend - tstart;
    sf_warning("elapsed=%f [s]",elapsed);
    #endif

}

void piintface_close(void)
/*< free allocated storage >*/
{
    tdpi_close();
    t2warp_close();
    azpwd_close();
    free(data);
    free(model);
    free(datat2);
    free(outputt2);
    free(output);
    free(pwddata);
    sf_warning("piintface: memory free");
}
