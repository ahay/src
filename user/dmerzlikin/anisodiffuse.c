/* Anisotropic diffusion by regularized inversion. */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
#include <assert.h>

#include <rsf.h>

static float *vx, *vy, eps, *vxtemp, *vytemp, *datatemp, *datatemp2;
static int niter, repeat;
static int n1, n2, n3, n12, n123, n23;

void anisogradient (bool adj, bool add, 
		      int nx, int ng, float* x, float* g)
/*< gradient operator not for external use >*/
{
    int i2,i3,i;

    //sf_warning("in the gradient x[100]=%g g[100]=%g",x[100],g[100]);

    assert(n23 == n3*n2 && nx == n23 && ng == n23);

    sf_adjnull (adj,add,nx,ng,x,g);

    /* currently simply browsing through the third dimension */
    //for (i3=0; i3 < n3; i3++) {

    	for (i3=1; i3 < n3-1; i3++) {  

		for (i2=1; i2 < n2-1; i2++) {
	    		
			i = i2+i3*n2;//+i3*n1*n2; /* map 2-D to 1-D */

	    		if (adj) {

				x[i+1]  += 0.5*vxtemp[i]*g[i];
	        		x[i-1]  -= 0.5*vxtemp[i]*g[i];
	        		x[i+n2] += 0.5*vytemp[i]*g[i];
	        		x[i-n2] -= 0.5*vytemp[i]*g[i];

	    		} else {

				g[i]     += 0.5*vxtemp[i]*(x[i+1] - x[i-1]) + 0.5*vytemp[i]*(x[i+n2] - x[i-n2]);

	    		}

		}/* i1 */

    	}/* i2 */

    //}/* i3 */    

    //sf_warning("done with the gradient");

}

void anisodiffuse_init(int nt, int nx, int ny /* data size */, 
		 float *fvx, float *fvy /* parallel to edges vector components */,
		 int fniter /* number of iterations */,
		 int frepeat /* number of smoothing iterations */,
		 float feps /* regularization */)
/*< initialize anisotropic diffusion >*/
{

    vx = fvx;
    vy = fvy;

    niter = fniter;
    repeat = frepeat;
  
    eps = feps;

    n1 = nt;
    n2 = nx;
    n3 = ny;
    n12 = n1*n2;
    n123 = n1*n2*n3;
    n23 = n2*n3;

    datatemp = sf_floatalloc(n23);
    datatemp2 = sf_floatalloc(n23);
    vxtemp = sf_floatalloc(n23);
    vytemp = sf_floatalloc(n23);

}

void anisodiffuse_lop(int nx, int ny, float* x, float* y)
/*< anisotropic diffusion loop >*/
{

    int i, i1, i2, i3;

    for (i1 = 0; i1 < n23; i1++){

	datatemp[i1] = 0.0;
	datatemp2[i1] = 0.0;
	vxtemp[i1] = 0.0;
	vytemp[i1] = 0.0;

    }

    /*for (i1 = 0; i1 < n123; i1++){

	y[i1] = 0.0;

    }*/


    //sf_warning("before repeat");

    for (i=0; i < repeat; i++) {
	/* Call a generic function for solving
	   .   F[m] =~ d
	   eps*D[m] =~ 0 
	*/

	for (i1 = 0; i1 < n1; i1++){


		//sf_warning("before resorting");
		
		/* resorting the data to apply anisotropic diffusion slice by slice */
		for (i2 = 0; i2 < n2; i2++){

			for (i3 = 0; i3 < n3; i3++){

				datatemp[i2 + i3*n2] = x[i1 + i2*n1 + i3*n1*n2];

				vxtemp[i2 + i3*n2] = vx[i1 + i2*n1 + i3*n1*n2];

				vytemp[i2 + i3*n2] = vy[i1 + i2*n1 + i3*n1*n2];

			}/* i1 */

		}/* i2 */
	
		//sf_warning("before solver");

		sf_solver_reg(sf_copy_lop /* F */,
			      sf_cgstep   /* CG step */,
			      anisogradient    /* D */,
			      n23         /* size of D[m] */,
			      n23         /* size of m */,
			      n23         /* size of d */,
			      datatemp2        /* m (output) */,
			      datatemp        /* d (input)  */,
			      niter       /* # iterations */,
			      eps         /* eps */,
			      "verb",false,"end");

		sf_cgstep_close(); /* free internal storage */

		//sf_warning("before resorting back");

		/* resorting back */
		for (i2 = 0; i2 < n2; i2++){

			for (i3 = 0; i3 < n3; i3++){

				y[i1 + i2*n1 + i3*n1*n2] = datatemp2[i2 + i3*n2];

				datatemp2[i2 + i3*n2] = 0.0;

			}/* i3 */

		}/* i2 */

	}/* time-slice loop */

    }/* repeat loop */

}

void anisodiffuse_close(void)
/*< free allocated storage >*/
{
    free(vxtemp);
    free(vytemp);
    free(datatemp);
    free(datatemp2);
    sf_warning("anisodiffuse: memory free");
}
