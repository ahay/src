/* 3D Kirchhoff Migration Interface */
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
#include "mig3.h"
#ifdef _OPENMP
#include <omp.h>
#endif

static int nt, nx, ny, n12, apt;
static float ot, ox, oy, dt, dx, dy, *vel, rho, angle;
static char antialias;
static bool omp;

void testmig3_init(int n1, int n2, int n3  /* data size */, 
		 float d1, float d2, float d3 /* sampling */,
                 float o1, float o2, float o3 /* origins */,
                 float *fvel /* velocity for Kirchhoff */,
                 float frho /* phase factor for Kirchhoff */,
                 char fantialias /* Kirchoff antialiasing method */,
	         bool fomp /* OpenMP */,
		 int fapt, float fangle /* apertures */)
/*< Initialize Chain of Operators  >*/
{

    omp = fomp;

    apt = fapt;
    angle = fangle;

    /* Copying Dimensions */
    nt = n1;
    nx = n2;
    ny = n3;
    ot = o1;
    dt = d1;
    ox = o2;
    dx = d2;
    oy = o3;
    dy = d3;

    /* Copying Kirchhoff Parameters */
    vel = fvel;
    rho = frho;
    antialias = fantialias;

    /* Allocating Memory */
    n12  = nt*nx*ny; 		

}

void testmig3_lop (bool adj, bool add, int fnx, int fny, float* x, float* y)
/*< Apply Kirchhoff Migration/Modeling >*/
{

    int i, ch=0;

    if((nt*nx*ny != n12) || (fnx != n12) || (fny != n12)) sf_error("Wrong dimensions");

    /* ADJNULL OPERATOR IS APPLIED IN MIG3LOP */

    #ifdef _OPENMP
    float tstart;
    tstart = omp_get_wtime();
#endif

    mig3_lop (adj,add, nt,nx,ny, dt,dx,dy, ot,ox,oy,
		y /* zero-offset */,
		x /* image */,
		vel /* velocity */,
		rho /* phase factor */,
		antialias /* antialiasing method */,
		omp /* OpenMP */,
		apt, angle /* apertures */);

    if(adj==false){/* check for problems in y */

		for(i=0; i<n12; i++){

			/* this is not the biggest float - checks if we approach it */
			if( ((fabsf(y[i]) > 3.0e+38)  || (isinf(y[i])) || (isnan(y[i]))) && !ch ) { sf_warning("testmig3 adj=false: problems in the zero-offset"); ch=1; }

		}/* for */

    }/* if */
	
    if(adj==true){/* check for problems in x */	

	for(i=0; i<n12; i++){

			if( ((fabsf(x[i]) > 3.0e+38)  || (isinf(x[i])) || (isnan(x[i]))) && !ch ) { sf_warning("testmig3 adj=true: problems in the image"); ch=1; }

		}/* for */

    }/* if */

    #ifdef _OPENMP
    float tend, elapsed;
    tend = omp_get_wtime();
    elapsed = tend - tstart;
    sf_warning("testmig3: elapsed=%f [s]",elapsed);
    #endif

}
