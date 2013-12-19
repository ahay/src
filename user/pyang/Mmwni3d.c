/* Minimum weighted-norm interpolation (MWNI) with shaping,
 implemented with conjugate gradient least squares (CGLS) method using
 shaping regularization
*/
/*
  Copyright (C) 2013 Pengliang Yang, Xi'an Jiaotong University, UT Austin

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
#include <rsf.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mwni3d.h"

int main(int argc, char* argv[])
{
    /* input and output variables */
    bool verb;
    int niter,dim, n1, n2, n3, rect[3], ndat[3];
    char key[6];
    float tol, *din, *mask;
    sf_file Fin, Fout, Fmask;/* mask and I/O files*/ 

    sf_init(argc,argv);	/* Madagascar initialization */
#ifdef _OPENMP
    omp_init(); 	/* initialize OpenMP support */
#endif

    /* setup I/O files */
    Fin=sf_input("in");	/* read the data to be interpolated */
    Fmask=sf_input("mask");  	/* read the 2-D mask for 3-D data */
    Fout=sf_output("out"); 	/* output the reconstructed data */
 
    if(!sf_getbool("verb",&verb))    	verb=false;
    /* verbosity */
    if (!sf_getint("niter",&niter)) 	niter=100;
    /* total number iterations */
    if (!sf_getfloat("tol",&tol)) 	tol=1.0e-6;
    /* iteration tolerance */

    /* Read the data size */
    if (!sf_histint(Fin,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(Fin,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(Fin,"n3",&n3)) sf_error("No n3= in input");
    ndat[0]=n1; ndat[1]=n2; ndat[2]=n3;   dim=3;
    for (int i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
    }

    /* allocate data and mask arrays */
    din=sf_floatalloc(n1*n2*n3); sf_floatread(din,n1*n2*n3,Fin);
    if (NULL != sf_getstring("mask")){
	mask=sf_floatalloc(n2*n3);
	sf_floatread(mask, n2*n3, Fmask);
    }

    mwni3d_init(dim, ndat, rect, mask, tol, verb);
    mwni3d_inversion(din, din, niter);
    mwni3d_close();

    sf_floatwrite(din, n1*n2*n3, Fout); /* output reconstructed seismograms */

    sf_close();
    exit(0);
}
