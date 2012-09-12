/* 2D frequency-domain migration with extend imaging condition. */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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
#include <umfpack.h>

#include "fdprep.h"

int main(int argc, char* argv[])
{
    bool verb;
    int npw, npml, pad1, pad2, n1, n2; 
    int ih, nh, is, ns, iw, nw, i, j;
    SuiteSparse_long n, nz, *Ti, *Tj;
    float eps, d1, d2, **vel, **image, dw, ow;
    double omega, *Tx, *Tz;
    SuiteSparse_long status, *Ap, *Ai, *Map;
    double *Ax, *Az, *Xx, *Xz, *Bx, *Bz;
    void *Symbolic, *Numeric;
    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
    sf_complex **srce, **recv;
    sf_file in, out, source, data;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getint("nh",&nh)) nh=1;
    /* horizontal space-lag */

    if (!sf_getint("npw",&npw)) npw=6;
    /* number of points per wave-length */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* epsilon for PML */

    /* read model */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input.");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input.");

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input.");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input.");

    vel = sf_floatalloc2(n1,n2);
    sf_floatread(vel[0],n1*n2,in);

    /* read source */
    if (NULL == sf_getstring("source"))
	sf_error("Need source=");
    source = sf_input("source");

    if (!sf_histint(source,"n3",&ns)) sf_error("No ns=.");
    if (!sf_histint(source,"n4",&nw)) sf_error("No nw=.");
    if (!sf_histfloat(source,"d4",&dw)) sf_error("No dw=.");
    if (!sf_histfloat(source,"o4",&ow)) sf_error("No ow=.");    

    srce = sf_complexalloc2(n1,n2);

    /* read receiver */
    if (NULL == sf_getstring("data"))
	sf_error("Need data=");
    data = sf_input("data");

    recv = sf_complexalloc2(n1,n2);

    /* write output header */
    sf_putint(out,"n3",2*nh-1);
    sf_putint(out,"n4",nw);

    /* allocate temporary memory */
    npml = npw*2;
    pad1 = n1+2*npml;
    pad2 = n2+2*npml;

    n = (pad1-2)*(pad2-2);
    nz = 5*(pad1-2)*(pad2-2)-2*(pad1-4)-2*(pad2-4)-8;

    Ti = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
    Tj = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
    Tx = (double*) sf_alloc(nz,sizeof(double));
    Tz = (double*) sf_alloc(nz,sizeof(double));

    Ap = (SuiteSparse_long*) sf_alloc(n+1,sizeof(SuiteSparse_long));
    Ai = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
    Map = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));

    Ax = (double*) sf_alloc(nz,sizeof(double));
    Az = (double*) sf_alloc(nz,sizeof(double));

    Bx = (double*) sf_alloc(n,sizeof(double));
    Bz = (double*) sf_alloc(n,sizeof(double));
    Xx = (double*) sf_alloc(n,sizeof(double));
    Xz = (double*) sf_alloc(n,sizeof(double));

    image = sf_floatalloc2(n1,n2);

    umfpack_zl_defaults (Control);

    /* loop over frequency */
    for (iw=0; iw < nw; iw++) {
	omega = (double) 2.*SF_PI*(ow+iw*dw);

	if (verb) sf_warning("Frequency %d of %d.",iw+1,nw);

	/* assemble matrix */
	fdprep(omega, eps, 
	       n1, n2, d1, d2, vel,
	       npml, pad1, pad2, n, nz, 
	       Ti, Tj, Tx, Tz);

	status = umfpack_zl_triplet_to_col (n, n, nz, 
					    Ti, Tj, Tx, Tz, 
					    Ap, Ai, Ax, Az, Map);

	/* LU */
	status = umfpack_zl_symbolic (n, n, 
				      Ap, Ai, Ax, Az, 
				      &Symbolic, Control, Info);

	status = umfpack_zl_numeric (Ap, Ai, Ax, Az, 
				     Symbolic, &Numeric, 
				     Control, Info);

	/* loop over shots */
	for (is=0; is < ns; is++) {

	    /* source wavefield */
	    sf_complexread(srce[0],n1*n2,source);

	    fdpad(npml,pad1,pad2, srce,Bx,Bz);

	    status = umfpack_zl_solve (UMFPACK_A, 
				       Ap, Ai, Ax, Az, 
				       Xx, Xz, Bx, Bz, 
				       Numeric, Control, Info);

	    fdcut(npml,pad1,pad2, srce,Xx,Xz);

	    /* receiver wavefield */
	    sf_complexread(recv[0],n1*n2,data);

	    fdpad(npml,pad1,pad2, recv,Bx,Bz);

	    status = umfpack_zl_solve (UMFPACK_A, 
				       Ap, Ai, Ax, Az, 
				       Xx, Xz, Bx, Bz, 
				       Numeric, Control, Info);

	    fdcut(npml,pad1,pad2, recv,Xx,Xz);

	    /* imaging condition */
	    for (ih=1-nh; ih < nh; ih++) {
		for (j=0; j < n2; j++) {
		    for (i=0; i < n1; i++) {
			if (j-abs(ih) < 0 || j+abs(ih) >= n2) {
			    image[j][i] = 0.;
			} else {
			    image[j][i] = creal(conjf(srce[j-ih][i])*recv[j+ih][i]);
			}
		    }
		}

		sf_floatwrite(image[0],n1*n2,out);
	    }
	}
    }

    exit(0);
}
