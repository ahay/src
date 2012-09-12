/* 2D Helmholtz solver by LU factorization. */
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
    int n1, n2, npw, npml, pad1, pad2, is, ns;
    SuiteSparse_long n, nz, *Ti, *Tj;
    float d1, d2, **v, freq, eps;
    double omega, *Tx, *Tz;
    SuiteSparse_long status, *Ap, *Ai, *Map;
    double *Ax, *Az, *Xx, *Xz, *Bx, *Bz;
    void *Symbolic, *Numeric;
    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
    sf_complex **f;
    char *save;
    sf_file in, out, source;
 
    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");
   
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */
    
    /* read input dimension */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input.");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input.");

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input.");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input.");

    /* read input */
    v = sf_floatalloc2(n1,n2);
    sf_floatread(v[0],n1*n2,in);
    
    if (!sf_getfloat("freq",&freq)) freq=2.5;
    /* frequency (Hz) */

    omega = (double) 2.*SF_PI*freq;

    if (!sf_getint("npw",&npw)) npw=6;
    /* number of points per wave-length */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* epsilon for PML */

    if (verb) sf_warning("Preparing...");

    npml = npw*2;
    pad1 = n1+2*npml;
    pad2 = n2+2*npml;

    n = (pad1-2)*(pad2-2);
    nz = 5*(pad1-2)*(pad2-2)-2*(pad1-4)-2*(pad2-4)-8;

    /* assemble matrix in triplet form */
    Ti = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
    Tj = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
    Tx = (double*) sf_alloc(nz,sizeof(double));
    Tz = (double*) sf_alloc(nz,sizeof(double));

    fdprep(omega, eps, 
	   n1, n2, d1, d2, v,
	   npml, pad1, pad2, n, nz, 
	   Ti, Tj, Tx, Tz);

    /* convert triplet to compressed-column form */
    Ap = (SuiteSparse_long*) sf_alloc(n+1,sizeof(SuiteSparse_long));
    Ai = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
    Map = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));

    Ax = (double*) sf_alloc(nz,sizeof(double));
    Az = (double*) sf_alloc(nz,sizeof(double));

    status = umfpack_zl_triplet_to_col (n, n, nz, 
					Ti, Tj, Tx, Tz, 
					Ap, Ai, Ax, Az, Map);

    free(Ti); free(Tj); free(Tx); free(Tz);

    if (verb) sf_warning("LU factorizing...");

    /* prepare LU */
    umfpack_zl_defaults (Control);
        
    /* LU factorization */
    status = umfpack_zl_symbolic (n, n, 
				  Ap, Ai, Ax, Az, 
				  &Symbolic, Control, Info);

    status = umfpack_zl_numeric (Ap, Ai, Ax, Az, 
				 Symbolic, &Numeric, 
				 Control, Info);

    (void) umfpack_zl_free_symbolic (&Symbolic);
    
    /* save Numeric */
    save = sf_getstring("lu");    
    if (save != NULL) 
	status = umfpack_zl_save_numeric (Numeric, save);
    
    if (verb) sf_warning("Solving...");

    /* solve linear system */
    if (NULL == sf_getstring("source"))
	sf_error("Need source=");
    source = sf_input("source");

    if (!sf_histint(source,"n3",&ns)) ns=1;

    sf_settype(out,SF_COMPLEX);
    sf_putint(out,"n3",ns);    

    f = sf_complexalloc2(n1,n2);
    Bx = (double*) sf_alloc(n,sizeof(double));
    Bz = (double*) sf_alloc(n,sizeof(double));
    Xx = (double*) sf_alloc(n,sizeof(double));
    Xz = (double*) sf_alloc(n,sizeof(double));

    for (is=0; is < ns; is++) {
	sf_complexread(f[0],n1*n2,source);

	fdpad(npml,pad1,pad2, f,Bx,Bz);
			
	status = umfpack_zl_solve (UMFPACK_A, 
				   Ap, Ai, Ax, Az, 
				   Xx, Xz, Bx, Bz, 
				   Numeric, Control, Info);


	fdcut(npml,pad1,pad2, f,Xx,Xz);

	sf_complexwrite(f[0],n1*n2,out);
    }

    (void) umfpack_zl_free_numeric (&Numeric);

    exit(0);
}
