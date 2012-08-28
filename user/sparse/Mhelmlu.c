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

double maxvel(int nm, float *vel)
/* find maximum velocity */
{
    int i;
    float val;
    
    val = vel[0];
    for (i=1; i < nm; i++) {
	if (vel[i] > val) val = vel[i];
    }

    return ((double) val);
}

int main(int argc, char* argv[])
{
    int n1, n2, npw, npml, pad1, pad2, i, j;
    int n_row, n_col, nz, count, index, *Ti, *Tj;
    float d1, d2, **v, **f, freq, eps;
    double omega, eta1, eta2, mvel, c1, c2;
    double *g1, *g2, **pad, *Tx, *Tz;
    double complex *s1, *s2, neib, cent;
    int status, *Ap, *Ai, *Map;
    double *Ax, *Az, *Xx, *Xz, *Bx, *Bz;
    void *Symbolic, *Numeric;
    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
    sf_complex *utemp;
    sf_file in, out, source;
 
    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");
   
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

    if (!sf_getint("npw",&npw)) npw=5;
    /* number of points per wave-length */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* epsilon for PML */

    /* prepare PML */
    npml = npw*2;
    pad1 = n1+2*npml;
    pad2 = n2+2*npml;

    eta1 = (double) npml*d1;
    eta2 = (double) npml*d2;

    mvel = maxvel(n1*n2,v[0]);
    c1 = -3.*log((double) eps)*mvel/(eta1*omega);
    c2 = -3.*log((double) eps)*mvel/(eta2*omega);

    g1 = (double*) sf_alloc(pad1,sizeof(double));
    for (i=0; i < pad1; i++) {
	if (i < npml) {
	    g1[i] = pow(c1*((npml-i)*d1/eta1),2.);
	} else if (i >= pad1-npml) {
	    g1[i] = pow(c1*((i-(pad1-npml-1))*d1/eta1),2.);
	} else {
	    g1[i] = 0.;
	}
    }

    s1 = (double complex*) sf_alloc(pad1,sizeof(double complex));
    for (i=0; i < pad1; i++) {
	s1[i] = 1./(1.+I*g1[i]);
    }

    g2 = (double*) sf_alloc(pad2,sizeof(double));
    for (j=0; j < pad2; j++) {
	if (j < npml) {
	    g2[j] = pow(c2*((npml-j)*d2/eta2),2.);
	} else if (j >= pad2-npml) {
	    g2[j] = pow(c2*((j-(pad2-npml-1))*d2/eta2),2.);
	} else {
	    g2[j] = 0.;
	}
    }

    s2 = (double complex*) sf_alloc(pad2,sizeof(double complex));
    for (j=0; j < pad2; j++) {
	s2[j] = 1./(1.+I*g2[j]);
    }
    
    /* extend model */
    pad = (double**) sf_alloc(pad2,sizeof(double*));
    pad[0] = (double*) sf_alloc(pad1*pad2,sizeof(double));
    for (j=1; j < pad2; j++) {
	pad[j] = pad[0]+j*pad1;
    }

    for (j=0; j < npml; j++) {
	for (i=0; i < npml; i++) {
	    pad[j][i] = v[0][0];
	}
	for (i=npml; i < pad1-npml; i++) {
	    pad[j][i] = v[0][i-npml];
	}
	for (i=pad1-npml; i < pad1; i++) {
	    pad[j][i] = v[0][n1-1];
	}
    }
    for (j=npml; j < pad2-npml; j++) {
	for (i=0; i < npml; i++) {
	    pad[j][i] = v[j-npml][0];
	}
	for (i=npml; i < pad1-npml; i++) {
	    pad[j][i] = v[j-npml][i-npml];
	}
	for (i=pad1-npml; i < pad1; i++) {
	    pad[j][i] = v[j-npml][n1-1];
	}
    }
    for (j=pad2-npml; j < pad2; j++) {
	for (i=0; i < npml; i++) {
	    pad[j][i] = v[n2-1][0];
	}
	for (i=npml; i < pad1-npml; i++) {
	    pad[j][i] = v[n2-1][i-npml];
	}
	for (i=pad1-npml; i < pad1; i++) {
	    pad[j][i] = v[n2-1][n1-1];
	}
    }
    
    /* assemble matrix in triplet form */
    n_row = n_col = (pad1-2)*(pad2-2);
    nz = 5*(pad1-2)*(pad2-2)-2*(pad1-4)-2*(pad2-4)-8;
    
    Ti = sf_intalloc(nz);
    Tj = sf_intalloc(nz);
    Tx = (double*) sf_alloc(nz,sizeof(double));
    Tz = (double*) sf_alloc(nz,sizeof(double));

    count = 0;
    for (j=1; j < pad2-1; j++) {
	for (i=1; i < pad1-1; i++) {
	    index = (j-1)*(pad1-2)+(i-1);
	    
	    cent = 0.+I*0.;

	    /* left */
	    neib = (s1[i]/s2[j]+s1[i-1]/s2[j])/(2.*d1*d1);
	    cent += -neib;

	    if (i != 1) {
		Ti[count] = index;
		Tj[count] = index-1;
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* right */
	    neib = (s1[i]/s2[j]+s1[i+1]/s2[j])/(2.*d1*d1);
	    cent += -neib;

	    if (i != pad1-2) {
		Ti[count] = index;
		Tj[count] = index+1;
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* down */
	    neib = (s2[j]/s1[i]+s2[j-1]/s1[i])/(2.*d2*d2);
	    cent += -neib;

	    if (j != 1) {
		Ti[count] = index;
		Tj[count] = index-(pad1-2);
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* up */
	    neib = (s2[j]/s1[i]+s2[j+1]/s1[i])/(2.*d2*d2);
	    cent += -neib;

	    if (j != pad2-2) {
		Ti[count] = index;
		Tj[count] = index+(pad1-2);
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* center */
	    cent += pow(omega/pad[j][i],2.)/(s1[i]*s2[j]);
	    
	    Ti[count] = index;
	    Tj[count] = index;
	    Tx[count] = creal(cent);
	    Tz[count] = cimag(cent);
	    
	    count++;
	}
    }

    /* convert triplet to compressed-column form */
    Ap = sf_intalloc(n_col+1);
    Ai = sf_intalloc(nz);
    Map = sf_intalloc(nz);

    Ax = (double*) sf_alloc(nz,sizeof(double));
    Az = (double*) sf_alloc(nz,sizeof(double));

    status = umfpack_zi_triplet_to_col (n_row, n_col, nz, 
					Ti, Tj, Tx, Tz, 
					Ap, Ai, Ax, Az, Map);

    /* LU factorization */
    status = umfpack_zi_symbolic (n_row, n_col, 
				  Ap, Ai, Ax, Az, 
				  &Symbolic, Control, Info);

    status = umfpack_zi_numeric (Ap, Ai, Ax, Az, 
				 Symbolic, &Numeric, 
				 Control, Info);

    /* read source */
    if (NULL == sf_getstring("source"))
	sf_error("Need source=");
    source = sf_input("source");

    f = sf_floatalloc2(n1,n2);
    sf_floatread(f[0],n1*n2,source);
    sf_fileclose(source);

    Bx = (double*) sf_alloc(n_col,sizeof(double));
    Bz = (double*) sf_alloc(n_col,sizeof(double));

    for (j=1; j < pad2-1; j++) {
	for (i=1; i < pad1-1; i++) {
	    if (i < npml || i >= pad1-npml || 
		j < npml || j >= pad2-npml) {
		Bx[(j-1)*(pad1-2)+(i-1)] = 0.;
	    } else {
		Bx[(j-1)*(pad1-2)+(i-1)] = f[j-npml][i-npml];
	    }

	    Bz[(j-1)*(pad1-2)+(i-1)] = 0.;
	}
    }    

    /* solve linear system */
    Xx = (double*) sf_alloc(n_col,sizeof(double));
    Xz = (double*) sf_alloc(n_col,sizeof(double));

    status = umfpack_zi_solve (UMFPACK_A, 
			       Ap, Ai, Ax, Az, 
			       Xx, Xz, Bx, Bz, 
			       Numeric, Control, Info);

    /* write output */
    sf_settype(out,SF_COMPLEX);

    utemp = sf_complexalloc(n1);

    for (j=npml; j < pad2-npml; j++) {
	for (i=npml; i < pad1-npml; i++) {
	    utemp[i-npml] = sf_cmplx((float) Xx[(j-1)*(pad1-2)+(i-1)], 
				     (float) Xz[(j-1)*(pad1-2)+(i-1)]);
	}

	sf_complexwrite(utemp,n1,out);
    }

    exit(0);
}
