/* L1 regression 0 ~= d - G * m
 *
 * adapted from sfbil1
 * */

/*
  Copyright (C) 2010 Politecnico di Milano
  
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
#include <rsfpwd.h>

int main(int argc, char* argv[])
{
    bool verb;
	int nd, n1, niter, liter, iter, id, ib, nb;
    float *n, *d,*dd, *r, *b,*btmp, **A, **Atmp, perc;
    double eb = 0.0,en = 0.0;
  //  double ad, bd, aa, bb, a0, b0, da, db, ab, det;
    sf_file inp, reg, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    reg = sf_input("reg");
    out = sf_output("out");

    if (!sf_getbool("verb",&verb)) verb=false;

    if (!sf_histint(inp,"n1",&nd)) sf_error("No n1= in input");
    if (!sf_histint(reg,"n1",&n1) || n1 != nd)
	sf_error("Need n1=%d in reg",nd);
    sf_histint(reg,"n2",&nb);

    sf_putint(out,"n1",nb);
    
    d = sf_floatalloc(nd);
    dd = sf_floatalloc(nd);

    n = sf_floatalloc(nd);
    r = sf_floatalloc(nd);
    Atmp = sf_floatalloc2(nd,nb); /*nb * nd*/
    A = sf_floatalloc2(nb,nd); /*nb * nd*/

    b = sf_floatalloc(nb);
    btmp = sf_floatalloc(nb);

    sf_floatread(d,nd,inp);
    sf_floatread(Atmp[0],nd*nb,reg);

    for (ib=0;ib<nb;ib++)
    	for (id=0;id<nd;id++)
			A[id][ib]=Atmp[ib][id];

    sf_fileclose(reg);

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of POCS iterations */

    if (!sf_getint("Liter",&liter)) liter=10;
    /* number of CG iterations */

    if (!sf_getfloat("perc",&perc)) perc=90.0;
    /* percentage for sharpening */

    sf_sharpen_init(nd,perc,0.5);

    /* initialize model (filter) with zero */
    for (ib=0; ib < nb; ib++) {
	b[ib]=0.0;
    }

    /* initialize with zero */
    for (id=0; id < nd; id++) {
	n[id] = 0.;
    }

    matmult_init(A);

    for (iter=0; iter < niter; iter++) {
	/* Solve |d - alpha * a - beta * b - n|_2 */
	/* -------------------------------------- */
    	for (id=0; id < nd; id++) {
    		 r[id] = dd[id] = d[id]-n[id];
     	}

    	for (ib=0; ib < nb; ib++) {
    		btmp[ib]= (-1) * b[ib];  /* -b */
    	}

    	matmult_lop(false,true,nb,nd,btmp,r); /* r = (d-n)-M*b;*/

    	for (id=0; id < nd; id++) {
        	n[id] += r[id];
        }

    	sf_solver(matmult_lop,sf_cgstep,nb,nd,b,dd,liter,"x0",b,"verb",false,"end");

    	for (ib=0; ib < nb; ib++) {
    		btmp[ib]= -btmp[ib] - b[ib];  /* -db */
    		eb+=b[ib]*b[ib];
    	}

    	matmult_lop(false,true,nb,nd,btmp,n); /* n[i] += r-M*b;*/

       	for (id=0; id < nd; id++) {
            en+= n[id]*n[id];
        }

    /* apply sharpening regularization*/
	/* Threshold n */
	/* ----------- */
	sf_sharpen(n);
	sf_weight_apply(nd,n);

	if(verb) sf_warning("%d %g %g",iter,b[0],b[1]);
    }
    
    sf_floatwrite(b,nb,out);

    exit(0);
}
