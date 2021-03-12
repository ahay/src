/* L1 1D matched filter */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
#include <rsfgee.h>

int main(int argc, char* argv[])
{
    bool verb;
	int nd,n1,nb, nm, nb2, niter, liter, iter, id,ib;
    float *d, *dd, *m, *o, *n, *b, *btmp,*r, eb,en,perc;
    sf_file inp=NULL, mult=NULL,filt=NULL, out=NULL;

    sf_init(argc,argv);
    inp = sf_input("in");
    filt = sf_output("filt");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");

	if (!sf_getint("nb",&nb)) nb=3;
	/* matched-filter order */

	if (NULL != sf_getstring("mult")) mult=sf_input("mult");
	else sf_error("No multiple file specified");

    if (!sf_histint(mult,"n1",&nm) || n1 != nm)
	sf_error("Need nm=%d is not n1=%d",nm,n1);

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    nd=n1+nb-1;
	nb2 = (nb-1)/2;

    d = sf_floatalloc(nd);
    dd = sf_floatalloc(nd);
    n = sf_floatalloc(nd);
    r = sf_floatalloc(nd);

    b = sf_floatalloc(nb);
    btmp = sf_floatalloc(nb);

    /* initialize with zero */
    for (id=0;id<nd;id++) dd[id]=n[id]=r[id]=d[id]=0.0;

    /*multiple vector*/
    m = sf_floatalloc(n1);
    /*output vector*/
    o = sf_floatalloc(nd);

    /* read data and multiple from file */
	sf_floatread(d+nb2,n1,inp);
    sf_floatread(m,n1,mult);

    sf_putint(filt,"n1",nb);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of POCS iterations */

    /* number of linear CG iterations */
    liter=5;

    if (!sf_getfloat("perc",&perc)) perc=90.0;
    /* percentage for sharpening */

    sf_sharpen_init(nd,perc,0.5);
    /* initialization of sharpening regularization*/

    /* initialize model (filter) with zero */
    for (ib=0; ib < nb; ib++) btmp[ib]=b[ib]=0.0;


    tcaf1_init(n1    /* multiple length */,
			   m    /* multiple [n1] */);


	sf_solver(tcaf1_lop,sf_cgstep,nb,nd,b,d,nb,"x0",b,"verb",verb,"end");
    //    L2 solution
    if (niter==1) {
    	sf_warning("\nL2 solution");
    }
    else {  //    L1 solution
    	sf_warning("\nL1 solution");

		for (iter=0; iter < niter; iter++) {
		/* Solve | (d - n) - M * b |_2 */
		/* -------------------------------------- */
			eb=en=0;
			for (id=0; id < nd; id++)
				r[id] = dd[id] = d[id]-n[id]; /* dd = d - n */

			for (ib=0; ib < nb; ib++)
				btmp[ib] = (-1) * b[ib];  /* -b */

			tcaf1_lop(false,true,nb,nd,btmp,r); /* r = dd - M * b;*/

			sf_solver(tcaf1_lop,sf_cgstep,nb,nd,b,dd,liter,"x0",b,"verb",false,"end");

			for (ib=0; ib < nb; ib++) {
				btmp[ib]= -btmp[ib] - b[ib];  /* -db */
				eb+=b[ib]*b[ib];
			}

			for (id=0; id < nd; id++)
					n[id] += r[id];

			tcaf1_lop(false,true,nb,nd,btmp,n); /* n[i] += r - M * db;*/

			for (id=0; id < nd; id++)
				en+= n[id]*n[id];

			/* apply sharpening regularization*/
			sf_sharpen(n);
			sf_weight_apply(nd,n);

		if (verb) sf_warning("iter=%d: b[%d]=%g eb=%g en=%g",iter,nb2,b[nb2],eb,en);

		}
	}
    sf_floatwrite(b,nb,filt);

    tcai1_init(nb    /* filter length */,
				b    /* filter samples [nb] */);

    tcai1_lop(false,false,n1,nd,m,o);
    sf_floatwrite(o+nb2,n1,out);

    free(d);
    free(dd);
    free(n);
    free(r);
    free(b);
    free(btmp);
    free(o);
    free(m);

    exit(0);
}
