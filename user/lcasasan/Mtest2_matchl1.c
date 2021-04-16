/* L1 1D matched filter */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

static int nb;
static int n1;
static float *mmult;

static void fit(bool adj, bool add, int nm, int nd, float *m, float *d)
{

    sf_adjnull(adj, add, nm, nd, m, d);
    sf_copy_lop(adj, true, nd, nd, m, d);


    tcaf1_init(n1    /* multiple length */,
			  mmult    /* multiple [n1] */);

	tcaf1_lop(adj,true,nb,nd,m+nd,d);
}



int main(int argc, char* argv[])
{
    bool verb;
	int nd, nm,nb2,nmult,niter, liter, iter, i1;
    float *d, *m, *o, perc;
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

    if (!sf_histint(mult,"n1",&nmult) || n1 != nmult)
	sf_error("Need nmult=%d is not n1=%d",nmult,n1);

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    nd=n1+nb-1;
	nm=nd+nb;

	nb2 = (nb-1)/2;

    d = sf_floatalloc(nd);
    m = sf_floatalloc(nm);
    o = sf_floatalloc(nd);
    mmult = sf_floatalloc(n1);

	sf_floatread(d+nb2,n1,inp);
    sf_floatread(mmult,n1,mult);

    sf_putint(filt,"n1",nb);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of POCS iterations */

    if (!sf_getint("liter",&liter)) liter=nb;
    /* number of CG iterations */

    if (!sf_getfloat("perc",&perc)) perc=90.0;
    /* percentage for sharpening */

    sf_sharpen_init(nd,perc,0.5);
    /* initialization of sharpening regularization*/

    /* initialize model (data+filter) with zero */
    for (i1=0; i1 < nm; i1++) {
	m[i1]=0.0;
    }

    for (iter=0; iter < niter; iter++) {
    	sf_solver(fit,sf_cgstep,nm,nd,m,d,liter,"x0",m,"verb",verb,"end");
    	sf_sharpen(m);
    	sf_weight_apply(nd,m);
    /* apply sharpening regularization*/
    }

    sf_floatwrite(m+nd,nb,filt);

    tcai1_init(nb    /* filter length */,
			   m+nd    /* filter samples [nb] */);

    tcai1_lop(false,false,n1,nd,mmult,o);
    sf_floatwrite(o+nb2,n1,out);

    free(m);
    free(o);
    free(d);
    free(mmult);

    exit(0);
}
