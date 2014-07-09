/* 5D sparsity-promoting LS interpolation with real acquisition geometry
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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
#include <complex.h>

#include "int4.h"
#include "fftn.h"

int main(int argc, char* argv[])
{
    bool verb;
    int n1, n2, n3, n4, n5, nm, nd, i, iter, niter, nouter, n[5];
    float o1, o2, o3, o4, o5, d1, d2, d3, d4, d5, eps;
    float *dat, *ax11, *ax12, *ax21, *ax22, *w;
    sf_complex *dobs, *mm, *dd;
    sf_file Fin, Fout, Fax11, Fax12, Fax21, Fax22;/* I/O files */ 

    sf_init(argc,argv);		/* Madagascar initialization */
    Fin=sf_input("in");		/* read the data to be interpolated */
    Fax11=sf_input("axis11");	/* 1st coordinate axis of component 1 */
    Fax12=sf_input("axis12");	/* 2nd coordinate axis of component 1 */
    Fax21=sf_input("axis21");	/* 1st coordinate axis of component 2 */
    Fax22=sf_input("axis22");	/* 1st coordinate axis of component 2 */
    Fout=sf_output("out"); 	/* output the reconstructed data */

    /* Read the data size, origin and the interval */
    if (!sf_histint(Fin,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(Fin,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(Fin,"n3",&n3)) sf_error("No n3= in input");
    if (!sf_histint(Fin,"n4",&n4)) sf_error("No n4= in input");
    if (!sf_histint(Fin,"n5",&n5)) sf_error("No n5= in input");
    if (!sf_histfloat(Fin,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(Fin,"o2",&o2)) sf_error("No o2= in input");
    if (!sf_histfloat(Fin,"o3",&o3)) sf_error("No o3= in input");
    if (!sf_histfloat(Fin,"o4",&o4)) sf_error("No o4= in input");
    if (!sf_histfloat(Fin,"o5",&o5)) sf_error("No o5= in input");
    if (!sf_histfloat(Fin,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(Fin,"d2",&d2)) sf_error("No d2= in input");
    if (!sf_histfloat(Fin,"d3",&d3)) sf_error("No d3= in input");
    if (!sf_histfloat(Fin,"d4",&d4)) sf_error("No d4= in input");
    if (!sf_histfloat(Fin,"d5",&d5)) sf_error("No d5= in input");

    if (!sf_getbool("verb",&verb))    	verb=false;/* verbosity */
    if (!sf_getint("niter",&niter)) 	niter=100; /* inner iterations */
    if (!sf_getint("nouter",&nouter)) 	nouter=5;  /* outer iterations */
    if (!sf_getfloat("eps",&eps)) 	eps=1.e-2; /* regularization parameter */
    n[0]=n1; 
    if (!sf_getint("sn11", &n[1])) sf_error("No n1 for interpolated component 1");
    if (!sf_getint("sn12", &n[2])) sf_error("No n2 for interpolated component 1");
    if (!sf_getint("sn21", &n[3])) sf_error("No n1 for interpolated component 2");
    if (!sf_getint("sn22", &n[4])) sf_error("No n1 for interpolated component 2");

    nd=n1*n2*n3*n4*n5;
    nm=n[0]*n[1]*n[2]*n[3]*n[4];

    /* allocate data and mask arrays */
    dat=sf_floatalloc(nd);
    ax11=sf_floatalloc(n2);
    ax12=sf_floatalloc(n3);
    ax21=sf_floatalloc(n4);
    ax22=sf_floatalloc(n5);
    w=sf_floatalloc(nm);
    dobs=sf_complexalloc(nd);
    mm=sf_complexalloc(nm);
    dd=sf_complexalloc(nm);

    sf_floatread(dat, nd, Fin);
    sf_floatread(ax11, n2, Fax11);
    sf_floatread(ax12, n3, Fax12);
    sf_floatread(ax21, n4, Fax21);
    sf_floatread(ax22, n5, Fax22);
    for(i=0; i<nd; i++) dobs[i]=dat[i];
    int4_init(o2, o3, o4, o5, d2, d3, d4, d5, n2, n3, n4, n5, n1, n[1]*n[2]*n[3]*n[4], ax11, ax12, ax21, ax22);
    fftn_init(5, n);
    for(i=0; i<nm; i++) w[i]=1.0;

    for(iter=0; iter<nouter; iter++){
	sf_csolver_prec(int4_lop, sf_ccgstep, fftn_lop, nm, nm, nd, dd, dobs, niter, eps, "mwt", w, "xp", mm, "verb", verb,"end");
	for(i=0; i<nm; i++) w[i]=cabsf(mm[i]); 
    }

    sf_ccgstep_close();
    int4_close();
    fftn_close();

    for(i=0; i<nm; i++) w[i]=crealf(dd[i]);
    sf_floatwrite(w, nm, Fout); /* output reconstructed seismograms */

    free(dat);
    free(ax11);
    free(ax12);
    free(ax21);
    free(ax22);
    free(w);
    free(dobs);
    free(mm);
    free(dd);

    exit(0);
}
