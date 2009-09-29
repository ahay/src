/* Local coherency and dip based on trace-by-trace complex statistical correlation. */
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

#include <math.h>
#include <float.h>
#include <rsf.h>

#define SIGN(a) (a > 0 ? 1 : (a < 0 ? -1 : 0))

#include "cplxstat.h"

int main(int argc, char* argv[])
{
    int i,j,k,iw;
    int nt,nx,nw,ntw;
    int nti,ntf,ntt;

    float dt,dx,ot,ox,c1,p1,c2,p2;

    float ***cc,**data,*t1,*t2,*t3; 

    sf_axis caxis,taxis,xaxis;
    sf_file in,out;

    sf_init (argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    /* input file parameters */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nt))   sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1= in input");

    if (!sf_histint(in,"n2",&nx))   sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2= in input");

    if (!sf_getint("nw",&nw)) sf_error("No nw= in input");
    /* half time-window size */

    /* time-window length */
    ntw = 2*nw + 1; 

    nti = nw;
    ntf = nt - nw;
    ntt = nt - 2*nw;

    /* output file parameters */
    taxis = sf_maxa(ntt,ot + nw*dt,dt);
    sf_oaxa(out,taxis,1);

    xaxis = sf_maxa(nx-2,ox + dx,dx);
    sf_oaxa(out,xaxis,2);

    caxis = sf_maxa(2,0,1);
    sf_oaxa(out,caxis,3);

    /* memory allocations */
    data = sf_floatalloc2(nt,nx);

    cc = sf_floatalloc3(ntt,nx-2,2);

    t1 = sf_floatalloc(ntw);
    t2 = sf_floatalloc(ntw);
    t3 = sf_floatalloc(ntw);

    /* input in degree */
    sf_floatread(data[0],nt*nx,in);


    /* vertical-window correlation : dip */
    for (i = 1; i < nx-1; i++) { 

	for (j = nti; j < ntf; j++) {

	    for (k = 0; k < ntw; k++) {

		iw = j - nw + k;

		t1[k] = data[i-1][iw];
		t2[k] = data[i][iw];
		t3[k] = data[i+1][iw];

	    }
       
	    c1 = 0.0; p1 = 0.0;
	    c2 = 0.0; p2 = 0.0;

            /* cross-correlation modulus and phase 1-2 */
	    circ_corr(t1,t2,ntw,&c1,&p1);

            /* cross-correlation modulus and phase 2-3 */
	    circ_corr(t2,t3,ntw,&c2,&p2);

	    cc[0][i-1][j] = (p1+p2)*90.0/SF_PI;       
	    cc[1][i-1][j] = 0.5*(c1+c2);

	}

    }

    /* horizontal-window correlation : frequency */


    /* Mixed for local attributes */

    sf_floatwrite(cc[0][0],ntt*(nx-2)*2,out);

    exit (0);

}
