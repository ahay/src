/* Plane-wave spray in 3-D. */
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

#include "predict.h"
#include "update.h"

int main (int argc, char *argv[])
{
    bool verb, up2, up3;
    unsigned char update;
    int n1,n2,n3, i1,i2,i3, ns2, ns3, ip, np2, np3, n23;
    int order, np, i4, n4, k2, k3, j2, j3, i, jp, j;
    float eps, ***u, **p1, **p2, **cost, *trace, *q2=NULL, *q3=NULL;
    sf_file inp, out, dip;

    sf_init(argc,argv);
    inp = sf_input("in");
    dip = sf_input("dip");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(inp,"n3",&n3)) sf_error("No n3= in input");
    n23 = n2*n3;
    n4 = sf_leftsize(inp,3);

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */
    
    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    if (!sf_getint("ns2",&ns2)) sf_error("Need ns1=");
    if (!sf_getint("ns3",&ns3)) sf_error("Need ns2=");
    /* spray radius */
    np2 = 2*ns2+1;
    np3 = 2*ns3+1;
    np = np2*np3;

    sf_putint(out,"n2",np);
    sf_shiftdim(inp, out, 2);

    cost = sf_floatalloc2(np2,np3);
    for (i3=0; i3 < np3; i3++) {
	for (i2=0; i2 < np2; i2++) {
	    cost[i3][i2] = hypotf(i2-ns2,i3-ns3);
	}
    }

    predict_init (n1, n2, eps*eps, order, 1, true);
    update_init(np2,np3,*cost);

    u = sf_floatalloc3(n1,np,n23);
    for (i3=0; i3 < n23; i3++) {
	for (ip=0; ip < np; ip++) {
	    for (i1=0; i1 < n1; i1++) {
		u[i3][ip][i1] = 0.;
	    }
	}
    }

    p1 = sf_floatalloc2(n1,n23);
    p2 = sf_floatalloc2(n1,n23);

    for (i=0; i < n23; i++) { 
	sf_floatread(p1[i],n1,dip);
    }

    for (i=0; i < n23; i++) { 
	sf_floatread(p2[i],n1,dip);
    }

    for (i4=0; i4 < n4; i4++) {
	for (i=0; i < n23; i++) { 
	    sf_floatread(u[i][ns3*np2+ns2],n1,inp);
	    
	    i2 = i%n2;
	    i3 = i/n2;

	    for (ip=0; ip < np; ip++) {
		update = get_update(ip,&up2,&up3,&jp);
		
		/* from jp to j */
		k2 = jp%np2;
		k3 = jp/np2;
		
		j2 = i2+k2-ns2;
		j3 = i3+k3-ns3;

		if (j2 < 0 || j2 >= n2 || 
		    j3 < 0 || j3 >= n3) continue;

		j = j2+j3*n2;
		trace = u[j][jp];

		if (update & 1) {		
		    if (up2) {
			if (j2==0) continue;
			j2 = j-1;
			q2 = p1[j2];
			k2 = jp-1;
		    } else {
			if (j2==n2-1) continue;
			j2 = j+1;
			q2 = p1[j];
			k2 = jp+1;
		    }
		}
		if (update & 2) {
		    if (up3) {
			if (j3==0) continue;
			j3 = j-n2;
			q3 = p2[j3];
			k3 = jp-np2;
		    } else {
			if (j3==n3-1) continue;
			j3 = j+n2;
			q3 = p2[j];
			k3 = jp+np2;
		    }
		}

		switch(update) {
		    case 0:			
			break;
		    case 1:
			predict1_step(up2,u[j2][k2],q2,trace);
			break;
		    case 2:
			predict1_step(up3,u[j3][k3],q3,trace);
			break;
		    case 3:
			predict2_step(up2,up3,u[j2][k2],u[j3][k3],
				      q2,q3,trace);
			break;
		}
	    }
	}

	for (i=0; i < n23; i++) {
	    for (ip=0; ip < np; ip++) {
		sf_floatwrite(u[i][ip],n1,out);
	    }
	}
    }

    exit (0);
}
