/* Make topography mask / receiver list / record list for first-arrival traveltime tomography */
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

#include "fastmarch.h"

int main(int argc, char* argv[])
{
    bool velocity, plane[3], skip;
    int n[SF_MAX_DIM], nt, dim, *flag, order, i, j, k, is, ns, nrecv, *mask, *recv;
    int temp[2], offset[2], left[2], right[2], count;
    float **source, d[SF_MAX_DIM], o[SF_MAX_DIM], air, *s, *t, *reco;
    char key[4];
    sf_file in, shot, out, record, topo;
    
    sf_init(argc,argv);
    in     = sf_input("in");
    out    = sf_output("out");
    record = sf_output("record");
    
    /* read model dimensions */
    dim = sf_filedims(in,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) sf_error("No %s= in input",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
	nt *= n[i]; plane[i] = false;
    }
    if (dim < 3) {
	n[2] = 1; o[2] = o[1]; d[2] = d[1]; plane[2] = false;
    }
    
    /* read model */
    s = sf_floatalloc(nt);
    sf_floatread(s,nt,in);

    /* read shot list */
    if (NULL == sf_getstring("shot"))
	sf_error("Need shot list shot=");
    shot = sf_input("shot");

    if (!sf_histint(shot,"n2",&ns)) ns=1;

    source = sf_floatalloc2(3,ns);
    sf_floatread(source[0],3*ns,shot);
    sf_fileclose(shot);

    if (!sf_getint("offset1",&offset[1])) offset[1]=0;
    /* receiver offset inline (on each side) */

    if (!sf_getint("offset2",&offset[2])) offset[2]=0;
    /* receiver offset crossline (on each side) */

    nrecv = ((2*offset[1]+1)*(2*offset[2]+1)>n[1]*n[2])?n[1]*n[2]
	:(2*offset[1]+1)*(2*offset[2]+1)-1;

    /* allocate memory for output */
    recv = sf_intalloc(nrecv);
    reco = sf_floatalloc(nrecv);

    if (NULL != sf_getstring("topo")) {
	topo = sf_output("topo");
	mask = sf_intalloc(nt);
    } else {
	topo = NULL;
	mask = NULL;
    }

    /* write output header */
    sf_putint(out,"n1",nrecv);
    sf_putint(out,"n2",ns);
    sf_putint(out,"n3",1);
    sf_settype(out,SF_INT);

    sf_putint(record,"n1",nrecv);
    sf_putint(record,"n2",ns);
    sf_putint(record,"n3",1);
    
    if (topo != NULL) sf_settype(topo,SF_INT);

    if (!sf_getfloat("air",&air)) air=0.5;
    /* air velocity for thresholding topography mask */

    if (!sf_getint("order",&order)) order=2;
    /* fast marching accuracy order */

    if (!sf_getbool("velocity",&velocity)) velocity=true;
    /* if y, the input is velocity; n, slowness squared */

    /* convert to slowness squared */
    if (velocity) {
	for (i=0; i < nt; i++) {
	    s[i] = 1./s[i]*1./s[i];
	}

	air = 1./air*1./air;
    }
    
    /* allocate temporary memory */
    t     = sf_floatalloc(nt);
    flag  = sf_intalloc(nt);
    
    /* make topography mask */
    if (topo != NULL) {
	for (k=0; k < n[2]; k++) {
	    for (j=0; j < n[1]; j++) {

		skip = true;
		for (i=0; i < n[0]; i++) {

		    if (s[k*n[1]*n[0]+j*n[0]+i] < air && skip) {
			mask[k*n[1]*n[0]+j*n[0]+i] = 1;
			skip = false;
		    } else {
			mask[k*n[1]*n[0]+j*n[0]+i] = 0;
		    }
		}
	    }
	}

	sf_intwrite(mask,nt,topo);
    }
    
    /* make receiver list and record list */
    fastmarch_init(n[2],n[1],n[0]);

    for (is=0; is < ns; is++) {
	fastmarch(t,s,flag,plane,
		  n[2],n[1],n[0],o[2],o[1],o[0],d[2],d[1],d[0],
		  source[is][2],source[is][1],source[is][0],1,1,1,order);

	temp[1]  = (source[is][1]-o[1])/d[1]+0.5;
	left[1]  = temp[1]-offset[1];
	right[1] = temp[1]+offset[1];
	
	temp[2]  = (source[is][2]-o[2])/d[2]+0.5;
	left[2]  = temp[2]-offset[2];
	right[2] = temp[2]+offset[2];
	
	count = 0;
	for (k=0; k < n[2]; k++) {
		for (j=0; j < n[1]; j++) {

		    if (left[1] <= j && j <= right[1] 
			&& left[2] <= k && k <= right[2]) {

			if (j == temp[1] && k == temp[2])
			    continue;
			
			skip = false;
			for (i=0; i < n[0]; i++) {
			    if (skip) break;
			    
			    if (s[k*n[1]*n[0]+j*n[0]+i] < air) {
				recv[count] = k*n[1]*n[0]+j*n[0]+i;
				reco[count] = t[k*n[1]*n[0]+j*n[0]+i];
				skip = true;
				count++;
			    }
			}
		    }
		}	    
	}
	
	/* negative flag for void space */
	if (count < nrecv) {
	    for (i=count; i < nrecv; i++) {
		recv[i] = -1;
		reco[i] = -1.0;
	    }
	}

	sf_intwrite(recv,nrecv,out);
	sf_floatwrite(reco,nrecv,record);
    }
    
    exit(0);
}
