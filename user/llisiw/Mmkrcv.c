/* Make topography mask / receiver list / record list for first-arrival traveltime tomography */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mkrcv.h"

int main(int argc, char* argv[])
{
    bool velocity, fix, plane;
    int n[SF_MAX_DIM], nt, dim, order, ip, np, i, j, k, is, ns, ic, nc, nrecv, *mask, **recv;
    int temp[2], offset[2], left[2], right[2], count;
    float **source, d[SF_MAX_DIM], o[SF_MAX_DIM], air, p, p0, dp, *s, **t, **data, **code;
    char key[4];
    sf_file in, shot, out, reco, topo, time;
    
    sf_init(argc,argv);
    in   = sf_input("in");
    out  = sf_output("out");
    reco = sf_output("reco");
    
    /* read model dimensions */
    dim = sf_filedims(in,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) sf_error("No %s= in input",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
	nt *= n[i];
    }
    if (dim < 3) {
	n[2] = 1; o[2] = o[1]; d[2] = d[1];
    }
    
    /* read model */
    s = sf_floatalloc(nt);
    sf_floatread(s,nt,in);

    if (!sf_getfloat("air",&air)) air=0.5;
    /* air velocity for thresholding topography */

    if (!sf_getbool("velocity",&velocity)) velocity=true;
    /* if y, the input is velocity; n, slowness squared */

    /* convert to slowness squared */
    if (velocity) {
	for (i=0; i < nt; i++) {
	    s[i] = 1./s[i]*1./s[i];
	}
	
	air = 1./air*1./air;
    }
    
    /* read shot list */
    if (NULL == sf_getstring("shot"))
	sf_error("Need shot list shot=");
    shot = sf_input("shot");

    if (!sf_histint(shot,"n2",&ns)) ns=1;

    source = sf_floatalloc2(3,ns);
    sf_floatread(source[0],3*ns,shot);
    sf_fileclose(shot);

    if (!sf_getint("order",&order)) order=2;
    /* fast marching accuracy order */

    if (!sf_getbool("fix",&fix)) fix=false;
    /* if y, fixed-spread; n, moving acquisition */

    if (!sf_getbool("plane",&plane)) plane=false;
    /* if y, plane-wave source; n, point source */

    if (!sf_getint("offset1",&offset[0])) offset[0]=0;
    /* receiver offset inline */

    if (!sf_getint("offset2",&offset[1])) offset[1]=0;
    /* receiver offset crossline */

    if (plane) {
	/* NOTE: plane-wave source only supports fixed-spread acquisition */
	fix = true;

	if (!sf_getint("np",&np)) np=1;
	/* ray-parameter number */
	
	if (!sf_getfloat("p0",&p0)) p0=0.;
	/* ray-parameter start */
	
	if (!sf_getfloat("dp",&dp)) dp=1.;
	/* ray-parameter increment */

	/* coding for plane-wave source */
	code = sf_floatalloc2(ns,np);
	
	for (ip=0; ip < np; ip++) {
	    p = p0+ip*dp;
	    
	    /* NOTE: 2D only! */
	    for (is=0; is < ns; is++) {
		code[ip][is] = (source[is][1]-source[p>=0.?0:ns-1][1])*p;
	    }
	}

	nc = np;
    } else {
	code = NULL;

	nc = ns;
    }

    /* determine maximum allowed receiver length */
    if (fix)
	/* fixed-spread: offset defines gap from domain border */
	nrecv = (n[1]-2*offset[0])*(n[2]-2*offset[1]);
    else
	/* moving acquisition: offset defines distance from source */
	nrecv = ((2*offset[0]+1)*(2*offset[1]+1)>n[1]*n[2])?n[1]*n[2]
	    :(2*offset[0]+1)*(2*offset[1]+1);

    /* allocate memory for output */
    recv = sf_intalloc2(nrecv,nc);
    data = sf_floatalloc2(nrecv,nc);

    if (NULL != sf_getstring("topo")) {
	topo = sf_output("topo");
	mask = sf_intalloc(nt);
    } else {
	topo = NULL;
	mask = NULL;
    }

    if (NULL != sf_getstring("time")) {
	time = sf_output("time");
    } else {
	time = NULL;
    }

    /* write output header */
    sf_putint(out,"n1",nrecv);
    sf_putint(out,"n2",nc);
    if (plane) {
	sf_putfloat(out,"o2",p0);
	sf_putfloat(out,"d2",dp);
    }
    sf_putint(out,"n3",1);
    sf_settype(out,SF_INT);

    sf_putint(reco,"n1",nrecv);
    sf_putint(reco,"n2",nc);
    if (plane) {
	sf_putfloat(reco,"o2",p0);
	sf_putfloat(reco,"d2",dp);
    }
    sf_putint(reco,"n3",1);
    
    if (topo != NULL) sf_settype(topo,SF_INT);
    
    if (time != NULL) sf_putint(time,"n4",nc);

    /* make topography mask */
    if (topo != NULL) {
	for (i=0; i < nt; i++) {
	    if (s[i] < air)
		mask[i] = 1;
	    else
		mask[i] = 0;
	}
	
	sf_intwrite(mask,nt,topo);
    }
    
    /* initialize */
    fastmarch_init(n,o,d,order);

    /* set background slowness */
    fastmarch_set(s);

    /* allocate temporary memory */
    t = sf_floatalloc2(nt,nc);

#ifdef _OPENMP
#pragma omp parallel for private(count,temp,left,right,k,j,i)
#endif
    for (ic=0; ic < nc; ic++) {

	if (plane)
	    fastmarch_plane(t[ic],ns,source,code[ic]);
	else
	    fastmarch_point(t[ic],source[ic]);
	
	if (fix) {
	    
	    count = 0;
	    for (k=offset[1]; k < n[2]-offset[1]; k++) {
		for (j=offset[0]; j < n[1]-offset[0]; j++) {
		    
		    if (!plane) {
			if ((j == (source[ic][1]-o[1])/d[1]+0.5) &&
			    (k == (source[ic][2]-o[2])/d[2]+0.5))
			    continue;
		    }
		    
		    for (i=0; i < n[0]; i++) {
			if (s[k*n[1]*n[0]+j*n[0]+i] < air) {
			    recv[ic][count] = k*n[1]*n[0]+j*n[0]+i;
			    data[ic][count] = t[ic][k*n[1]*n[0]+j*n[0]+i];
			    count++;
			    break;
			}
		    }
		}
	    }
	    
	} else {
	    
	    temp[0]  = (source[ic][1]-o[1])/d[1]+0.5;
	    left[0]  = temp[0]-offset[0];
	    right[0] = temp[0]+offset[0];
	    
	    temp[1]  = (source[ic][2]-o[2])/d[2]+0.5;
	    left[1]  = temp[1]-offset[1];
	    right[1] = temp[1]+offset[1];
	    
	    count = 0;
	    for (k=0; k < n[2]; k++) {
		for (j=0; j < n[1]; j++) {
		    
		    if (left[0] <= j && j <= right[0] 
			&& left[1] <= k && k <= right[1]) {
			
			for (i=0; i < n[0]; i++) {
			    if (s[k*n[1]*n[0]+j*n[0]+i] < air) {
				recv[ic][count] = k*n[1]*n[0]+j*n[0]+i;
				data[ic][count] = t[ic][k*n[1]*n[0]+j*n[0]+i];
				count++;
				break;
			    }
			}
		    }
		}	    
	    }
	    
	}
	
	/* negative flag for void space */
	if (count < nrecv) {
	    for (i=count; i < nrecv; i++) {
		recv[ic][i] = -1;
		data[ic][i] = 0.;
	    }
	}
	
    }
    
    sf_intwrite(recv[0],nrecv*nc,out);
    sf_floatwrite(data[0],nrecv*nc,reco);
    if (NULL != time) sf_floatwrite(t[0],nt*nc,time);

    exit(0);
}
