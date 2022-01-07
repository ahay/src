/* Select seismic data cross well position. */
/*
  Copyright (C) 2012 Jilin University
  
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

int main (int argc, char* argv[]) 
{
    int n1, n2, n3, nw, wtemp, ndcs;
    int i, i1, i2, i3, j=0, j1, j2, lr, ud, n23;
    int o2, d2, o3, d3, *nwx, *nwy;
    
    float *sdata, *wx, *wy, *swdata, **cost, **path;
    sf_file in, out, well=NULL, spath=NULL;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");

    if (NULL != sf_getstring("path")) {
	spath = sf_output("path");
    } else {
	spath = NULL;
    }

    if (NULL != sf_getstring ("well")) {
	well = sf_input("well");
	if (!sf_histint(well,"n2",&nw)) sf_error("no n2= in well");
	if (!sf_histint(well,"n1",&wtemp)) sf_error("no n1= in well");
	if (2 != wtemp) sf_error("n1 of well should be 2");
    } else {
	well = NULL;
	if (!sf_getint("nw",&nw)) nw=0;
	/* number of well points */
	if (0 == nw) sf_error("Need well location");
    }
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    n23 = n2*n3;
    if (!sf_histint(in,"o2",&o2)) sf_error("No o2= in input");
    if (!sf_histint(in,"o3",&o3)) sf_error("No o3= in input");
    if (!sf_histint(in,"d2",&d2)) sf_error("No d2= in input");
    if (!sf_histint(in,"d3",&d3)) sf_error("No d3= in input");
    
    sdata = sf_floatalloc(n1*n2*n3);
    nwx = sf_intalloc(nw);
    nwy = sf_intalloc(nw);
    wx = sf_floatalloc(nw);
    wy = sf_floatalloc(nw);
    swdata = sf_floatalloc(n1*n2*n3);
    path = sf_floatalloc2(2,n23);

    sf_floatread(sdata,n1*n2*n3,in);

    if (NULL != well) {
	for (i=0; i< nw; i++) {
	    sf_floatread(wx,1,well);
	    sf_floatread(wy,1,well);
	}
    } else {
	if (!sf_getfloats("wx",wx,nw)) sf_error("Need wx=");
	/* well x coordinates */
	if (!sf_getfloats("wy",wy,nw)) sf_error("Need wy=");
	/* well y coordinates */
    } 

    for (i=0; i < nw; i++) {
	if (wx[i] < o2 || wx[i] > o2+(n2-1)*d2 ||
	    wy[i] < o3 || wy[i] > o3+(n3-1)*d3) {
	    sf_error("Well location is out of range of seismic data");
	} else {
	    nwx[i] = (int)(wx[i]-o2)/d2;
	    nwy[i] = (int)(wy[i]-o3)/d3;
	}
    }

    cost = sf_floatalloc2(n2,n3);
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    cost[i3][i2] = 1.;
	}
    }

    dijkstra_init(n2,n3,cost,cost);
/*  dijkstra_source(nwx[0],nwy[0]);
    while (dijskstra_step(&j1,&j2,&ud,&lr)) {;} 
*/
    ndcs = 0;
    for (i=0; i < nw-1; i++) {
	dijkstra_source(nwx[i+1],nwy[i+1]);
	while (dijskstra_step(&j1,&j2,&ud,&lr)) {;} 
	dijkstra_start(nwx[i],nwy[i]);

	j1 = nwx[i];
	j2 = nwy[i];
	path[ndcs][0]=j1*d2+o2;
	path[ndcs][1]=j2*d3+o3;
	for (i1=0; i1 < n1; i1++) {
	    swdata[ndcs*n1+i1] = sdata[j2*n1*n2+j1*n1+i1];
	} 
	
	for (j=1; dijkstra_next(&ud,&lr); j++) {
	    if (ud > 0) {
		j1 -= ud;
	    } else {
		j1 -= ud;
	    }
	    if (lr > 0) {
		j2 -= lr;
	    } else {
		j2 -= lr;
	    }
	    ndcs ++;
	    path[ndcs][0]=j1*d2+o2;
	    path[ndcs][1]=j2*d3+o3;
	    for (i1=0; i1 < n1; i1++) {
		swdata[ndcs*n1+i1] = sdata[j2*n1*n2+j1*n1+i1];
	    } 
	}
    }

    if (NULL != spath) {
	sf_putint(spath,"n1",2);
	sf_putint(spath,"n2",ndcs+1);
	sf_putint(spath,"n3",1);
	sf_floatwrite(path[0],2*(ndcs+1),spath);
	sf_fileclose(spath);
    }
    
    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",ndcs+1);
    sf_putint(out,"n3",1);
    sf_floatwrite(swdata,n1*(ndcs+1),out);

    exit (0);
}

/* 	$Id$	 */
