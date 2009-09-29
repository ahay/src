/* Dijkstra shortest-path algorithm in 2-D */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include "dijkstra.h"

int main (int argc, char *argv[])
{
    int n1,n2,ref1,ref2,n12,j,j1,j2,ud,lr,i,nf;
    int *fin1,*fin2,**jj;
    float **p,**q;
    char **paths;
    sf_file out, cost, path;

    sf_init(argc,argv);
    cost = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(cost,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(cost,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;

    if (!sf_getint("ref1",&ref1)) ref1=0;
    if (!sf_getint("ref2",&ref2)) ref2=0;
    /* source point */

    p = sf_floatalloc2(n1,n2);
    q = sf_floatalloc2(n1,n2);

    sf_floatread(p[0],n12,cost);
    sf_floatread(q[0],n12,cost);

    sf_putint(out,"n3",1);

    dijkstra_init(n1,n2,p,q);
    dijkstra_source(ref1,ref2);
    while (dijskstra_step(&j1,&j2,&ud,&lr)) {;}
    dijkstra_cost(out);

    if (!sf_getint("nf",&nf)) nf=0;
    /* number of final points */

    if (0 == nf) exit(0);
 
    fin1 = sf_intalloc(nf);
    fin2 = sf_intalloc(nf);
    
    if (!sf_getints("fin1",fin1,nf)) sf_error("Need fin1=");
    if (!sf_getints("fin2",fin2,nf)) sf_error("Need fin2=");
    /* final points */

    jj = sf_intalloc2(2,n12);

    paths = (char**) sf_alloc(nf,sizeof(char*));
    if (!sf_getstrings("paths",paths,nf)) sf_error("Need paths=");
    
    for (i=0; i < nf; i++) { /* loop over final points */
	dijkstra_start(fin1[i],fin2[i]);
	j1 = fin1[i];
	j2 = fin2[i];
	jj[0][0]=j1;
	jj[0][1]=j2;
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
	    jj[j][0]=j1;
	    jj[j][1]=j2;
	}

	path = sf_output(paths[i]);
	sf_settype(path,SF_INT);
	sf_putint(path,"n1",2);
	sf_putint(path,"n2",j);
	sf_putint(path,"n3",1);
    
	sf_intwrite(jj[0],j*2,path);
	sf_fileclose(path);
    }


    exit(0);
}
