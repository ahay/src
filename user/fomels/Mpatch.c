/* Testing patching. 

Takes: w1=n1 w2=n2 w3=n3 p1=1 p2=1 p3=1

w1,w2,w3 is window size, p1,p2,p3 is number of patches. 
*/
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
#include <stdio.h>
#include <unistd.h>

#include <rsf.h>

#include "dip3.h"
#include "mask6.h"
#include "ocpatch.h"
#include "oc.h"

int main (int argc, char *argv[])
{
    int w123,p123,n123, j, ip,iw;
    int n[3], nw[3], w[3]; 
    off_t nall;
    float *u, *tent, *tmp;
    char key[3], *dipname, *wallname;
    sf_file in, out;
    FILE *wall, *dip;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    if (!sf_histint(in,"n1",&n[0])) sf_error("Need n1= in input");
    if (!sf_histint(in,"n2",&n[1])) n[1]=1;
    if (!sf_histint(in,"n3",&n[2])) n[2]=1;
    n123 = n[0]*n[1]*n[2];

    for (j=0; j <3; j++) {
	snprintf(key,3,"w%d",j+1);
	if (!sf_getint(key,w+j) || w[j] > n[j]) w[j] = n[j];
	
	snprintf(key,3,"p%d",j+1);
	if (!sf_getint(key,nw+j)) {
	    if (n[j] > w[j]) {
		nw[j] = 1 + 1.5*n[j]/w[j]; /* 50% overlap */
	    } else {
		nw[j] = 1;
	    }
	}
    }
    w123 = w[0]*w[1]*w[2];
    p123 = nw[0]*nw[1]*nw[2];

    u = sf_floatalloc(w123);

    sf_warning("window: %d %d %d",w[0],w[1],w[2]);
    sf_warning("patches: %d %d %d",nw[0],nw[1],nw[2]);

    nall = n123*sizeof(float);
    
    sf_unpipe(in,nall);
    wall = sf_tempfile(&wallname,"w+b");
    dip = sf_tempfile(&dipname,"w+b");

    tent = sf_floatalloc(w123);
    tmp = sf_floatalloc(w123);
    
    sf_tent2 (3, w, tent);
    
    ocpatch_init(3,w123,p123,nw,n,w);
    oc_zero(nall,wall);
    oc_zero(nall,dip);
    
    /* loop over patches */
    for (ip=0; ip < p123; ip++) {
	sf_warning("patch %d of %d",ip+1,p123);

	/* read data */
	ocpatch_flop (ip,false,in,u);
		
	/* write weight */
	ocpatch_lop (ip,false,wall,tmp);
	for (iw=0; iw < w123; iw++) {
	    tmp[iw] += tent[iw];
	}
	ocpatch_lop (ip,true,wall,tmp);
		
	/* write dip */
	ocpatch_lop (ip,false,dip,tmp);
	for (iw=0; iw < w123; iw++) {
	    tmp[iw] += tent[iw]*u[iw];
	}
	ocpatch_lop (ip,true,dip,tmp);
    }
	    
    oc_divide(nall,dip,wall,out);
    unlink(dipname);
    unlink(wallname);
    
    exit (0);
}

/* 	$Id: Mdip.c 1071 2005-03-20 20:08:23Z fomels $	 */
