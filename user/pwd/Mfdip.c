/* 3D fast dip estimation by plane wave destruction */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>

#include "fdip.h"
#include "mask6.h"

int main(int argc,char**argv)
{
    int dim,n[SF_MAX_DIM], rect[3], n123, n4, nr, ir, j, liter;
    bool verb, **m;
    float *u,*p;
    char key[4];
    sf_file in, out, mask;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    dim = sf_filedims(in,n);
    if (dim < 2) n[1]=1;
    if (dim < 3) n[2]=1;
    n123 = n[0]*n[1]*n[2];
    nr = 1;
    for (j=3; j < dim; j++) 
	nr *= n[j];
	
    if (1 == n[2]) {
	n4=0;
    } else {
	if (!sf_getint("n4",&n4)) n4=2;
	/* what to compute in 3-D. 0: in-line, 1: cross-line, 2: both */ 
	if (n4 > 2) n4=2;
	if (2==n4) {
	    sf_putint(out,"n4",2);
	    for (j=3; j < dim; j++) {
		snprintf(key,4,"n%d",j+2);
		sf_putint(out,key,n[j]);
	    }
	}
    }

    if (!sf_getint("liter",&liter)) liter=20;
    /* number of linear iterations */

    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    /* dip smoothness on 1st axis */
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    /* dip smoothness on 2nd axis */
    if (!sf_getint("rect3",&rect[2])) rect[2]=1;
    /* dip smoothness on 3rd axuis */
    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */


    u=sf_floatalloc(n123);
    if(n4==0 || n4==1)	p=sf_floatalloc(n123);
    else	p=sf_floatalloc(n123*2);

    if (NULL != sf_getstring("mask")) {
	m = sf_boolalloc2(n123,2);
	mask = sf_input("mask");
    } else {
	m = (bool**) sf_alloc(2,sizeof(bool*));
	m[0] = m[1] = NULL;
	mask = NULL;
    }

    fdip_init(n[0], n[1], n[2], rect, liter, verb);
	
    for(ir=0;ir<nr;ir++){
	if(verb)	sf_warning(" %d/%d;", ir, nr);
	if (NULL != mask) {
	    sf_floatread(u,n123,mask);
	    mask32 (false, 1, 1, 1, n[0], n[1], n[2], u, m);
	}
	sf_floatread(u,n123,in);
	fdip( u, p, m, n4);
	if( n4 == 2 )	sf_floatwrite(p, n123*2, out);
	else	sf_floatwrite(p, n123, out);
    }
    if(verb) sf_warning(" finished.\n");
    fdip_close();
    free(u);
    free(p);
    exit(0);
}



