/* Linearized complex eikonal equation */
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
#include <time.h>

#include <rsf.h>

#include "cpxeiko.h"

int main(int argc, char* argv[])
{
    bool *m;
    int dim, i, n[SF_MAX_DIM], it, nt, *temp, niter;
    float d[SF_MAX_DIM], *tr, *ti, *w, *rhs, *output, *x0, *temp1, *temp2;
/*    double dot1[2], dot2[2]; */
    char key[4];
    sf_file realtime, imagtime, out, mask;
    
    sf_init(argc,argv);
    realtime = sf_input("in");
    imagtime = sf_input("imagtime");
    out = sf_output("out");

    dim = sf_filedims(realtime,n);

    nt = 1;
    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(realtime,key,d+i)) sf_error("No %s= in input",key);
	nt *= n[i];
    }

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    m = sf_boolalloc(nt);
    w = sf_floatalloc(nt);
    tr = sf_floatalloc(nt);
    ti = sf_floatalloc(nt);
    x0 = sf_floatalloc(nt);
    rhs = sf_floatalloc(nt);
    temp = sf_intalloc(nt);
    temp1 = sf_floatalloc(nt);
    temp2 = sf_floatalloc(nt);
    output = sf_floatalloc(nt);

    if (NULL != sf_getstring("mask"))
    {
	mask = sf_input("mask");
	if (SF_INT != sf_gettype(mask))
	    sf_error("Need int mask input.");
	sf_intread(temp,nt,mask);
	sf_fileclose(mask);
	for (it=0; it < nt; it++) {
	    if (temp[it] == 0)
		m[it] = true;
	    else
		m[it] = false;
	}
    } else {
	for (it=0; it < nt; it++) {
	    m[it] = false;
	}
    }
    /* known mask */

    for (it=0; it < nt; it++) {
	x0[it] = 0.;
    }
    /* known value: psuedo slowness update is zero */

    sf_floatread(tr,nt,realtime);
    sf_floatread(ti,nt,imagtime);

    cpxeiko_init(dim,n,nt,d,tr,ti);
    /* initiate operator */
    
    cpxeiko_forw(false,ti,w);
    /* current psuedo slowness */

    cpxeiko_forw(false,tr,temp1);
    cpxeiko_forw(true, ti,temp2);
    

    for (it=0; it < nt; it++) {
	rhs[it] = -temp1[it]-temp2[it];
    }
    /* right-hand side */

    sf_solver (cpxeiko_operator,sf_cgstep,nt,nt,output,rhs,niter,"known",m,"x0",x0,"verb",true,"end");

/*
    init_genrand((unsigned long) time(NULL));

    sf_dot_test(cpxeiko_operator,nt,nt,dot1,dot2);
    sf_warning("compare %g and %g",dot1[0],dot1[1]);
    sf_warning("compare %g and %g",dot2[0],dot2[1]);
*/

//    sf_solver(cpxeiko_operator,sf_cgstep,nt,nt,output,rhs,niter,"verb",true,"end");

    sf_cgstep_close();
    cpxeiko_close();

    sf_floatwrite(output,nt,out);
    
    exit(0);
}
