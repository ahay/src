/* Kirchhoff 2-D modeling for linear slowness squared. */
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

int main(int argc, char* argv[]) 
{
    int nx, nt, ns, nr, is, ir;
    float *rfl, *crv, *recv, *shot, *trace;
    float dx, x0, dt, t0, ds, s0, dr, r0;
    sf_file refl, curv, modl, shots, recvs;

    sf_init(argc,argv);
    curv = sf_input("in");
    modl = sf_output("out");

    if (SF_FLOAT != sf_gettype(curv)) sf_error("Need float input");
    if (!sf_histint(curv,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(curv,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(curv,"o1",&x0)) sf_error("No o1= in input");

    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    /* time samples */
    if (!sf_getfloat("dt",&dt)) dt=0.004;
    /* time sampling */
    if (!sf_getfloat("t0",&t0)) t0=0.;
    /* time origin */
    
    if (NULL != sf_getstring("shots")) {
	shots = sf_input("shots");
	
	if (!sf_histint(shots,"n1",&ns)) sf_error("No n1= in shots");
    } else {
	shots = NULL;

	if (!sf_getint("ns",&ns)) ns=nx;
	/* number of shots */
	if (!sf_getfloat("s0",&s0)) s0=x0;
	/* first shot */
	if (!sf_getfloat("ds",&ds)) ds=dx;
	/* shot increment */
    }

    shot = sf_floatalloc(ns);

    if (NULL != shots) {
	sf_floatread(shot,ns,shots);
    } else {
	for (is=0; is < ns; is++) {
	    shot[is] = s0+is*ds;
	}
    }

       
    exit(0);
}
    
