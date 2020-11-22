/* 2-D Kirchhoff time migration with antialiasing with adjoint flag. */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

#include "mig2.h" 

int main(int argc, char* argv[])
{
    int nt, nx, ntr, apt;
    char *antialias;
    float *trace, *out, rho, angle;
    float dx, ox, dt,ot;
    float *vel;
    bool adj, doomp, ps, hd, dd;
    sf_file in, mig, fvel;

    sf_init (argc,argv);
    in = sf_input("in");
    mig = sf_output("out");
    fvel = sf_input("vel");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    ntr = sf_leftsize(in,1);

    if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    
    if (!sf_histint(in,"n2",&nx))   sf_error("Need n2=");// sf_putint(mig,"n2",nx);
    /* new inline dimension size */
    if (!sf_histfloat(in,"d2",&dx)) sf_error("Need d2=");// sf_putfloat(mig,"d2",dx);
    if (!sf_histfloat(in,"o2",&ox)) sf_error("Need o2=");// sf_putfloat(mig,"o2",ox);

    if (ntr != nx) sf_error("Dimensions are inconsistent");

    //if (!sf_getint("n1",&n1))   sf_error("Need n1="); sf_putint(mig,"n1",n1);
    //n1 = nt;
    /* new time dimension size */
    
    if (!sf_getbool("adj",&adj)) adj=true;
    /* adjoint flag */

    if (!sf_getbool("ps",&ps)) ps=true;
    /* spherical divergence */

    if (!sf_getbool("doomp",&doomp)) doomp=true;
    /* perform OpenMP optimization */

    if (NULL == (antialias = sf_getstring("antialias"))) antialias="triangle";
    /* antialiasing type [triangle,flat,steep,none] */

    /* Do we use nt from input or the one we specify for the new time axis? */
    if (!sf_getfloat("rho",&rho)) rho = 1.-1./nt;
    /* Leaky integration constant */

    if (!sf_getbool("hd",&hd)) hd=true;
    /* half derivative */

    if (!sf_getint("apt",&apt)) apt = nx;
    /* integral aperture */

    if (!sf_getfloat("angle",&angle)) angle = 90.0;
    /* angle aperture */

    if (!sf_getbool("dd",&dd)) dd = true;
    /* differentiation in the data domain */
    
    /* migration velocity */
    
    trace = sf_floatalloc(nt*nx);
    out = sf_floatalloc(nt*nx);

    if (adj) {

    	sf_floatread(trace,nt*nx,in);

    } else {

    	sf_floatread(out,nt*nx,in);

    }

    vel = sf_floatalloc(nt*nx);
    sf_floatread(vel,nt*nx,fvel);

    //sf_warning("be4 mig2_lop");

    mig2_lop (adj,false, nt,nx, dt,dx, ot,ox, trace, out, vel, rho, hd, antialias[0],doomp,apt,angle,ps,dd);

    //sf_warning("4ter mig2_lop");

    if (adj == true){

	/*for (i=0; i<nt*nx; i++){

		out[i] = trace[i];

	}*/ 

    	sf_floatwrite(out,nt*nx,mig);
 
   } else { 

	/*for (i=0; i<nt*nx; i++){

		trace[i] = out[i];

	}*/ 

   	sf_floatwrite(trace,nt*nx,mig);

   }     

    exit(0);
}
