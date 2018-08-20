/* 3-D Kirchhoff time migration with antialiasing with adjoint flag. */
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

#include "mig3.h" 

int main(int argc, char* argv[])
{
    int nt, nx, ny, n1, ntr, apt;
    char *antialias;
    float *trace, *out, rho, angle;
    float dx, dy, ox, oy, dt,ot;
    float *vel;
    bool adj, doomp;
    sf_file in, mig, velFile;

    sf_init (argc,argv);
    in = sf_input("in");
    mig = sf_output("out");
    velFile = sf_input("vel");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    ntr = sf_leftsize(in,1);

    if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    
    if (!sf_histint(in,"n2",&nx))   sf_error("Need n2=");// sf_putint(mig,"n2",nx);
    /* new inline dimension size */
    if (!sf_histfloat(in,"d2",&dx)) sf_error("Need d2=");// sf_putfloat(mig,"d2",dx);
    if (!sf_histfloat(in,"o2",&ox)) sf_error("Need o2=");// sf_putfloat(mig,"o2",ox);

    if (!sf_histint(in,"n3",&ny))   sf_error("Need n3=");// sf_putint(mig,"n3",ny);
    /* new crossline dimension size */
    if (!sf_histfloat(in,"d3",&dy)) sf_error("Need d3=");// sf_putfloat(mig,"d3",dy);
    if (!sf_histfloat(in,"o3",&oy)) sf_error("Need o3=");// sf_putfloat(mig,"o3",oy);

    if (ntr != nx*ny) sf_error("Dimensions are inconsistent");

    //if (!sf_getint("n1",&n1))   sf_error("Need n1="); sf_putint(mig,"n1",n1);
    //n1 = nt;
    /* new time dimension size */
    
    if (!sf_getbool("adj",&adj)) adj=true;
    /* adjoint flag */

    if (!sf_getbool("doomp",&doomp)) doomp=false;
    /* perform OpenMP optimization */

    if (NULL == (antialias = sf_getstring("antialias"))) antialias="triangle";
    /* antialiasing type [triangle,flat,steep,none] */

    /* Do we use nt from input or the one we specify for the new time axis? */
    if (!sf_getfloat("rho",&rho)) rho = 1.-1./nt;
    /* Leaky integration constant */

    if (!sf_getint("apt",&apt)) apt = nx;
    /* integral aperture */

    if (!sf_getfloat("angle",&angle)) angle = 90.0;
    /* angle aperture */

    //if (!sf_getfloat("vel",&vel)) sf_error("Need vel=");
    /* migration velocity */
    
    trace = sf_floatalloc(nt*nx*ny);
    out = sf_floatalloc(nt*nx*ny);

    if (adj) {

    	sf_floatread(trace,nt*nx*ny,in);

    } else {

    	sf_floatread(out,nt*nx*ny,in);

    }

    /* allocating and reading velocity */
    vel = sf_floatalloc(nt*nx*ny);
    sf_floatread(vel,nt*nx*ny,velFile);

    mig3_lop (adj,false, nt,nx,ny, dt,dx,dy, ot,ox,oy, trace, out, vel, rho, antialias[0],doomp,apt,angle);

    if (adj == true){

    	sf_floatwrite(out,nt*nx*ny,mig);
 
   } else { 

   	sf_floatwrite(trace,nt*nx*ny,mig);

   }     

    exit(0);
}
