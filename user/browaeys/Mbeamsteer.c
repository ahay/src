/* Beam steering. */
/* Adapted from Steve Cole, Stanford University, 1995 */
/*
  Copyright (C) 2008 University of Texas at Austin
  Original author: Steve Cole, Stanford University, 1995
  
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



/*sfbeamsteer sem=0 nslo=40 nazim=73 slomin=0. slomax=0.0033333 azmin=0. 
              azmax=360. xref=0 yref=0 mode=0 d1= d2= d3= npad=500 
              icoord=0 cord= < in.h > out.h */

/* from aux:	cord integer n1:n1coord=2 */

/*
        integer infd,outfd
	integer auxout,auxfd

	if (n1out .eq. 1) then
	to history:	integer n1:nslo, n2:nazim, n3:n1out
	to history:	real o1:slomin, d1:dslo, o2:azmin, d2:dazim
	else
	to history:	integer n1:n1out, n2:nslo, n3:nazim
	to history:	real d1:d1out
	to history:	real o2:slomin, d2:dslo, o3:azmin, d3:dazim
	endif
	to history:	integer mode,sem
	to aux:		live integer n1:n2, n2:n3 

	infd = input()
	outfd = output()
	if (icoord .ne. 0) auxfd = auxin('cord')

	call hclose()

	auxfd = auxout('live')
*/




#include <rsf.h>
#include <float.h>
#include <math.h>

#include "bsteer.h"

#define MAX(a,b) (a > b ? a : b)

int main(int argc, char* argv[])
{
    int i1,i2,i3,nlive,n1,n2,n3;
    float d1,d2,d3,o2,o3,dazim,dslo,d1out;
    float slomin,slomax,azmin,azmax,xref,yref,pmax,sum;
    int nazim,nslo;
    int npad,nrec,n1pad,lwind,n1out;
    int mode,sem,icoord;
    float ***data,***coord;
    int **live;
    float ***semb;

    sf_file in,out;
   
    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");


    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    /* time sampling interval */
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    /* spatial sampling interval in x */
    if (!sf_histfloat(in,"d3",&d3)) d3=1.;
    /* spatial sampling interval in y */

    if (!sf_histfloat(in,"o2",&o2)) o2=0.;
    if (!sf_histfloat(in,"o3",&o3)) o3=0.;


    if (!sf_getint("sem",&sem)) sem=0;
    /* =0 result is stack along beam steering trajectory.
       =1 result is semblance computed over beam steering trajectory */

    if (!sf_getint("mode",&mode)) mode=0;
    /* =0  beams are computed as a function of apparent slowness and azimuth angle. 
       =1  changes to px and py. Nslo, slomin, slomax, nazim, azmin, azmax are 
           interpreted as npx, pxmin, pxmax, npy, pymin, pymax. 
           This is effectively a slant stack. In this case the default parameter 
           values are not appropriate, so you need to supply the correct values */

    if (!sf_getint("npad",&npad)) npad=500;
    /* Program pads the input traces with this many samples at each end.
       This is done to avoid having to check for beams going outside of 
       available data. Does not affect the output. */

    if (!sf_getint("lwind",&lwind)) lwind=1;
    /* Semblance values computed over time windows lwind samples long */

    if (!sf_getint("icoord",&icoord)) icoord=0;
    /* =0 receiver coordinates are computed from o2,d2,o3,d3. */
      
    if (!sf_getfloat("xref",&xref)) xref = 0.;
    /* Coordinates of the point at which the beams are computed */

    if (!sf_getfloat("yref",&yref)) yref = 0.;
    /* Default is center of the array if no coordinates are supplied.
       If a coordinate file is supplied, need to specify this.*/

    if (!sf_getint("nazim",&nazim)) nazim=20;
    /* Number of azimuth angles */
    if (!sf_getfloat("azmax",&azmax)) azmax = 360.0;
    if (!sf_getfloat("azmin",&azmin)) azmin = 0.0;
    if (!sf_getfloat("dazim",&dazim)) dazim = (azmax - azmin) / (nazim-1);

    if (!sf_getint("nslo",&nslo)) nslo=20;
    /* Number of slownesses in output */
    if (!sf_getfloat("slomax",&slomax)) slomax = 0.003333;
    if (!sf_getfloat("slomin",&slomin)) slomin = 0.0;
    if (!sf_getfloat("dslo",&dslo)) dslo = (slomax - slomin) / (nslo-1);


    n1pad = n1 + 2*npad;
    nrec = n2*n3;
    if (mode == 1) pmax = MAX(slomax,azmax);
    n1out = n1/lwind;
    d1out = d1*lwind;

    /* Input file and auxiliary */
    data = sf_floatalloc3(n1pad,n2,n3);
    live = sf_intalloc2(n2,n3);
    coord = sf_floatalloc3(3,n2,n3);

    /* Output file */ 
    semb = sf_floatalloc3(nazim,nslo,n1out);

    /* clear the array */
    for (i3 = 0; i3 < n3; i3++) {
	for (i2 = 0; i2 < n2; i2++) {
	    for (i1 = 0; i1 < n1pad; i1++) {
		data[i1][i2][i3] = 0.;
	    }
	}
    }

    /* read in the data */
    sf_floatread(data[0][0],n1pad*n2*n3,in);

    /* determine if this is a dead or live trace (saves time later) */
    nlive = 0;
    for (i3 = 0; i3 < n3; i3++) {
	for (i2 = 0; i2 < n2; i2++) {
	    sum = 0.;
	    for (i1 = npad-1; i1 < npad+n1; i1++) sum = sum + data[i1][i2][i3];

	    if (sum == 0.) {
		live[i2][i3] = 0;
	    } else {
		live[i2][i3] = 1;
		nlive = nlive + 1;
	    }
	}
    }

    /* beam steering */
    bsteer(data,n1pad,n2,n3,nslo,nazim,slomin,dslo,azmin,dazim,d1,d2,d3,live,mode,n1,npad,nlive,sem,pmax,lwind,n1out,icoord,coord,xref,yref,semb);

    /* output beam */
    sf_floatwrite (semb[0][0],nazim*nslo*n1out,out);

    exit(0);
}
