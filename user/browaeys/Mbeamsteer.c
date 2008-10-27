/* Beam steering for 2D surface seismic array. */
/*
  Copyright (C) 2008 University of Texas at Austin
  Adapted from Steve Cole, Stanford University, 1995

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
#include <float.h>
#include <math.h>

#include "bsteer.h"

#define MAX(a,b) (a > b ? a : b)

int main(int argc, char* argv[])
{
    int n1,n2,n3,i1,i2,i3,nlive;
    float o1,o2,o3,d1,d2,d3;

    int nazim,nslo;
    float dazim,dslo,d1out;
    float slomin,slomax,azmin,azmax;

    float xref,yref;

    int npad,n1pad,lwind,n1out,mode;
    float pmax,sum;

    float ***data, ***semb;
    int **live;

    sf_axis axout,axslo,axazi;
    sf_file in,out;
   
    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&n1pad)) sf_error("No n1pad= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    /* time sampling interval */
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    /* spatial sampling interval in x */
    if (!sf_histfloat(in,"d3",&d3)) d3=1.;
    /* spatial sampling interval in y */

    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;
    if (!sf_histfloat(in,"o3",&o3)) o3=0.;

    if (!sf_getint("mode",&mode)) mode=1;
    /* =0 beams are computed as a function of apparent slowness and azimuth angle. 
       =1 changes from (slo,azim) to (px,py).  Parameters (nslo,slomin,slomax)
          and (nazim,azmin,azmax) are (npx,pxmin,pxmax) and (npy,pymin,pymax). 
          This is effectively a slant stack. */

    if (!sf_getint("npad",&npad)) sf_error("No npad= in input");
    /* Number previously used to pad input traces with this 
       many samples at each end, to avoid having to check for 
       beams going outside of available data. */

    if (!sf_getint("n1",&n1)) sf_error("No n1= in input");
    /* Number n1 of unpadded input traces.*/

    if (!sf_getint("lwind",&lwind)) lwind=20;
    /* Stack values computed over time windows lwind samples long. */
      
    if (!sf_getfloat("xref",&xref)) sf_error("No xref= in input");
    if (!sf_getfloat("yref",&yref)) sf_error("No yref= in input");
    /* Coordinates of the point at which the beams are computed. */

    if (!sf_getint("nslo",&nslo)) nslo=21;
    /* number of px or slownesses. */
    if (!sf_getfloat("slomax",&slomax)) slomax =  4.0;
    if (!sf_getfloat("slomin",&slomin)) slomin = -4.0;
    if (!sf_getfloat("dslo",&dslo)) dslo = (slomax-slomin)/(nslo-1);

    if (!sf_getint("nazim",&nazim)) nazim=21;
    /* number of py or azimuth. */
    if (!sf_getfloat("azmax",&azmax)) azmax =  4.0;
    if (!sf_getfloat("azmin",&azmin)) azmin = -4.0;
    if (!sf_getfloat("dazim",&dazim)) dazim = (azmax-azmin)/(nazim-1);

    if (mode == 1) {
	pmax = MAX(slomax,azmax);
    } else {
	pmax = 1000000.0;
    }

    if ( !(n1 == (n1pad - 2*npad)) ) sf_error("Unconsistent n1, npad and n1pad");

    /* Number of time windows in output */
    /* n1out=1 when stacking over all time samples */ 
    n1out = n1/lwind;
    d1out = d1*lwind;

    /* output */ 
    semb = sf_floatalloc3(n1out,nslo,nazim);

    axazi = sf_maxa(nazim,dazim,azmin);
    axslo = sf_maxa(nslo,dslo,slomin);
    axout = sf_maxa(n1out,d1out,o1);
    
    sf_oaxa(out,axazi,1);
    sf_oaxa(out,axslo,2);
    sf_oaxa(out,axout,3);

    /* input */
    data = sf_floatalloc3(n3,n2,n1pad);

    /* auxiliary */
    live = sf_intalloc2(n3,n2);

    /* clear the array */
    for (i3 = 0; i3 < n3; i3++) {
	for (i2 = 0; i2 < n2; i2++) {
	    for (i1 = 0; i1 < n1pad; i1++) {
		data[i3][i2][i1] = 0.;
	    }
	}
    }

    /* read in the (paddded) data */
    sf_floatread(data[0][0],n3*n2*n1pad,in);

    /* determine if this is a dead or live trace (saves time later) */
    nlive = 0;
    for (i3 = 0; i3 < n3; i3++) {
	for (i2 = 0; i2 < n2; i2++) {
	    sum = 0.;
	    for (i1 = npad-1; i1 < npad+n1-1; i1++) sum += data[i3][i2][i1];
	    if (sum == 0.) {
		live[i3][i2] = 0;
	    } else {
		live[i3][i2] = 1;
		nlive += 1;
	    }
	}
    }

    /* beam steering */
    bsteer(data,n1pad,n2,n3,nslo,nazim,slomin,dslo,azmin,dazim,d1,d2,d3,o2,o3,live,mode,n1,npad,nlive,pmax,lwind,n1out,xref,yref,semb);

    /* output beam */
    sf_floatwrite (semb[0][0],n1out*nslo*nazim,out);

    exit(0);
}
