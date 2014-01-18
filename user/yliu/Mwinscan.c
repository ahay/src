/*   Picking scanned data window trace by trace with fixed t0
*/
/*
  Copyright (C) 2014 Jilin University
  
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

#include <math.h>

#include <rsf.h>

int main(int argc,char *argv[])
{
    int ix,iw,i,n,tn,xn,midp,winsz,starp;
    float dt,dx,t0,x0,v0,dv,f;
    float *orig,*winds;
    sf_file cmp,outf;
    
    sf_init (argc, argv); 
    cmp = sf_input("in");
    outf = sf_output("out");
    
    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&tn)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histint(cmp,"n2",&xn)) sf_error("No n2= in input");
    if (!sf_histfloat(cmp,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o2",&x0)) sf_error("No o2= in input");
    
    if(!sf_getint("winsz",&winsz)) winsz=200;
    /* for each trace,the width of window. unit:sample point */
    if(!sf_getfloat("v0",&v0)) v0=1000;
    /* init Vel for velocity scan */
    if(!sf_getfloat("deltav",&dv)) dv=20;
    /* step lenth for velocity scan */
    if(!sf_getfloat("t0",&t0)) t0=0.5;
    /* t0 fixed */
    if(!sf_getint("n",&n)) n=100;
    /* numbers of velscan*/
    
    orig = sf_floatalloc(tn*xn);
    
    winds = sf_floatalloc(winsz*xn*n);
    
    sf_putint(outf, "n1", winsz*xn);
    sf_putint(outf, "n2", n);
    sf_putfloat(outf, "d2", dv);
    sf_putfloat(outf,"o2",v0);
    /*modify output dimension*/
    
    sf_floatread(orig,tn*xn,cmp);
    for(i = 0; i < n; i++){
	for(ix = 0; ix < xn; ix++){
	    f = t0*t0+powf(ix*dx,2)/powf(v0,2);
	    midp = sqrtf(f)/dt;
	    starp = midp-(winsz/2);
	    
	    for(iw = 0; iw < winsz; iw++){
		if((starp+ix) > tn){
		    winds[iw + winsz*ix + i*winsz*xn]= 0 ;
		}else{
		    winds[iw + winsz*ix + i*winsz*xn]=orig[starp+iw + ix*tn];
		}
	    }
	}
	v0 += dv;
    }
    
    sf_floatwrite(winds,winsz*xn*n,outf);
    
    exit(0);
    
}

/* 	$Id$	 */
