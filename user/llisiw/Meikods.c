/* Fast marching with source perturbation. */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

#include "eikods.h"

int main (int argc,char* argv[]) 
{
    int b1, b2, b3, n1, n2, n3, i, nshot, ndim, is,order,n123, *p, l;
    float br1, br2, br3, o1, o2, o3, d1, d2, d3, slow;
    float **s, *t, *v, *dl1, *ds1, *dl2=NULL, *ds2=NULL;
    char *sfile;
    bool isvel, sweep, plane[3];
    sf_file vel, time, shots, tdl1, tds1, tdl2=NULL, tds2=NULL;

    sf_init (argc, argv);
    vel = sf_input("in");
    time = sf_output("out");

    if (SF_FLOAT != sf_gettype(vel)) 
	sf_error("Need float input");
    if(!sf_histint(vel,"n1",&n1)) sf_error("No n1= in input");
    if(!sf_histint(vel,"n2",&n2)) sf_error("No n2= in input");
    if(!sf_histint(vel,"n3",&n3)) n3=1;

    if(!sf_histfloat(vel,"d1",&d1)) sf_error("No d1= in input");
    if(!sf_histfloat(vel,"d2",&d2)) sf_error("No d2= in input");
    if(!sf_histfloat(vel,"d3",&d3)) d3=d2;

    if(!sf_histfloat(vel,"o1",&o1)) o1=0.;
    if(!sf_histfloat(vel,"o2",&o2)) o2=0.;
    if(!sf_histfloat(vel,"o3",&o3)) o3=0.;

    if(!sf_getbool("vel",&isvel)) isvel=true;
    /* if y, the input is velocity; n, slowness squared */

    if(!sf_getint("order",&order)) order=2;
    /* [1,2] Accuracy order */

    if (!sf_getbool("sweep",&sweep)) sweep=false;
    /* if y, use fast sweeping instead of fast marching */

    if(!sf_getfloat("br1",&br1)) br1=d1;    
    if(!sf_getfloat("br2",&br2)) br2=d2; 
    if(!sf_getfloat("br3",&br3)) br3=d3;
    /* Constant-velocity box around the source (in physical dimensions) */
 
    if(!sf_getbool("plane1",&plane[2])) plane[2]=false;
    if(!sf_getbool("plane2",&plane[1])) plane[1]=false;
    if(!sf_getbool("plane3",&plane[0])) plane[0]=false;
    /* plane-wave source */

    if(!sf_getint("b1",&b1)) b1= plane[2]? n1: (int) (br1/d1+0.5); 
    if(!sf_getint("b2",&b2)) b2= plane[1]? n2: (int) (br2/d2+0.5); 
    if(!sf_getint("b3",&b3)) b3= plane[0]? n3: (int) (br3/d3+0.5); 
    /* Constant-velocity box around the source (in samples) */

    if( b1<1 ) b1=1;  
    if( b2<1 ) b2=1;  
    if( b3<1 ) b3=1;

    sfile = sf_getstring("shotfile");
    /* File with shot locations (n2=number of shots, n1=3) */

    if(NULL != sfile) {
	shots = sf_input("shotfile");

	if (SF_FLOAT != sf_gettype(shots)) 
	    sf_error("Need float shotfile");
	if(!sf_histint(shots,"n2",&nshot)) 
	    sf_error("No n2= in shotfile");
	if(!sf_histint(shots,"n1",&ndim) || ndim != 3) 
	    sf_error("Need n1=3 in shotfile");
  
	s = sf_floatalloc2 (ndim,nshot);
	sf_floatread(s[0],nshot*ndim,shots);
	sf_fileclose(shots);
    
	sf_putint (time,"n4",nshot);
	free (sfile);
    } else {
	nshot = 1;
	ndim = 3;
    
	s = sf_floatalloc2 (ndim,nshot);     

	if(!sf_getfloat("zshot",&s[0][0])  ) s[0][0]=0.; 
	/* Shot location (used if no shotfile) */
	if(!sf_getfloat("yshot",&s[0][1])) s[0][1]=o2 + 0.5*(n2-1)*d2;
	if(!sf_getfloat("xshot",&s[0][2])) s[0][2]=o3 + 0.5*(n3-1)*d3;

	sf_warning("Shooting from zshot=%g yshot=%g xshot=%g",
		   s[0][0],s[0][1],s[0][2]);
    }

    n123 = n1*n2*n3;

    t  = sf_floatalloc (n123);
    v  = sf_floatalloc (n123);
    p  = sf_intalloc   (n123);

    sf_floatread(v,n123,vel);
    if (isvel) {
	/* transform velocity to slowness squared */
	for(i = 0; i < n123; i++) {
	    slow = v[i];
	    v[i] = 1./(slow*slow);
	}
    } 
    
    /* first-order derivative */	
    if (NULL == sf_getstring("tdl1"))
	sf_error("Output derivative tdl1= missing.");
    if (NULL == sf_getstring("tds1"))
	sf_error("Output derivative tds1= missing.");
    tdl1 = sf_output("tdl1");
    tds1 = sf_output("tds1");
    dl1  = sf_floatalloc(n123);
    ds1  = sf_floatalloc(n123);
    if (nshot != 1) sf_putint(tdl1,"n4",nshot);
    if (nshot != 1) sf_putint(tds1,"n4",nshot);

    /* second-order derivative */    
    if (NULL != sf_getstring("tdl2")) {
	tdl2 = sf_output("tdl2");
	if (nshot != 1) sf_putint(tdl2,"n4",nshot);
    }
    if (NULL != sf_getstring("tds2")) {
	tds2 = sf_output("tds2");
	if (nshot != 1) sf_putint(tds2,"n4",nshot);
    }

    if (tdl2 != NULL || tds2 != NULL) {
	dl2  = sf_floatalloc(n123);
	ds2  = sf_floatalloc(n123);
    }
    
    if(!sf_getint("l",&l)) l=1;
    /* source perturbation direction */

    if (!sweep) eikods_init(n3,n2,n1);
 
    /* loop over shots */
    for( is = 0; is < nshot; is++) {
	sf_warning("shot %d of %d;",is+1,nshot);
	if (sweep) {
	    continue;
	} else {
	    eikods(t,v,p, plane,
		   n3,n2,n1,
		   o3,o2,o1,
		   d3,d2,d1,
		   s[is][2],s[is][1],s[is][0], 
		   b3,b2,b1,
		   order,l,
		   dl1,ds1,dl2,ds2);
	}	

	sf_floatwrite(t,n123,time);
	sf_floatwrite(dl1,n123,tdl1);
	sf_floatwrite(ds1,n123,tds1);
	if (tdl2 != NULL) sf_floatwrite(dl2,n123,tdl2);
	if (tds2 != NULL) sf_floatwrite(ds2,n123,tds2);
    }
    sf_warning(".");

    exit (0);
}

/* 	$Id: Meikonal.c 7107 2011-04-10 02:04:14Z ivlad $	 */
