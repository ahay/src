/* VTI eikonal solver (3-D). */
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

#include <math.h>

#include<rsf.h>

#include "wkbjti.h"

int main (int argc,char* argv[]) 
{
    int b1, b2, b3, n1, n2, n3, i, nshot, ndim, is,n123, /* **p, */ j;
    float br1, br2, br3, o1, o2, o3, d1, d2, d3;
    float **s, **t, **v, **vv,**vs,**q,**a,**sg,**bet;
    char *sfile, *file;
    bool plane[3];
    sf_file vvf, bv,  eta, time, shots;

    sf_init (argc, argv);
    bv = sf_input("in");
    time = sf_output("out");

    if (SF_FLOAT != sf_gettype(bv)) 
	sf_error("Need float input");
    if (NULL != (file = sf_getstring("vv"))) {
	vvf = sf_input(file);
	free(file);
    } else {
	vvf = NULL;
    }
    if (NULL != (file = sf_getstring("eta"))) {
	eta = sf_input(file);
	free(file);
    } else {
	eta = NULL;
    }
    if(!sf_histint(bv,"n1",&n1)) sf_error("No n1= in input");
    if(!sf_histint(bv,"n2",&n2)) sf_error("No n2= in input");
    if(!sf_histint(bv,"n3",&n3)) n3=1;
    if (n1 <2 || n2<2) sf_error("n1 and n2 muxt be bigger than 1");

    if(!sf_histfloat(bv,"d1",&d1)) sf_error("No d1= in input");
    if(!sf_histfloat(bv,"d2",&d2)) sf_error("No d2= in input");
    if(!sf_histfloat(bv,"d3",&d3)) d3=d2;

    if(!sf_histfloat(bv,"o1",&o1)) o1=0.;
    if(!sf_histfloat(bv,"o2",&o2)) o2=0.;
    if(!sf_histfloat(bv,"o3",&o3)) o3=0.;

    /* if y, the input is background time; n, Velocity */

    /* The value of the constant eta */

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

    t  = sf_floatalloc2(n1,n2);
    v  = sf_floatalloc2(n1,n2);
    vv  = sf_floatalloc2(n1,n2);
    vs  = sf_floatalloc2(n1,n2);
    q  = sf_floatalloc2(n1,n2);
    a  = sf_floatalloc2(n1,n2);
    bet = sf_floatalloc2(n1,n2);
    sg = sf_floatalloc2(n1,n2); 
    /* p  = sf_intalloc2(n1,n2); */

    sf_floatread(v[0],n123,bv);
	

    if (NULL != vvf) {
	sf_floatread(vv[0],n123,vvf);
	sf_fileclose(vvf);
	/* transform velocity to slowness squared */
	for(j = 0; j < n2; j++)
	  for(i = 0; i < n1; i++)
	    vv[j][i] = vv[j][i]*vv[j][i];
    } else {
      for(j = 0; j < n2; j++)
	for(i = 0; i < n1; i++)
	    vv[j][i] = v[j][i]*v[j][i];
    }

    if(NULL != eta) {
	sf_floatread(q[0],n123,eta);
	sf_fileclose(eta);
	/* transform eta to q */
    } else { /* assume elliptic anisotropy */
      for(j = 0; j < n2; j++)
	for(i = 0; i < n1; i++)
	    q[j][i] = 0.;
    }

    for(j = 0; j < n2; j++)
      for(i = 0; i < n1; i++){
	v[j][i] = v[j][i]*v[j][i]*(1+2*q[j][i]);
	vs[j][i] = 0.25*vv[j][i];
	q[j][i] = vv[j][i]-2*vs[j][i];
      };

    /* loop over shots and find traveltime for source perturbation*/
    for( is = 0; is < nshot; is++) {
      wkbjti(s[is][1],s[is][0],n1,d1,o1,n2,d2,o2,v,vv,vs,q,t,a,sg,bet);
	sf_floatwrite (t[0],n123,time);
    }

    /* close input */
    sf_fileclose(bv);
    
    exit (0);
}

/* 	$Id: MwkbjTI.c 4136 2009-02-07 17:20:32Z sfomel $	 */
