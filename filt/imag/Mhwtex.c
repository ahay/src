/* Huygens wavefront tracing traveltimes */
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
#include <rsf.h>
#include "hwt2d.h"

int main (int argc, char *argv[])
{
    bool verb;

    axa az,ax,at,ag;
    int       it,ig;
    
    sf_file Fv,Fs,Fw;

    float **vv; /* velocity       */
    pt2d   *wm; /* wavefront it-1 */
    pt2d   *wo; /* wavefront it   */
    pt2d   *wp; /* wavefront it+1 */

    pt2d    Ro;    /* point  on wft it-1 */
    pt2d Pm,Po,Pp; /* points on wft it   */
    pt2d    Qo;    /* point  on wft it+1 */
    
/*------------------------------------------------------------*/

    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false;

    /* velocity file */
    Fv = sf_input ("in");
    iaxa(Fv,&az,1); az.l="z"; if(verb) raxa(az);
    iaxa(Fv,&ax,2); ax.l="x"; if(verb) raxa(ax);

    vv=sf_floatalloc2(az.n,ax.n); sf_floatread(vv[0],az.n*ax.n,Fv);

    /* source location = initial wavefront file*/
    Fs = sf_input ("sou");
    iaxa(Fs,&ag,1); ag.l="g"; if(verb) raxa(ag);


    /* time axis */
    if(! sf_getint  ("nt",&at.n)) at.n=100;
    if(! sf_getfloat("ot",&at.o)) at.o=0;
    if(! sf_getfloat("dt",&at.d)) at.d=0.001;
    at.l="t";





/*------------------------------------------------------------*/

    /* wavefronts file */
    Fw = sf_output("out");
    oaxa(Fw,&ag,1); if(verb) raxa(ag);
    oaxa(Fw,&at,2); if(verb) raxa(at);

    /* set the output to complex */
    sf_putint(Fw,"esize",8);
    sf_settype(Fw,SF_COMPLEX);

/*------------------------------------------------------------*/

    /* wavefronts */
    wm = (pt2d*) sf_alloc(ag.n,sizeof(*wm));
    wo = (pt2d*) sf_alloc(ag.n,sizeof(*wo));
    wp = (pt2d*) sf_alloc(ag.n,sizeof(*wp));

    for( ig=0; ig<ag.n; ig++) {
	wm[ig].x=wo[ig].x=wp[ig].x=0;
	wm[ig].z=wo[ig].z=wp[ig].z=0;
	wm[ig].v=wo[ig].v=wp[ig].v=0;
    }

/*------------------------------------------------------------*/
    
    /* init HWT */
    hwt2d_init(vv,az,ax,at,ag);

/*------------------------------------------------------------*/

    /* read initial wavefront (it=0) */
    it=0;
    pt2dread1(Fs,wm,ag.n,2);
    for( ig=0; ig<ag.n; ig++) {
	Po = wm[ig];
	Po.v=hwt2d_getv(Po);
	wm[ig] = Po;
    }
    pt2dwrite1(Fw,wm,ag.n,2); /* write wavefront it=0 */

/*------------------------------------------------------------*/

    /* compute second wavefront (it=1) by orthogonal rays */
    it=1;

    Po=wm[0]; 
    Pp=wm[2]; 
    Ro=hwt2d_orth(Po,Po,Pp);
    wo[0] = Ro;

    for( ig=1; ig<ag.n-1; ig++) {

	Pm = wm[ig-1];
	Po = wm[ig  ];
	Pp = wm[ig+1];
	
	/* orthogonal rays */
	Ro = hwt2d_orth(Pm,Po,Pp);

	wo[ig] = Ro;
    }

    Pm=wm[ag.n-3]; 
    Po=wm[ag.n-1]; 
    Ro=hwt2d_orth(Pm,Po,Po);
    wo[ag.n-1] = Ro;

    pt2dwrite1(Fw,wo,ag.n,2); /* write wavefront it=1 */

/*------------------------------------------------------------*/

    for (it=2; it<at.n; it++) {
	if(verb) fprintf(stderr,"it=%d\n",it);

	if(ag.n>3) {
	    /* boundary */
	    ig=0;      wp[ig] = hwt2d_raytr(wm[ig],wo[ig]);

	    for (ig=1; ig<ag.n-1; ig++) {		
		Pm = wo[ig-1];
		Po = wo[ig  ];  Qo = wm[ig];
		Pp = wo[ig+1];
		
		if(hwt2d_cusp(Qo,Pm,Po,Pp)) {
		    Ro = hwt2d_raytr(Qo,   Po   );
		} else {
		    Ro = hwt2d_wfttr(Qo,Pm,Po,Pp);
		}
		wp[ig] = Ro;
	    }

	    /* boundary */
	    ig=ag.n-1; wp[ig] = hwt2d_raytr(wm[ig],wo[ig]);
	} else {
	    for (ig=0; ig<ag.n; ig++) {
		Po = wo[ig];  
		Qo = wm[ig];
		Ro = hwt2d_raytr(Qo,Po);
		wp[ig] = Ro;
	    }
	}

	/* write wavefront it */
	pt2dwrite1(Fw,wp,ag.n,2);

	/* step in time */
	for( ig=0; ig<ag.n; ig++) {
	    wm[ig] = wo[ig];
	    wo[ig] = wp[ig];
	}
    }
/*------------------------------------------------------------*/    
    exit (0);
}



