/* 2-D metric tensor */
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

#define LOOP1( a) for(iq1=0;iq1<aq1.n;iq1++){ {a} }
#define LOOP2( a) for(iq2=0;iq2<aq2.n;iq2++){ {a} }
#define LOOP21(a) LOOP2(LOOP1( {a}))

#define PD1(v) ( v[iq2  ][iq1+1] - v[iq2  ][iq1-1]) / (2*aq1.d)
#define PD2(v) ( v[iq2+1][iq1  ] - v[iq2-1][iq1  ]) / (2*aq2.d)

int main (int argc, char *argv[])
{
    bool verb;  /* verbosity flag */
    sf_file Fm; /* file with coordinate system mapping */
    sf_file Fg; /* file with metric tensor coefficients */

    pt2d *mm; /* coordinates (x2,x1)[q1] */
    float **zz,**xx;               /* coordinates */
    float **g11,**g12,**g21,**g22; /* coefficients */

    axa aq1,aq2,ac;
    int iq1,iq2;
    /*------------------------------------------------------------*/
    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false;

    /* mapping */
    Fm = sf_input ("in");
    iaxa(Fm,&aq1,1); if(verb) raxa(aq1);
    iaxa(Fm,&aq2,2); if(verb) raxa(aq2);

    mm = (pt2d*) sf_alloc(aq1.n,sizeof(*mm));
    zz = sf_floatalloc2(aq1.n,aq2.n);
    xx = sf_floatalloc2(aq1.n,aq2.n);

    ac.n=4;
    ac.o=0;
    ac.d=1;

    /* coefficients */
    Fg = sf_output("out");
    oaxa(Fg,&aq1,1); if(verb) raxa(aq1);
    oaxa(Fg,&aq2,2); if(verb) raxa(aq2);
    oaxa(Fg,&ac ,3); if(verb) raxa(ac );

    /* set the output to float */
    sf_putint(Fg,"esize",4);
    sf_settype(Fg,SF_FLOAT);

    sf_close();

    g11 = sf_floatalloc2(aq1.n,aq2.n); LOOP21( g11[iq2][iq1]=0.; );
    g12 = sf_floatalloc2(aq1.n,aq2.n); LOOP21( g12[iq2][iq1]=0.; );
    g21 = sf_floatalloc2(aq1.n,aq2.n); LOOP21( g21[iq2][iq1]=0.; );
    g22 = sf_floatalloc2(aq1.n,aq2.n); LOOP21( g22[iq2][iq1]=0.; );
    /*------------------------------------------------------------*/
    
    /* read CS mapping */
    for( iq2=0; iq2<aq2.n; iq2++) {
	if(verb) sf_warning("read %d \n",iq2);
	pt2dread1(Fm,mm,aq1.n,2);

	LOOP1( zz[iq2][iq1] = mm[iq1].z;
	       xx[iq2][iq1] = mm[iq1].x; );
    }

    /*------------------------------------------------------------*/
    /* compute metric tensor */
    for (iq2=1; iq2<aq2.n-1; iq2++) {
	for (iq1=1; iq1<aq1.n-1; iq1++) {
	    
	    g11[iq2][iq1] = PD1(zz) * PD1(zz) + PD1(xx) * PD1(xx);
	    g12[iq2][iq1] = PD1(zz) * PD2(zz) + PD1(xx) * PD2(xx);
	    g21[iq2][iq1] = PD2(zz) * PD1(zz) + PD2(xx) * PD1(xx);
	    g22[iq2][iq1] = PD2(zz) * PD2(zz) + PD2(xx) * PD2(xx);
	}
    }
    LOOP1( g11[0][iq1] = g11[1][iq1]; g11[aq2.n-1][iq1]= g11[aq2.n-2][iq1];
	   g12[0][iq1] = g12[1][iq1]; g12[aq2.n-1][iq1]= g12[aq2.n-2][iq1];
	   g21[0][iq1] = g21[1][iq1]; g21[aq2.n-1][iq1]= g21[aq2.n-2][iq1];
	   g22[0][iq1] = g22[1][iq1]; g22[aq2.n-1][iq1]= g22[aq2.n-2][iq1];
	);
    LOOP2( g11[iq2][0] = g11[iq2][1]; g11[iq2][aq1.n-1]= g11[iq2][aq1.n-2];
	   g12[iq2][0] = g12[iq2][1]; g12[iq2][aq1.n-1]= g12[iq2][aq1.n-2];
	   g21[iq2][0] = g21[iq2][1]; g21[iq2][aq1.n-1]= g21[iq2][aq1.n-2];
	   g22[iq2][0] = g22[iq2][1]; g22[iq2][aq1.n-1]= g22[iq2][aq1.n-2];
     );

    sf_floatwrite(g11[0],aq1.n*aq2.n,Fg);
    sf_floatwrite(g12[0],aq1.n*aq2.n,Fg);
    sf_floatwrite(g21[0],aq1.n*aq2.n,Fg);
    sf_floatwrite(g22[0],aq1.n*aq2.n,Fg);

    /*------------------------------------------------------------*/
    exit (0); 

}
