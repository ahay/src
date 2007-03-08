/* Local slant stacks (3D) */
/*
  Copyright (C) 2006 Colorado School of Mines
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    bool verb;
    float sig;
    int  vers;

    sf_file Fs=NULL; /* input   */
    sf_file Fr=NULL; /* input   */
    sf_file Fi=NULL; /* output  */

    float  ***us=NULL; /* source wavefield */
    float  ***ur=NULL; /* receiver wavefield */

    float  ***ts=NULL; /* source slant-stack */
    float  ***tr=NULL; /* source slant-stack */

    float   **ii=NULL; /* image */

    float    *gg=NULL; /*  taper */
    float ****ww=NULL; /* weight */
    int   ****k1=NULL; /* slant-stack index on axis 1*/
    int   ****k2=NULL; /* slant-stack index on axis 2*/
    int   ****k3=NULL; /* slant-stack index on axis 3*/

    /* cube axes */
    sf_axis a1,a2,a3,aa,bb,ll,ak;
    int     n1,n2,n3,na,nb,nl;
    int     i1,i2,i3,ia,ib,il;
    int     j1,j2,j3;
    int     m1,m2,m3;
    int     h1,h2,h3;
    int     g1,g2,g3;
    float wo;
    float oa,da,a;
    float ob,db,b;
    float ol, dl,l,l1,l2,l3;
    float o1,d1,f1;
    float o2,d2,f2;
    float o3,d3,f3;

    const int nc=8;
    int       ic;
    int ompchunk; 

/*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */
    if(! sf_getfloat("sig",&sig))    sig=1.0;
    if(! sf_getint("vers",&vers)) vers=0;

    Fs = sf_input ("in" );
    Fr = sf_input ("ur" );
    Fi = sf_output("out");

    /* angle axis (in degrees) */
    if(! sf_getint  ("na",&na)) na=1;
    if(! sf_getfloat("oa",&oa)) oa=0.0;
    if(! sf_getfloat("da",&da)) da=1.0;
    aa=sf_maxa(na,oa,da);
    sf_setlabel(aa,"a"); 
    sf_setunit (aa,"");

    /* angle axis (in degrees) */
    if(! sf_getint  ("nb",&nb)) nb=1;
    if(! sf_getfloat("ob",&ob)) ob=0.0;
    if(! sf_getfloat("db",&db)) db=1.0;
    bb=sf_maxa(nb,ob,db);
    sf_setlabel(bb,"b"); 
    sf_setunit (bb,"");

    /* length axis (in samples) */
    if(! sf_getint  ("nl",&nl)) nl=0;
    if(! sf_getfloat("dl",&dl)) dl=1.;
    if(! sf_getfloat("ol",&ol)) ol=0.;
    ll=sf_maxa(nl,ol,dl);
    sf_setlabel(ll,""); 
    sf_setunit (ll,"");

    /* input axes */
    a1 = sf_iaxa(Fs,1);
    a2 = sf_iaxa(Fs,2);
    a3 = sf_iaxa(Fs,3);

    n1 = sf_n(a1); o1=sf_o(a1); d1=sf_d(a1);
    n2 = sf_n(a2); o2=sf_o(a2); d2=sf_d(a2);
    n3 = sf_n(a3); o3=sf_o(a3); d3=sf_d(a3);

    if(verb) {
	sf_raxa(a1);
	sf_raxa(a2);
	sf_raxa(a3);
	sf_raxa(aa);
	sf_raxa(bb);
	sf_raxa(ll);
    }    

    /* set output axes */
    ak=sf_maxa(1,0,1); 

    /* setup output header */
    sf_oaxa(Fi,a1,1);
    sf_oaxa(Fi,a2,2);
    sf_oaxa(Fi,ak,3);

/*------------------------------------------------------------*/

    /* allocate arrays */
    us=sf_floatalloc3(n1,n2,n3); /* source   wavefield */
    ur=sf_floatalloc3(n1,n2,n3); /* receiver wavefield */

    ts=sf_floatalloc3(n1,n2,n3); /* source   slant-stack */ 
    tr=sf_floatalloc3(n1,n2,n3); /* receiver slant-stack */ 

    ii=sf_floatalloc2(n1,n2); /* image */ 

    gg=sf_floatalloc (   2*nl+1);
    ww=sf_floatalloc4(nc,2*nl+1,na,nb);
    k1=sf_intalloc4  (nc,2*nl+1,na,nb);
    k2=sf_intalloc4  (nc,2*nl+1,na,nb);
    k3=sf_intalloc4  (nc,2*nl+1,na,nb);
    
/*------------------------------------------------------------*/
    /* taper */
    for(il=0;il<2*nl+1;il++) {
	l = ol + (il-nl)*dl;
	l /= (nl/2);
	l/= sig;
	gg[il] = exp(-l*l);
    }

/*------------------------------------------------------------*/
    /* compute bilinear interpolation indices and weights */
    for(ib=0;ib<nb;ib++) {     // b=0 means constant coordinate 1
	b  = ob + ib * db;
	b *= SF_PI/180.;
	
	for(ia=0;ia<na;ia++) { // a=0 means constant coordinate 2
	    a  = oa + ia * da;
	    a *= SF_PI/180.;
	    
	    for(il=0;il<2*nl+1;il++){
		l = ol + (il-nl)*dl;   
		
/*		l1 = l*sin(b);        // z*/
/*		l2 = l*cos(b)*cos(a); // x  */
/*		l3 = l*cos(b)*sin(a); // t*/

		l1 = l * cos(a) * sin(b);
		l2 = l * cos(a) * cos(b);
		l3 = l * sin(a);

		/* 
		   b=0: SS in   x-t plane 
		   a=0: SS in z-x   plane
		*/

		f3 = l3-floor(l3);
		f2 = l2-floor(l2);
		f1 = l1-floor(l1);

		// 000
		k3[ib][ia][il][0] = (floor)(l3);
		k2[ib][ia][il][0] = (floor)(l2);
		k1[ib][ia][il][0] = (floor)(l1);
		ww[ib][ia][il][0] = (1-f1)*(1-f2)*(1-f3);

		// 001
		k3[ib][ia][il][1] = k3[ib][ia][il][0];
		k2[ib][ia][il][1] = k2[ib][ia][il][0];
		k1[ib][ia][il][1] = k1[ib][ia][il][0]+1;
		ww[ib][ia][il][1] = (  f1)*(1-f2)*(1-f3);

		// 010
		k3[ib][ia][il][2] = k3[ib][ia][il][0];
		k2[ib][ia][il][2] = k2[ib][ia][il][0]+1;
		k1[ib][ia][il][2] = k1[ib][ia][il][0];  
		ww[ib][ia][il][2] = (1-f1)*(  f2)*(1-f3);
		
		// 011
		k3[ib][ia][il][3] = k3[ib][ia][il][0];
		k2[ib][ia][il][3] = k2[ib][ia][il][0]+1;
		k1[ib][ia][il][3] = k1[ib][ia][il][0]+1;  
		ww[ib][ia][il][3] = (  f1)*(  f2)*(1-f3);
		
		// 100
		k3[ib][ia][il][4] = k3[ib][ia][il][0]+1;
		k2[ib][ia][il][4] = k2[ib][ia][il][0];
		k1[ib][ia][il][4] = k1[ib][ia][il][0];
		ww[ib][ia][il][4] = (1-f1)*(1-f2)*(  f3);

		// 101
		k3[ib][ia][il][5] = k3[ib][ia][il][0]+1;
		k2[ib][ia][il][5] = k2[ib][ia][il][0];
		k1[ib][ia][il][5] = k1[ib][ia][il][0]+1;
		ww[ib][ia][il][5] = (  f1)*(1-f2)*(  f3);

		// 110
		k3[ib][ia][il][6] = k3[ib][ia][il][0]+1;
		k2[ib][ia][il][6] = k2[ib][ia][il][0]+1;
		k1[ib][ia][il][6] = k1[ib][ia][il][0];  
		ww[ib][ia][il][6] = (1-f1)*(  f2)*(  f3);
		
		// 111
		k3[ib][ia][il][7] = k3[ib][ia][il][0]+1;
		k2[ib][ia][il][7] = k2[ib][ia][il][0]+1;
		k1[ib][ia][il][7] = k1[ib][ia][il][0]+1;  
		ww[ib][ia][il][7] = (  f1)*(  f2)*(  f3);

		for(ic=0;ic<nc;ic++) {
		    ww[ib][ia][il][ic] *= gg[il];
		} // ic

	    } // a
	} // b
    } //l
    
/*------------------------------------------------------------*/

    sf_floatread(us[0][0],n1*n2*n3,Fs);        /* read   source wavefield */
    sf_floatread(ur[0][0],n1*n2*n3,Fr);        /* read receiver wavefield */
    
    /* loop over angles */
    if(verb) fprintf(stderr,"  b   a\n");
    if(verb) fprintf(stderr,"%3d %3d\n",nb-1,na-1);

    for(    ib=0;ib<nb;ib++) {
	for(ia=0;ia<na;ia++) {
	    if(verb) fprintf(stderr,"%3d %3d",ib,ia);
	    for(        i3=0; i3<n3; i3++) {
		for(    i2=0; i2<n2; i2++) {
		    for(i1=0; i1<n1; i1++) {	
			ts[i3][i2][i1] = 0;
			tr[i3][i2][i1] = 0;
		    }
		}
	    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(il,ic,i1,i2,i3,m1,m2,m3,j1,j2,j3,wo) shared(ia,ib,n1,n2,n3,nl,ts,tr,ww,us,ur,k1,k2,k3,h1,h2,h3,g1,g2,g3)
#endif
	    for(    il=0;il<2*nl+1;il++) {
		for(ic=0;ic<nc;    ic++) {
		    m1=k1[ib][ia][il][ic]; h1=SF_ABS(m1); g1=n1-h1; 
		    m2=k2[ib][ia][il][ic]; h2=SF_ABS(m2); g2=n2-h2; 
		    m3=k3[ib][ia][il][ic]; h3=SF_ABS(m3); g3=n3-h3; 
		    wo=ww[ib][ia][il][ic];
		    
		    for(        i3=h3; i3<g3; i3++) { j3=i3+m3;
			for(    i2=h2; i2<g2; i2++) { j2=i2+m2;
			    for(i1=h1; i1<g1; i1++) { j1=i1+m1;
				ts[i3][i2][i1] += us[j3][j2][j1] * wo;
				tr[i3][i2][i1] += ur[j3][j2][j1] * wo;
			    } // 1 loop
			}     // 2 loop
		    }         // 3 loop

		}             // c loop
	    }                 // l loop

	    for(        i3=0; i3<n3; i3++) {
		for(    i2=0; i2<n2; i2++) {
		    for(i1=0; i1<n1; i1++) {	
			ii[i2][i1] += ts[i3][i2][i1] * tr[i3][i2][i1];
		    }
		}
	    }
	    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b");	
	}                     // a loop
    }                         // b loop
    if(verb) fprintf(stderr,"\n");

    sf_floatwrite(ii[0],n1*n2,Fi);	 /* write image */

    exit (0);
}
