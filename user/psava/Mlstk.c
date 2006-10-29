/* Local slant stacks (2D) */
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

    sf_file Fi=NULL; /* input   */
    sf_file Fo=NULL; /* output  */

    float  **dd=NULL; /* data  */
    float  **ss=NULL; /* slant-stack */
    float   *gg=NULL; /*  taper */
    float ***ww=NULL; /* weight */
    int   ***k1=NULL; /* slant-stack index on axis 1*/
    int   ***k2=NULL; /* slant-stack index on axis 2*/

    /* cube axes */
    sf_axis a1,a2,aa,ll;
    int     n1,n2,na,nl;
    int     i1,i2,ia,il;
    float oa,da,a;
    float ol, dl,l,l1,l2;
    float o1,d1,f1;
    float o2,d2,f2;

    int ic;
    int ompchunk; 

/*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */

    Fi = sf_input ("in" );
    Fo = sf_output("out");

    /* angle axis (in degrees) */
    if(! sf_getint  ("na",&na)) na=1;
    if(! sf_getfloat("oa",&oa)) oa=0.0;
    if(! sf_getfloat("da",&da)) da=1.0;
    aa=sf_maxa(na,oa,da);
    sf_setlabel(aa,""); 
    sf_setunit (aa,"");

    /* length axis (in samples) */
    if(! sf_getint  ("nl",&nl)) nl=1;
    if(! sf_getfloat("dl",&dl)) dl=1.;
    if(! sf_getfloat("ol",&ol)) ol=0.;
    ll=sf_maxa(nl,ol,dl);
    sf_setlabel(ll,""); 
    sf_setunit (ll,"");

    /* input axes */
    a1 = sf_iaxa(Fi,1); sf_setlabel(a1,""); sf_setunit(a1,"");
    a2 = sf_iaxa(Fi,2); sf_setlabel(a2,""); sf_setunit(a2,""); 
    n1 = sf_n(a1); o1=sf_o(a1); d1=sf_d(a1);
    n2 = sf_n(a2); o2=sf_o(a2); d2=sf_d(a2);

    if(verb) {
	sf_raxa(a1);
	sf_raxa(a2);
	sf_raxa(aa);
	sf_raxa(ll);
    }    

    /* setup output header */
    sf_oaxa(Fo,a1,1);
    sf_oaxa(Fo,a2,2);
    sf_oaxa(Fo,aa,3); /* angle */

/*------------------------------------------------------------*/

    /* allocate arrays */
    dd=sf_floatalloc2(n1,n2); /* data */
    ss=sf_floatalloc2(n1,n2); /* slant-stack */ 

    gg=sf_floatalloc (  2*nl);
    ww=sf_floatalloc3(4,2*nl,na);
    k1=sf_intalloc3  (4,2*nl,na);
    k2=sf_intalloc3  (4,2*nl,na);
    
/*------------------------------------------------------------*/
    /* taper */
    for(il=0;il<2*nl;il++) {
	l = ol + (il-nl)*dl;
	l /= (nl/2);
	gg[il] = exp(-l*l);
    }

/*------------------------------------------------------------*/
    /* compute indices and weight */
    for(ia=0;ia<na;ia++) {
	a  = oa + ia * da;
	a *= SF_PI/180.;
	
	for(il=0;il<2*nl;il++){
	    l = ol + (il-nl)*dl;   
	    
	    l1 = l*sin(a);
	    l2 = l*cos(a);	    
	    
	    k1[ia][il][0] = (floor)(l1);
	    k1[ia][il][1] = k1[ia][il][0] +1;
	    k1[ia][il][2] = k1[ia][il][0];
	    k1[ia][il][3] = k1[ia][il][0] +1;
	    
	    k2[ia][il][0] = (floor)(l2);
	    k2[ia][il][1] = k2[ia][il][0];
	    k2[ia][il][2] = k2[ia][il][0] +1;
	    k2[ia][il][3] = k2[ia][il][0] +1;
	    
	    f1 = l1-k1[ia][il][0];
	    f2 = l2-k2[ia][il][0];
	    
	    ww[ia][il][0] = (1-f1)*(1-f2);
	    ww[ia][il][1] =    f1 *(1-f2);
	    ww[ia][il][2] = (1-f1)*   f2 ;
	    ww[ia][il][3] =    f1 *   f2 ;

	    for(ic=0;ic<4;ic++) {
		ww[ia][il][ic] *= gg[il];
	    }
	}
    }
    
/*------------------------------------------------------------*/

    /* read input */
    sf_floatread(dd[0],n1*n2,Fi);
    
    /* loop over angles */
    if(verb) fprintf(stderr,"\n");
    for(ia=0;ia<na;ia++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",ia);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(i1,i2,il,ic) shared(ia,n1,n2,nl,ss,ww,dd,k1,k2)
#endif
	for(    i2=nl; i2<n2-nl; i2++) {
	    for(i1=nl; i1<n1-nl; i1++) {	
		
		ss[i2][i1]=0;
		for(il=0;il<2*nl;il++) {
		    for(ic=0;ic<4;ic++) {
			ss[i2][i1] += ww[ia][il][ic] *
			    dd[ i2+k2[ia][il][ic] ][ i1+k1[ia][il][ic] ];
		    }
		}
		
	    }
	}

	/* write output */
	sf_floatwrite(ss[0],n1*n2,Fo);
    }
    if(verb) fprintf(stderr,"\n");

    exit (0);
}
