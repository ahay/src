/* */
/*
  Copyright (C) 2008 Colorado School of Mines
  
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
#include "tool3.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    bool verb;
    bool stat;

    /* I/O files */
    sf_file Fccc=NULL;
    sf_file Fspk=NULL;  /* template for derivative operator */    
    sf_file Fzdel=NULL; /* z derivative */
    sf_file Fxdel=NULL; /* x derivative */
    
    sf_axis ax,ay,az;      /* cube axes */
    sf_axis sx,sy,sz;      /* cube axes */
    
    float **c11,**c13,**c33,**c55; /* input file */
    float  mc11, mc13, mc33, mc55;
    float **xdel,**zdel;                /* output derivative */

    int ix,iz,nn;
    wfs2d wfs;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("stat",&stat)) stat=true;  /* stationary operator */

    /*------------------------------------------------------------*/
    /* input files and axes */
    Fccc = sf_input ("ccc");      
    az = sf_iaxa(Fccc,1);
    ax = sf_iaxa(Fccc,2);
    ay = sf_maxa(1,0.,1.);

    Fspk = sf_input ("in");
    sz = sf_iaxa(Fspk,1);
    sx = sf_iaxa(Fspk,2);
    sy = sf_maxa(1,0.,1.);

    /*------------------------------------------------------------*/
    /* output files*/
    Fzdel = sf_output ("zdel");
    sf_oaxa(Fzdel,sz,1);
    sf_oaxa(Fzdel,sx,2);
    sf_oaxa(Fzdel,sy,3);

    if(!stat) {
	sf_oaxa(Fzdel,az,4);
	sf_oaxa(Fzdel,ax,5);
	sf_oaxa(Fzdel,ay,6);
    }

    Fxdel = sf_output ("xdel");
    sf_oaxa(Fxdel,sz,1);
    sf_oaxa(Fxdel,sx,2);
    sf_oaxa(Fxdel,sy,3);
    
    if(!stat) {
	sf_oaxa(Fxdel,az,4);
	sf_oaxa(Fxdel,ax,5);
	sf_oaxa(Fxdel,ay,6);
    }

    /*------------------------------------------------------------*/
    if(verb) { 
	sf_raxa(sz); 
	sf_raxa(sx); 
	sf_raxa(sy);

	sf_raxa(az); 
	sf_raxa(ax); 
	sf_raxa(ay); 
    }

    /*------------------------------------------------------------*/
    /* read model params */
    c11=sf_floatalloc2(sf_n(az),sf_n(ax));
    c33=sf_floatalloc2(sf_n(az),sf_n(ax));
    c55=sf_floatalloc2(sf_n(az),sf_n(ax));
    c13=sf_floatalloc2(sf_n(az),sf_n(ax));
    
    nn = sf_n(az)*sf_n(ax);
    sf_floatread(c11[0],nn,Fccc);
    sf_floatread(c33[0],nn,Fccc);
    sf_floatread(c55[0],nn,Fccc);
    sf_floatread(c13[0],nn,Fccc);

    /*------------------------------------------------------------*/
    /* allocate arrays for derivatives */
    zdel=sf_floatalloc2(sf_n(sz),sf_n(sx));
    xdel=sf_floatalloc2(sf_n(sz),sf_n(sx));
    
    /*------------------------------------------------------------*/
    /* compute derivatives */
    wfs = wfsep_init(sx,sz);

    if(stat) {
	mc11 = sf_quantile(nn/2,nn,c11[0]);
	mc33 = sf_quantile(nn/2,nn,c33[0]);
	mc55 = sf_quantile(nn/2,nn,c55[0]);
	mc13 = sf_quantile(nn/2,nn,c13[0]);

	if(verb) {
	    sf_warning("c11=%g c33=%g c55=%g c13=%g",
		       mc11,  mc33,  mc55,  mc13);
	}

	wfsep(zdel,xdel,sx,sz,
	      mc11,mc33,mc55,mc13,
	      wfs);
	
	sf_floatwrite(zdel[0],sf_n(sz)*sf_n(sx),Fzdel);
	sf_floatwrite(xdel[0],sf_n(sz)*sf_n(sx),Fxdel);
	
    } else {

	for(    ix=0;ix<sf_n(ax);ix++){
	    for(iz=0;iz<sf_n(az);iz++){
		if(verb) fprintf(stderr,"%5d %5d",iz,ix);
		
		wfsep(zdel,xdel,sx,sz,
		      c11[ix][iz],c33[ix][iz],c55[ix][iz],c13[ix][iz],
		      wfs);
		
		sf_floatwrite(zdel[0],sf_n(sz)*sf_n(sx),Fzdel);
		sf_floatwrite(xdel[0],sf_n(sz)*sf_n(sx),Fxdel);
		
		if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b");
	    }
	}
    }

    /*------------------------------------------------------------*/
    wfsep_close(wfs);

    exit (0);
}
