/* Conventional IC 2D */

/*
  Copyright (C) 2013 Colorado School of Mines
  
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

int main(int argc, char* argv[])
{
    bool verb,isreversed;

    sf_file Fs,Fr,Fi;    /* I/O files */
    sf_axis a1,a2,a3,aa; /* cube axes */
    int     i1,i2,i3;
    int     n1,n2,n3;

    float **us=NULL,**ur=NULL,**ii=NULL;

    float scale;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
    
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("isreversed",&isreversed)) isreversed=false; /* received wavefield */

    Fs = sf_input ("in" );
    Fr = sf_input ("uu" );
    Fi = sf_output("out");

    /*------------------------------------------------------------*/
    /* read axes */
    a1=sf_iaxa(Fs,1); if(verb) sf_raxa(a1);
    a2=sf_iaxa(Fs,2); if(verb) sf_raxa(a2);
    a3=sf_iaxa(Fs,3); if(verb) sf_raxa(a3);

    aa=sf_maxa(1,0,1); 
    sf_setlabel(aa,""); 
    sf_setunit (aa,""); 

    n1 = sf_n(a1);
    n2 = sf_n(a2);
    n3 = sf_n(a3);
    scale = 1./n3;

    /* write axes */
    sf_oaxa(Fi,a1,1);
    sf_oaxa(Fi,a2,2);
    sf_oaxa(Fi,aa,3);
    
    /*------------------------------------------------------------*/
    /* allocate arrays */
    ii = sf_floatalloc2(n1,n2); 
    us = sf_floatalloc2(n1,n2);
    ur = sf_floatalloc2(n1,n2);
    for    (i2=0; i2<n2; i2++) {
	for(i1=0; i1<n1; i1++) {
	    ii[i2][i1]=0.;
	    us[i2][i1]=0.;
	    ur[i2][i1]=0.;
	}
    }

    /*------------------------------------------------------------*/
    if(isreversed){ /* receiver wavefield is reversed */
	for (i3=0; i3<n3; i3++) {
	    if(verb) fprintf(stderr,"\b\b\b\b\b%d",i3);
	    
	    sf_floatread  (us[0],n1*n2,Fs);
	    sf_floatread  (ur[0],n1*n2,Fr);
	    
	    for    (i2=0; i2<n2; i2++) {
		for(i1=0; i1<n1; i1++) {
		    ii[i2][i1] += us[i2][i1]*ur[i2][i1];
		}
	    }
	} 
	if(verb) fprintf(stderr,"\n");
	
    } else { /* receiver wavefield is NOT reversed */
	for (i3=0; i3<n3; i3++) {
	    if(verb) fprintf(stderr,"\b\b\b\b\b%d",(n3-i3-1));
	    
	    sf_floatread  (us[0],n1*n2,Fs);

	    sf_seek(Fr,(n3-i3-1)*n1*n2*sizeof(float),SEEK_SET);
	    sf_floatread  (ur[0],n1*n2,Fr);
	    
	    for    (i2=0; i2<n2; i2++) {
		for(i1=0; i1<n1; i1++) {
		    ii[i2][i1] += us[i2][i1]*ur[i2][i1];
		}
	    }
	}
	if(verb) fprintf(stderr,"\n");

    }
       
    /*------------------------------------------------------------*/
    /* scale image */
    for    (i2=0; i2<n2; i2++) {
	for(i1=0; i1<n1; i1++) {
	    ii[i2][i1] *=scale;
	}
    }

    /*------------------------------------------------------------*/
    /* write image */
    sf_floatwrite(ii[0],n2*n1,Fi);
    
    /*------------------------------------------------------------*/
    free(*ii); free(ii);
    free(*us); free(us);
    free(*ur); free(ur);
    /*------------------------------------------------------------*/
    
    exit (0);
}
