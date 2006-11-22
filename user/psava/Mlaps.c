/* Lagged-products */
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

    sf_file Fs,Fr,Fi; /* I/O files */
    sf_axis a1,a2,aa; /* cube axes */

    int     n1,n2, nh1,nh2;
    int     i1,i2, ih1,ih2;
    int     j1,j2;
    int     k1,k2;

    float **ii=NULL,**us=NULL,**ur=NULL; /* arrays */

    int ompchunk; 

    int lo1,hi1;
    int lo2,hi2;

/*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */

    Fs = sf_input ("in" ); /*   source wavefield */
    Fr = sf_input ("ur" ); /* receiver wavefield */
    Fi = sf_output("out"); /* image */

    /* read axes */
    a1=sf_iaxa(Fs,1); sf_setlabel(a1,"a1"); if(verb) sf_raxa(a1);
    a2=sf_iaxa(Fs,2); sf_setlabel(a2,"a2"); if(verb) sf_raxa(a2);

    if(! sf_getint("nh1",&nh1)) nh1=0;
    if(! sf_getint("nh2",&nh2)) nh2=0;
    sf_warning("nh1=%d nh2=%d",2*nh1+1,2*nh2+1);

    /* set output axes */
    aa=sf_maxa(2*nh1+1,-nh1*sf_d(a1),sf_d(a1));
    sf_setlabel(aa,"h1");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,3);
    
    aa=sf_maxa(2*nh2+1,-nh2*sf_d(a2),sf_d(a2)); 
    sf_setlabel(aa,"h2");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,4);
    
    n1 = sf_n(a1);
    n2 = sf_n(a2);

    /* allocate work arrays */
    us=sf_floatalloc2(n1,n2);
    ur=sf_floatalloc2(n1,n2);
    ii=sf_floatalloc2(n1,n2);

    /* read input */
    sf_floatread(us[0],n1*n2,Fs);
    sf_floatread(ur[0],n1*n2,Fr);

    if(verb) fprintf(stderr," h1  h2\n");
    if(verb) fprintf(stderr,"%3d %3d\n",2*nh1,2*nh2);
    for(        ih2=-nh2; ih2<nh2+1; ih2++) { lo2=SF_ABS(ih2); hi2=n2-SF_ABS(ih2);
	for(    ih1=-nh1; ih1<nh1+1; ih1++) { lo1=SF_ABS(ih1); hi1=n1-SF_ABS(ih1);
	    if(verb) fprintf(stderr,"%3d %3d",nh2+ih2,nh1+ih1);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(i2,i1,j2,j1,k2,k1) shared(lo2,lo1,hi2,hi1,ih2,ih1,ii,us,ur)
#endif		
	    for(    i2=lo2; i2<hi2; i2++) { j2=i2-ih2; k2=i2+ih2;
		for(i1=lo1; i1<hi1; i1++) { j1=i1-ih1; k1=i1+ih1;
		    ii[i2][i1] = us[j2][j1] 
			*        ur[k2][k1];
		} // nz
	    } // nx
	    
	    /* write output */
	    sf_floatwrite(ii[0],n1*n2,Fi);    
	    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	} // nh1
    } // nh2
    if(verb) fprintf(stderr,"\n");
    
    exit (0);
}
