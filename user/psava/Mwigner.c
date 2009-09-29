/* Assymptotic Wigner distribution in space-time */
/*
  Copyright (C) 2007 Colorado School of Mines
  
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

#ifdef _OPENMP
#include <omp.h>
#endif


/* 
 * input:  uu(n1,n2,n3)
 * output: ww(n1,n2,n3,nw)
 */

int main(int argc, char* argv[])
{
    bool verb;

    sf_file Fu,Fw; /* I/O files */
    sf_axis a1,a2,a3; /* cube axes */

    float ***uu=NULL, ***ww=NULL;
    
    int nh1,nh2,nh3;
    int ih1,ih2,ih3;
    int  n1, n2, n3;
    int  i1, i2, i3;
    int  j1, j2, j3;
    int  k1, k2, k3;

    int ompchunk;

    int lo1,hi1;
    int lo2,hi2;
    int lo3,hi3;

/*
    float h1,h2,h3,wh;
    float wk;
    sf_complex w;
*/

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */
    
    if(! sf_getint("nh1",&nh1)) nh1=0;
    if(! sf_getint("nh2",&nh2)) nh2=0;
    if(! sf_getint("nh3",&nh3)) nh3=0;
    sf_warning("nh1=%d nh2=%d nh3=%d",2*nh1+1,2*nh2+1,2*nh3+1);
    
/*
    if(! sf_getfloat("wk",&wk)) wk=0.0;
*/

    Fu = sf_input ("in" ); /*  input field */
    Fw = sf_output("out"); /* wigner distribution */

    /* read axes */
    a1=sf_iaxa(Fu,1); sf_setlabel(a1,"a1"); if(verb) sf_raxa(a1);
    a2=sf_iaxa(Fu,2); sf_setlabel(a2,"a2"); if(verb) sf_raxa(a2);
    a3=sf_iaxa(Fu,3); sf_setlabel(a3,"a3"); if(verb) sf_raxa(a3);

    n1 = sf_n(a1);
    n2 = sf_n(a2);
    n3 = sf_n(a3);

    /* allocate work arrays */
    uu=sf_floatalloc3(n1,n2,n3);
    ww=sf_floatalloc3(n1,n2,n3);

    /* read input */
    sf_floatread(uu[0][0],n1*n2*n3,Fu);

    if(verb) fprintf(stderr," h3  h2  h1\n");
    if(verb) fprintf(stderr,"%3d %3d %3d\n",2*nh3,2*nh2,2*nh1);
    for(        ih3=-nh3; ih3<nh3+1; ih3++) { lo3=SF_ABS(ih3); hi3=n3-lo3;
	for(    ih2=-nh2; ih2<nh2+1; ih2++) { lo2=SF_ABS(ih2); hi2=n2-lo2;
	    for(ih1=-nh1; ih1<nh1+1; ih1++) { lo1=SF_ABS(ih1); hi1=n1-lo1;
		if(verb) fprintf(stderr,"%3d %3d %3d",nh3+ih3,nh2+ih2,nh1+ih1);

/*		sf_warning('%d %d ',lo3,hi3);*/
		
/*		h1 = ih1 * sf_d(a1); h1*=h1;*/
/*		h2 = ih2 * sf_d(a2); h2*=h2;*/
/*		h3 = ih3 * sf_d(a3); h3*=h3;*/
/*		wh = sqrt(h1+h2+h3);*/
/*		w = cexpf( - sf_cmplx(0.,wh*wk));*/
		
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(i3,i2,i1,j3,j2,j1,k3,k2,k1) shared(ih3,ih2,ih1,lo3,hi3,lo2,hi2,lo1,hi1,uu,ww)
#endif	
		for(        i3=lo3; i3<hi3; i3++) { j3=i3-ih3; k3=i3+ih3;
		    for(    i2=lo2; i2<hi2; i2++) { j2=i2-ih2; k2=i2+ih2;
			for(i1=lo1; i1<hi1; i1++) { j1=i1-ih1; k1=i1+ih1;
			    ww[i3][i2][i1] += uu[j3][j2][j1] 
				*             uu[k3][k2][k1];
/*			    ww[i3][i2][i1] +=1;*/
			} /* n1 */
		    }     /* n2 */
		}         /* n3 */

		if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");

	    } /* nh1 */
	}     /* nh2 */
    }         /* nh3 */
    
    /* write output */
    sf_floatwrite(ww[0][0],n1*n2*n3,Fw);  

    exit (0);
}
