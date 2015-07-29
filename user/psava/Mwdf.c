/* Assymptotic Wigner distribution */
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
    int *ii,ik;

    int nh1,nh2,nh3;
    int nb1,nb2,nb3;

    int ih1,ih2,ih3;
    int ib1,ib2,ib3;
    int  n3;
    int  i3;
    int  j1, j2, j3;
    int  k1, k2, k3;

    int ompchunk;

    int lo1,hi1;
    int lo2,hi2;
    int lo3,hi3;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */

    Fu = sf_input ("in" ); /*  input field */
    Fw = sf_output("out"); /* wigner distribution */

    /* read axes */
    a1=sf_iaxa(Fu,1); sf_setlabel(a1,"a1"); if(verb) sf_raxa(a1);
    a2=sf_iaxa(Fu,2); sf_setlabel(a2,"a2"); if(verb) sf_raxa(a2);
    a3=sf_iaxa(Fu,3); sf_setlabel(a3,"a3"); if(verb) sf_raxa(a3);
    n3 = sf_n(a3);

    if(! sf_getint("nh1",&nh1)) nh1=0;
    if(! sf_getint("nh2",&nh2)) nh2=0;
    if(! sf_getint("nh3",&nh3)) nh3=0;
    if(n3<=1) nh3=0;
    if(verb) sf_warning("nh1=%d nh2=%d nh3=%d",2*nh1+1,2*nh2+1,2*nh3+1);

    nb1 = sf_n(a1);
    nb2 = sf_n(a2);
    nb3=2*nh3+1;

    uu=sf_floatalloc3(nb1,nb2,nb3);
    ww=sf_floatalloc3(nb1,nb2,nb3); 

    ii = sf_intalloc(nb3);
    for(ib3=0;ib3<nb3;ib3++) {
	ii[ib3]=ib3;
    }

    /*------------------------------------------------------------*/
    /* low end on axis 3 */
    /*------------------------------------------------------------*/
    for(         ib3=0; ib3<nb3; ib3++) {
	for(     ib2=0; ib2<nb2; ib2++) {
	    for( ib1=0; ib1<nb1; ib1++) {
		ww[ib3][ib2][ib1] = 0.;
	    }
	}
    }

    sf_floatread(uu[0][0],nb1*nb2*nb3,Fu);

    for(        ih3=-nh3; ih3<nh3+1; ih3++) { lo3=SF_ABS(ih3); hi3=nb3-lo3;
	for(    ih2=-nh2; ih2<nh2+1; ih2++) { lo2=SF_ABS(ih2); hi2=nb2-lo2;
	    for(ih1=-nh1; ih1<nh1+1; ih1++) { lo1=SF_ABS(ih1); hi1=nb1-lo1;
		
		for(        ib3=lo3; ib3<nh3+1;ib3++) { j3=ii[ib3-ih3]; k3=ii[ib3+ih3];
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ib1,ib2,j2,j1,k2,k1) \
    shared(ib3,j3,k3,ih2,ih1,lo2,hi2,lo1,hi1,uu,ww)
#endif
		    for(    ib2=lo2; ib2<hi2;  ib2++) { j2=ib2-ih2; k2=ib2+ih2;
			for(ib1=lo1; ib1<hi1;  ib1++) { j1=ib1-ih1; k1=ib1+ih1;
			    
			    ww[ib3][ib2][ib1] += uu[j3][j2][j1] 
				*                uu[k3][k2][k1];
			} /* nb1  */
		    }     /* nb2  */
		}         /* nb3  */

	    } /* nh1 */
	}     /* nh2 */
    }         /* nh3 */
    
    for(ih3=0;ih3<nh3+1;ih3++) {
	sf_floatwrite(ww[ih3][0],nb1*nb2,Fw);
    }

    /*------------------------------------------------------------*/
    /* loop on axis 3 */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr," n3\n");
    if(verb) fprintf(stderr,"%5d\n",n3-1);
    for(i3=nh3+1;i3<n3-nh3;i3++) {
	if(verb) fprintf(stderr,"%5d",i3);

	/* zero WDF */
	ib3=nh3;
	for(     ib2=0; ib2<nb2; ib2++) {
	    for( ib1=0; ib1<nb1; ib1++) {
		ww[ib3][ib2][ib1] = 0.;
	    }
	}

	/* circulate index to slices */
	ik = ii[0];
	for(ib3=0;ib3<nb3-1;ib3++) {
	    ii[ib3]=ii[ib3+1];
	}
	ii[nb3-1]=ik;

	/* read new slice */
	sf_floatread(uu[ ii[nb3-1] ][0],nb1*nb2,Fu);

	for(        ih3=-nh3; ih3<nh3+1; ih3++) { 
	    for(    ih2=-nh2; ih2<nh2+1; ih2++) { lo2=SF_ABS(ih2); hi2=nb2-lo2;
		for(ih1=-nh1; ih1<nh1+1; ih1++) { lo1=SF_ABS(ih1); hi1=nb1-lo1;
		    
		    ib3=nh3; j3=ii[ib3-ih3]; k3=ii[ib3+ih3];
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ib1,ib2,j2,j1,k2,k1) \
    shared(ib3,j3,k3,ih2,ih1,lo2,hi2,lo1,hi1,uu,ww)
#endif
		    for(    ib2=lo2; ib2<hi2; ib2++) { j2=ib2-ih2; k2=ib2+ih2;
			for(ib1=lo1; ib1<hi1; ib1++) { j1=ib1-ih1; k1=ib1+ih1;
			    
			    ww[ib3][ib2][ib1] += uu[j3][j2][j1] 
				*                uu[k3][k2][k1];
			} /* nb1  */
		    }     /* nb2  */
		    
		} /* nh1 */
	    }     /* nh2 */
	}         /* nh3 */
	
	sf_floatwrite(ww[nh3][0],nb1*nb2,Fw);

	if(verb) fprintf(stderr,"\b\b\b\b\b");
    }

    /*------------------------------------------------------------*/
    /* high-end on axis 3*/
    /*------------------------------------------------------------*/
    for(        ih3=-nh3; ih3<nh3+1; ih3++) { lo3=SF_ABS(ih3); hi3=nb3-lo3;
	for(    ih2=-nh2; ih2<nh2+1; ih2++) { lo2=SF_ABS(ih2); hi2=nb2-lo2;
	    for(ih1=-nh1; ih1<nh1+1; ih1++) { lo1=SF_ABS(ih1); hi1=nb1-lo1;
		
		for(        ib3=nh3+1; ib3<hi3; ib3++) { j3=ii[ib3-ih3]; k3=ii[ib3+ih3];
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) \
    private(ib1,ib2,j2,j1,k2,k1) \
    shared(ib3,j3,k3,ih2,ih1,lo2,hi2,lo1,hi1,uu,ww)
#endif
		    for(    ib2=lo2; ib2<hi2; ib2++) { j2=ib2-ih2; k2=ib2+ih2;
			for(ib1=lo1; ib1<hi1; ib1++) { j1=ib1-ih1; k1=ib1+ih1;
			    
			    ww[ib3][ib2][ib1] += uu[j3][j2][j1] 
				*                uu[k3][k2][k1];
			} /* nb1  */
		    }     /* nb2  */
		}         /* nb3  */
		
	    } /* nh1 */
	}     /* nh2 */
    }         /* nh3 */
    
    for(ih3=nh3+1;ih3<2*nh3+1;ih3++) {
	sf_floatwrite(ww[ih3][0],nb1*nb2,Fw);
    }

    /*------------------------------------------------------------*/

    exit (0);
}
