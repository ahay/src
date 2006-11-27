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

    sf_file Fs,Fr,Fi;    /* I/O files */
    sf_axis a1,a2,a3,aa; /* cube axes */

    int     n1,n2,n3, nh1,nh2,nh3;
    int     i1,i2,i3, ih1,ih2,ih3;
    int     j1,j2,    jh1,jh2;
    int     k1,k2;

    float ****ii=NULL,**us=NULL,**ur=NULL; /* arrays */

    int ompchunk; 

    int lo1,hi1;
    int lo2,hi2;
    int lo3,hi3;

    int ictype;

/*------------------------------------------------------------*/

    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */
    if(! sf_getint("ictype",&ictype)) ictype=0;        /* I.C. type */
    
    Fs = sf_input ("in" ); /*   source wavefield */
    Fr = sf_input ("ur" ); /* receiver wavefield */
    Fi = sf_output("out"); /* image */

    /* read axes */
    a1=sf_iaxa(Fs,1); sf_setlabel(a1,"a1"); if(verb) sf_raxa(a1);
    a2=sf_iaxa(Fs,2); sf_setlabel(a2,"a2"); if(verb) sf_raxa(a2);
    a3=sf_iaxa(Fs,3); sf_setlabel(a3,"a3"); if(verb) sf_raxa(a3);

    n1 = sf_n(a1);
    n2 = sf_n(a2);
    n3 = sf_n(a3);

    if(! sf_getint("nh1",&nh1)) nh1=0;
    if(! sf_getint("nh2",&nh2)) nh2=0;
    if(! sf_getint("nh3",&nh3)) nh3=0;
    sf_warning("nh1=%d nh2=%d nh3=%d",2*nh1+1,2*nh2+1,2*nh3+1);

    /* set output axes */
    aa=sf_maxa(2*nh1+1,-nh1*sf_d(a1),sf_d(a1));
    sf_setlabel(aa,"h1");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,3);
    
    aa=sf_maxa(2*nh2+1,-nh2*sf_d(a2),sf_d(a2)); 
    sf_setlabel(aa,"h2");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,4);

    aa=sf_maxa(2*nh3+1,-nh3*sf_d(a3),sf_d(a3)); 
    sf_setlabel(aa,"h3");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,5);

    /* allocate work arrays */
    us=sf_floatalloc2(n1,n2);
    ur=sf_floatalloc2(n1,n2);
    ii=sf_floatalloc4(n1,n2,2*nh1+1,2*nh2+1);

    if(verb) fprintf(stderr," n3   h3\n");
    if(verb) fprintf(stderr,"%4d %3d \n",n3-1,2*nh3);

    for(ih3=-nh3; ih3<nh3+1; ih3++) { lo3=SF_ABS(ih3); hi3=n3-SF_ABS(ih3);

	/* seek in input */
	sf_seek(Fs,(lo3+ih3)*n1*n2*sizeof(float),SEEK_SET);
	sf_seek(Fr,(lo3-ih3)*n1*n2*sizeof(float),SEEK_SET);
	
	/* zero output */
	for(        ih2=-nh2; ih2<nh2+1; ih2++) { jh2=nh2+ih2;
	    for(    ih1=-nh1; ih1<nh1+1; ih1++) { jh1=nh1+ih1;
		for(    i2=0; i2<n2; i2++) { 
		    for(i1=0; i1<n1; i1++) { 
			ii[jh2][jh1][i2][i1] = 0;
		    } // n1
		} // n2
	    } // nh1
	} // nh2

	for(i3=lo3; i3<hi3; i3++) {
	    /* read input */
	    sf_floatread(us[0],n1*n2,Fs);
	    sf_floatread(ur[0],n1*n2,Fr);

	    lo2=nh2; hi2=n2-nh2;
	    lo1=nh1; hi1=n1-nh1;
	    
	    if(verb) fprintf(stderr,"%4d %3d",i3,nh3+ih3);

/*#ifdef _OPENMP*/
/*#pragma omp parallel for schedule(dynamic,ompchunk) private(i1,i2,ih1,ih2,jh1,jh2,j1,j2,k1,k2) shared(lo1,lo2,hi1,hi2,nh1,nh2,ii,us,ur)*/
/*#endif	  */
/*	    for(    ih2=0;  ih2<2*nh2+1; ih2++) { j2=i2-(ih2-nh2); k2=i2+(ih2-nh2);*/
/*		for(ih1=0;  ih1<2*nh1+1; ih1++) { j1=i1-(ih1-nh1); k1=i1+(ih1-nh1);*/
/*		    for(             i2=lo2; i2<hi2;      i2++) { */
/*			for(         i1=lo1; i1<hi1;      i1++) {*/
/*			    ii[ih2][ih1][i2][i1] += us[j2][j1] */
/*				*                   ur[k2][k1];*/
/*			}*/
/*		    }*/
/*		}*/
/*	    }*/
/**/

	    if(SF_ABS(nh2)>0) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,ompchunk) private(ih2,ih1,lo2,lo1,hi2,hi1,jh2,jh1,i2,i1,j2,j1,k2,k1) shared(nh2,nh1,ii,us,ur)
#endif		
		for(        ih2=-nh2; ih2<nh2+1; ih2++) { lo2=SF_ABS(ih2); hi2=n2-lo2; jh2=nh2+ih2;
		    for(    ih1=-nh1; ih1<nh1+1; ih1++) { lo1=SF_ABS(ih1); hi1=n1-lo1; jh1=nh1+ih1;  
			for(    i2=lo2; i2<hi2; i2++) { j2=i2-ih2; k2=i2+ih2;
			    for(i1=lo1; i1<hi1; i1++) { j1=i1-ih1; k1=i1+ih1;
				ii[jh2][jh1][i2][i1] += us[j2][j1] 
				    *                   ur[k2][k1];
			    } // n1
			} // n2  
		    } // nh1
		} // nh2
	    } else {
		
		
	    }

	    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
	} // n3
	
	/* write output */
	sf_floatwrite(ii[0][0][0],n1*n2*(2*nh1+1)*(2*nh2+1),Fi);    
    } // nh3
    if(verb) fprintf(stderr,"\n");
	
    exit (0);
}
