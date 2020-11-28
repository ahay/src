/* OpenMP stack on axis 1,2 or 3 */

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

int main(int argc, char* argv[])
{
    bool verb;
    int  axis;
    bool norm;

    sf_file Fi,Fo;       /* I/O files */
    sf_axis a1,a2,a3,aa; /* cube axes */
    int     i1,i2,i3;

    int nbuf,ibuf;

    float ***dat=NULL;
    float  **stk=NULL;
    int    **fld=NULL;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
    
    /* OMP parameters */
#ifdef _OPENMP
    int ompnth;
    ompnth=omp_init();
	if(!ompnth)
		abort();
#endif

	if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getint ("axis",&axis)) axis=2;     /* stack axis */
    if(! sf_getbool("verb",&norm)) norm=true;  /* verbosity flag */
    if(! sf_getint ("nbuf",&nbuf)) nbuf=1;     /* buffer size */

    Fi = sf_input ("in" ); /*  input data */
    Fo = sf_output("out"); /* output data */

    /* read axes */
    a1=sf_iaxa(Fi,1); if(verb) sf_raxa(a1);
    a2=sf_iaxa(Fi,2); if(verb) sf_raxa(a2);
    a3=sf_iaxa(Fi,3); if(verb) sf_raxa(a3);

    aa=sf_maxa(1,0,1); 
    sf_setlabel(aa,""); 
    sf_setunit (aa,""); 

    nbuf = SF_MIN(nbuf,sf_n(a3));

    switch(axis) {
	case 3:
	    sf_oaxa(Fo,a1,1);
	    sf_oaxa(Fo,a2,2);
	    sf_oaxa(Fo,aa,3);
	    stk=sf_floatalloc2(sf_n(a1),sf_n(a2)); 
	    if(norm) {
		fld=sf_intalloc2(sf_n(a1),sf_n(a2)); 
	    }
	    for    (i2=0; i2<sf_n(a2); i2++) {
		for(i1=0; i1<sf_n(a1); i1++) {
		    stk[i2][i1] =0.0;
		    if(norm) fld[i2][i1]=0;
		}
	    }
	    break;
	case 2:
	    sf_oaxa(Fo,a1,1);
	    sf_oaxa(Fo,a3,2);
	    sf_oaxa(Fo,aa,3);
	    stk=sf_floatalloc2(sf_n(a1),nbuf); 
	    if(norm) {
		fld=sf_intalloc2(sf_n(a1),nbuf); 
	    }
	    break;
	case 1:
	default:
	    sf_oaxa(Fo,a2,1);
	    sf_oaxa(Fo,a3,2);
	    sf_oaxa(Fo,aa,3);
	    stk=sf_floatalloc2(sf_n(a2),nbuf); 
	    if(norm) {
		fld=sf_intalloc2(sf_n(a2),nbuf); 
	    }
	    break;
    }

    dat =sf_floatalloc3(sf_n(a1),sf_n(a2),nbuf);

    i3=sf_n(a3);
    for (; i3 > 0; i3 -= nbuf) {
	if (nbuf > i3) nbuf=i3;
	if(verb) sf_warning("nsiz=%ld nbuf=%ld",i3,nbuf);

	sf_floatread(dat[0][0],sf_n(a1)*sf_n(a2)*nbuf,Fi);

	switch(axis) {
	    case 3:
		if(norm) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2)			   \
    shared( nbuf,a1,a2,stk,fld,dat)
#endif
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
			for    (i2=0; i2<sf_n(a2); i2++) {
			    for(i1=0; i1<sf_n(a1); i1++) {
				stk[i2][i1] += dat[ibuf][i2][i1];
				if(0.0 != dat[ibuf][i2][i1]) fld[i2][i1]++;
			    }
			}
		    }
		} else {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2)			   \
    shared( nbuf,a1,a2,stk,dat)
#endif
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
			for    (i2=0; i2<sf_n(a2); i2++) {
			    for(i1=0; i1<sf_n(a1); i1++) {
				stk[i2][i1] += dat[ibuf][i2][i1];
			    }
			}
		    }
		}
		break;
		
	    case 2:
		if(norm) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2)			   \
    shared( nbuf,a1,a2,stk,fld,dat)
#endif
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
			for(i1=0; i1<sf_n(a1); i1++) {
			    stk[ibuf][i1]=0;
			    fld[ibuf][i1]=0;
			}
			
			for    (i2=0; i2<sf_n(a2); i2++) {
			    for(i1=0; i1<sf_n(a1); i1++) {
				stk[ibuf][i1] += dat[ibuf][i2][i1];
				if(0.0 != dat[ibuf][i2][i1]) fld[ibuf][i1]++;
			    }
			}
			
			for(i1=0; i1<sf_n(a1); i1++) {
			    stk[ibuf][i1] /= fld[ibuf][i1];
			}
			
		    } /* ibuf */
		    
		} else {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2)			   \
    shared( nbuf,a1,a2,stk,dat)
#endif
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
			for(i1=0; i1<sf_n(a1); i1++) {
			    stk[ibuf][i1]=0;
			}
			
			for    (i2=0; i2<sf_n(a2); i2++) {
			    for(i1=0; i1<sf_n(a1); i1++) {
				stk[ibuf][i1] += dat[ibuf][i2][i1];
			    }
			}
		    } /* ibuf */

		}		
		sf_floatwrite(stk[0],sf_n(a1)*nbuf,Fo);
		break;
		
	    case 1:
	    default:		

		if(norm) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2)			   \
    shared( nbuf,a1,a2,stk,fld,dat)
#endif
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
			for(i2=0; i2<sf_n(a2); i2++) {
			    stk[ibuf][i2]=0;
			    fld[ibuf][i2]=0;
			}
			
			for    (i2=0; i2<sf_n(a2); i2++) {
			    for(i1=0; i1<sf_n(a1); i1++) {
				stk[ibuf][i2] += dat[ibuf][i2][i1];
				if(0.0 != dat[ibuf][i2][i1]) fld[ibuf][i2]++;
			    }
			}
			
			for(i2=0; i2<sf_n(a2); i2++) {
			    stk[ibuf][i2] /= fld[ibuf][i2];
			}
			
		    } /* ibuf */
		} else {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2)			   \
    shared( nbuf,a1,a2,stk,dat)
#endif
		    for(ibuf=0; ibuf<nbuf; ibuf++) {
			for(i2=0; i2<sf_n(a2); i2++) {
			    stk[ibuf][i2]=0;
			}
			
			for    (i2=0; i2<sf_n(a2); i2++) {
			    for(i1=0; i1<sf_n(a1); i1++) {
				stk[ibuf][i2] += dat[ibuf][i2][i1];
			    }
			}
		    }
		}
		sf_floatwrite(stk[0],sf_n(a2)*nbuf,Fo);   
		break;
		
	} /* n3 */

    }
    if(verb) fprintf(stderr,"\n");    
    
    if(axis==3) {
	if(norm) {
	    for    (i2=0; i2<sf_n(a2); i2++) {
		for(i1=0; i1<sf_n(a1); i1++) {
		    stk[i2][i1] /= fld[i2][i1];
		}
	    }
	}
	
	sf_floatwrite(stk[0],sf_n(a2)*sf_n(a1),Fo);
    }
    
    /*------------------------------------------------------------*/
    ;            free(*stk); free(stk);
    ;            free(*fld); free(fld);
    free(**dat); free(*dat); free(dat);
    /*------------------------------------------------------------*/
    
    /*------------------------------------------------------------*/


    exit (0);
}
