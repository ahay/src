/* OpenMP time-domain cross-correlation on axes 1,2,3 */

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

    sf_file Fs,Fr,Fi;    /* I/O files */
    sf_axis a1,a2,a3,aa; /* cube axes */
    int     i1,i2,i3;

    int nbuf,ibuf;

    float ***us=NULL;
    float ***ur=NULL;
    float  **ii=NULL;

    float scale;

    int ompnth=0;
#ifdef _OPENMP
    int ompath=1; 
#endif

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

#ifdef _OPENMP
    if(! sf_getint("ompnth",  &ompnth))     ompnth=0;
    /* OpenMP available threads */
    
#pragma omp parallel
    ompath=omp_get_num_threads();
    if(ompnth<1) ompnth=ompath;
    omp_set_num_threads(ompnth);
    sf_warning("using %d threads of a total of %d",ompnth,ompath);
#endif
    
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getint ("axis",&axis)) axis=2;     /* stack axis */
    if(! sf_getint ("nbuf",&nbuf)) nbuf=1;     /* buffer size */

    Fs = sf_input ("in" );
    Fr = sf_input ("uu" );
    Fi = sf_output("out");

    /* read axes */
    a1=sf_iaxa(Fs,1); if(verb) sf_raxa(a1);
    a2=sf_iaxa(Fs,2); if(verb) sf_raxa(a2);
    a3=sf_iaxa(Fs,3); if(verb) sf_raxa(a3);

    aa=sf_maxa(1,0,1); 
    sf_setlabel(aa,""); 
    sf_setunit (aa,""); 

    nbuf = SF_MIN(nbuf,sf_n(a3));

    switch(axis) {
	case 3:
	    sf_oaxa(Fi,a1,1);
	    sf_oaxa(Fi,a2,2);
	    sf_oaxa(Fi,aa,3);
	    ii=sf_floatalloc2(sf_n(a1),sf_n(a2)); 
	    scale = 1./sf_n(a3);
	    break;
	case 2:
	    sf_oaxa(Fi,a1,1);
	    sf_oaxa(Fi,a3,2);
	    sf_oaxa(Fi,aa,3);
	    ii=sf_floatalloc2(sf_n(a1),nbuf); 
	    scale = 1./sf_n(a2);
	    break;
	case 1:
	default:
	    sf_oaxa(Fi,a2,1);
	    sf_oaxa(Fi,a3,2);
	    sf_oaxa(Fi,aa,3);
	    ii=sf_floatalloc2(sf_n(a2),nbuf); 
	    scale = 1./sf_n(a1);
	    break;
    }

    us = sf_floatalloc3(sf_n(a1),sf_n(a2),nbuf);
    ur = sf_floatalloc3(sf_n(a1),sf_n(a2),nbuf);

    i3=sf_n(a3);
    for (; i3 > 0; i3 -= nbuf) {
	if (nbuf > i3) nbuf=i3;
	if(verb) sf_warning("nsiz=%ld nbuf=%ld",i3,nbuf);

	sf_floatread(us[0][0],sf_n(a1)*sf_n(a2)*nbuf,Fs);
	sf_floatread(ur[0][0],sf_n(a1)*sf_n(a2)*nbuf,Fr);

	switch(axis) {
	    case 3:
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2)			   \
    shared( nbuf,a1,a2,ii,us,ur,scale)
#endif
		for(ibuf=0; ibuf<nbuf; ibuf++) {
		    for    (i2=0; i2<sf_n(a2); i2++) {
			for(i1=0; i1<sf_n(a1); i1++) {
			    ii[i2][i1] += us[ibuf][i2][i1]*ur[ibuf][i2][i1];
			}
		    }
		} // ibuf
		
		break;
		
	    case 2:
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2)			   \
    shared( nbuf,a1,a2,ii,us,ur,scale)
#endif
		for(ibuf=0; ibuf<nbuf; ibuf++) {
		    for(i1=0; i1<sf_n(a1); i1++) {
			ii[ibuf][i1]=0;
		    }
		    
		    for    (i2=0; i2<sf_n(a2); i2++) {
			for(i1=0; i1<sf_n(a1); i1++) {
			    ii[ibuf][i1] += us[ibuf][i2][i1]*ur[ibuf][i2][i1];
			}
		    }
		    
		    for(i1=0; i1<sf_n(a1); i1++) {
			ii[ibuf][i1] *= scale;
		    }
		} // ibuf
				
		sf_floatwrite(ii[0],sf_n(a1)*nbuf,Fi);
		break;
		
	    case 1:
	    default:
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(ibuf,i1,i2)			   \
    shared( nbuf,a1,a2,ii,us,ur,scale)
#endif
		for(ibuf=0; ibuf<nbuf; ibuf++) {
		    for(i2=0; i2<sf_n(a2); i2++) {
			ii[ibuf][i2]=0;
		    }
		    
		    for    (i2=0; i2<sf_n(a2); i2++) {
			for(i1=0; i1<sf_n(a1); i1++) {
			    ii[ibuf][i2] += us[ibuf][i2][i1]*ur[ibuf][i2][i1];
			}
		    }

		    for(i2=0; i2<sf_n(a2); i2++) {
			ii[ibuf][i2] *= scale;
		    }
		} // ibuf

		sf_floatwrite(ii[0],sf_n(a2)*nbuf,Fi);   
		break;
		
	} /* n3 */

    }
    if(verb) fprintf(stderr,"\n");    
    
    if(axis==3) {
	for    (i2=0; i2<sf_n(a2); i2++) {
	    for(i1=0; i1<sf_n(a1); i1++) {
		ii[i2][i1] *=scale;
	    }
	}
	sf_floatwrite(ii[0],sf_n(a2)*sf_n(a1),Fi);
    }
    
    /*------------------------------------------------------------*/
    ;           free(*ii); free(ii);
    free(**us); free(*us); free(us);
    free(**ur); free(*ur); free(ur);
    /*------------------------------------------------------------*/
    
    exit (0);
}
