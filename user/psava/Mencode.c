/* shot encoding with arbitrary delays */

/*
  Copyright (C) 2009 Colorado School of Mines

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

    sf_file Fi = NULL;	/* input */
    sf_file Fd = NULL;	/* delays */
    sf_file Fo = NULL;	/* output*/
    
    float     **idat = NULL;
    float    ***odat = NULL;
    float     **dels = NULL;

    sf_complex****tmp = NULL;
    sf_complex ***bak = NULL;

    /* cube axes */
    sf_axis at,ax,as,ae;
    int     it,ix,is,ie,ne;

    int ompchunk, ompith=0, ompnth=1;
#ifdef _OPENMP
    int ompath;
#endif

    ompfft3d ompfft=NULL; 	 /* FT structure*/
    ompsft3d ompsft=NULL;	 /* SF structure*/

    /*-----------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);    

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;
    /* OpenMP data chunk size */
#ifdef _OPENMP
    if(! sf_getint("ompnth",  &ompnth))     ompnth=0;
    /* OpenMP available threads */

#pragma omp parallel
    ompath=omp_get_num_threads();
    if(ompnth<1) ompnth=ompath;
    omp_set_num_threads(ompnth);
    sf_warning("using %d threads of a total of %d",ompnth,ompath);
#endif

    if(! sf_getbool("verb",&verb)) verb=false;	/* verbosity flag */
    
    Fi = sf_input ("in");
    Fd = sf_input ("del"); /* delay file */
    Fo = sf_output("out");
        
    /* input axes */
    at = sf_iaxa(Fi,1);
    ax = sf_iaxa(Fi,2);
    as = sf_iaxa(Fi,3);
    ae = sf_iaxa(Fd,2);
    if (verb){
	sf_raxa(at);
	sf_raxa(ax);
	sf_raxa(as);
    }

    /* output axes */
    sf_oaxa(Fo,at,1);
    sf_oaxa(Fo,ax,2);
    sf_oaxa(Fo,ae,3);

    /*-----------------------------------------------------------------*/
    /* allocate arrays */
    idat = sf_floatalloc2(sf_n(at),sf_n(ax));
    odat = sf_floatalloc3(sf_n(at),sf_n(ax),sf_n(ae));

    /* delays */
    dels = sf_floatalloc2(sf_n(as),sf_n(ae));
    sf_floatread (dels[0],sf_n(as)*sf_n(ae),Fd); 

    /* init FFT */
    ompfft = sf_ompfft3a1_init(sf_n(at),sf_n(ax),1,ompnth);
    /* init SFT */
    ompsft = sf_ompsft3_init(sf_n(at),0.0,sf_d(at),ompnth);

    tmp = sf_complexalloc4(sf_n(at),sf_n(ax),1,ompnth);
    bak = sf_complexalloc3(sf_n(at),sf_n(ax),1);

    ne = sf_n(ae);

    /* zero output */
    for         (ie=0;ie<ne;ie++){
	for     (ix=0;ix<sf_n(ax);ix++){
	    for (it=0;it<sf_n(at);it++){
		odat[ie][ix][it] = 0.0;
	    }
	}
    }
	
    for(is=0; is<sf_n(as); is++) { /* loop over sources */
	if(verb) fprintf(stderr,"%03d: ",is);	

	sf_floatread (idat[0],sf_n(at)*sf_n(ax),Fi); 

	for     (ix=0;ix<sf_n(ax);ix++){
	    for (it=0;it<sf_n(at);it++){
		bak[0][ix][it] = sf_cmplx(idat[ix][it],0);
	    }
	}
	sf_ompfft3a1(false,(kiss_fft_cpx***) bak,ompfft,0);

#ifdef _OPENMP
#pragma omp parallel for \
    private(ie,ix,it,ompith) \
    shared (is,ne,ax,at,tmp,bak,dels,ompsft,ompfft,odat)
#endif
	for         (ie=0;ie<ne;ie++){
#ifdef _OPENMP
            ompith=omp_get_thread_num();
#pragma omp critical
#endif
	    if(verb) fprintf(stderr,"%02d.%d ",ie,ompith);

	    for     (ix=0;ix<sf_n(ax);ix++){
		for (it=0;it<sf_n(at);it++){
		    tmp[ompith][0][ix][it] = bak[0][ix][it];
		}
	    }

	    sf_ompsft3_reset(sf_n(at),-dels[ie][is],sf_d(at),ompsft,ompith);  
	    sf_ompsft3a1(tmp[ompith],ompsft,ompfft,ompith);
	    sf_ompfft3a1( true,(kiss_fft_cpx***) tmp[ompith],ompfft,ompith);

	    for     (ix=0;ix<sf_n(ax);ix++){
		for (it=0;it<sf_n(at);it++){
		    odat[ie][ix][it] += crealf(tmp[ompith][0][ix][it]);
		}
	    }
	    
	} /* e loop */
	if(verb) fprintf(stderr,"\n");    

    } /* s loop */

    sf_floatwrite(odat[0][0],sf_n(at)*sf_n(ax)*ne,Fo); 
 
    sf_ompfft3a1_close(ompfft);
    sf_ompsft3_close  (ompsft);    
    /*------------------------------------------------------------*/


    exit(0);
}		
