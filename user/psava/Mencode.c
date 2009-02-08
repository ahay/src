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

#include "ftutil.h"

int main(int argc, char* argv[])
{
    bool verb;

    sf_file Fi = NULL;	/* input */
    sf_file Fd = NULL;	/* delays */
    sf_file Fo = NULL;	/* output*/
    
    float     **idat = NULL;
    float     **odat = NULL;
    float      *dels = NULL;

    sf_complex ***tmp = NULL;

    /* cube axes */
    sf_axis at,ax,as;
    int     it,ix,is;

    fft3d fft=NULL;	 /* FT structure*/
    sft3d sft=NULL;	 /* SF structure*/
    
    /*-----------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);    
    if(! sf_getbool("verb",&verb)) verb=false;	/* verbosity flag */
    
    Fi = sf_input ("in");
    Fd = sf_input ("del");	/* auxiliary input file */
    Fo = sf_output("out");
        
    /* input axes */
    at = sf_iaxa(Fi,1);
    ax = sf_iaxa(Fi,2);
    as = sf_iaxa(Fi,3);
    if (verb){
	sf_raxa(at);
	sf_raxa(ax);
	sf_raxa(as);
    }
    /* output axes */
    sf_oaxa(Fo,at,1);
    sf_oaxa(Fo,ax,2);
    sf_oaxa(Fo,as,3);

    /*-----------------------------------------------------------------*/
    /* allocate arrays */
    idat = sf_floatalloc2(sf_n(at),sf_n(ax));
    odat = sf_floatalloc2(sf_n(at),sf_n(ax));

    /* delays */
    dels = sf_floatalloc (sf_n(as));
    sf_floatread (dels,sf_n(as),Fd); 

    /* init FFT */
    fft = fft3a1_init    (sf_n(at),sf_n(ax),1); 
    tmp= sf_complexalloc3(sf_n(at),sf_n(ax),1);

    for(is=0; is<sf_n(as); is++) { /* loop over sources */
	sf_warning("is=%d del=%g",is,dels[is]);

	sf_floatread (idat[0],sf_n(at)*sf_n(ax),Fi); 
	for     (ix=0;ix<sf_n(ax);ix++){
	    for (it=0;it<sf_n(at);it++){
		tmp[0][ix][it] = sf_cmplx(idat[ix][it],0);
	    }
	}
	
	fft3a1(false,(kiss_fft_cpx***) tmp,fft);

	sft = sft3_init(sf_n(at),-dels[is],sf_d(at));
	sft3a1(tmp,sft,fft);
	sft3_close(sft);

	fft3a1( true,(kiss_fft_cpx***) tmp,fft);

	for     (ix=0;ix<sf_n(ax);ix++){
	    for (it=0;it<sf_n(at);it++){
		odat[ix][it] = crealf(tmp[0][ix][it]);
	    }
	}
	sf_floatwrite(odat[0],sf_n(at)*sf_n(ax),Fo); 
    }

    exit(0);
}		
