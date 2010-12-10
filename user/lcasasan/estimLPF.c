/* Spitz linear filtering for trace interpolation 
   compute linear prediction filter taps at each frequency
*/
/*
  Copyright (C) 2010 Politecnico di Milano
  
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
#include "fxPEF.h"

#include "estimLPF.h"

static int nb; // PEF order, number of filter taps
static int L;
static int Ntot;
static int Nint;
static int nh;
static bool verb;
static sf_complex* xx_up;
static sf_complex* xx_dw;
static sf_complex* y;

void LPFestim_init(int nh1, /*window in space direction*/
					int order,
					int ntraces,
					bool verbosity)
/*<spitz filtering initialization>*/
{
	nh=nh1;	
    //hanning_win (nWIN, 1, win);

	nb=order;
	L = ntraces+1;
	verb=verbosity;
 	
	Ntot=L*(nh-1)+1; 
	Nint=Ntot-nh;
	/* memory ALLOCATIONS*/	
	
	xx_up 	=	sf_complexalloc( nh - 1 );
	xx_dw 	=	sf_complexalloc( nh - 1 );
	y 		= 	sf_complexalloc( 2 * (nh - 1 ) );
}

void LPFestim_apply(sf_complex* fftL /*in*/ ,
					sf_complex* bb   /*out*/ )
/*<spitz filtering>*/
{
	int lag=1;
	int ih,ib;

	/* now we compute the linear PEF from the low frequencies in fftL */

	for (ih=0; ih<nh-1;ih++){	
		xx_up[ih]=fftL[ih];		
		y[ih]=(fftL[ih+1]);					
		#ifdef SF_HAS_COMPLEX_H
		xx_dw[ih]=conjf(fftL[nh-ih-1]);	
		y[ih+(nh-1)]=conjf(fftL[nh-ih-2]);	
		#else
		xx_dw[ih]=sf_conjf(fftL[nh-ih-1]);	
		y[ih+(nh-1)]=sf_conjf(fftL[nh-ih-2]);	
		#endif	
	}
	
    /*sf_warning("##### data vector #######");
	for (int ih=0;ih<2*(nh-1);ih++) cprint(yy[ih]);
	for (int ih=0;ih<2*(nh-1);ih++) cprint(y[ih]); */
	
	/*initalizing PEF lop*/
	fxPEF_init(lag,nh-1,xx_up,xx_dw);
	
	/*fxPEF_lop(adj, add, nb, 2*(nh-(nb-1)-1), bb, data);

	sf_warning("#####debug: uscita del PEF LOP %d",2*(nh-(nb-1)-1));
	if (adj) {
		for (int ib=0; ib<nb;ib++){	
		cprint (bb[ib]);
		}
	} else { 
		for (int iy=0; iy<2*(nh-(nb-1)-1);iy++){	
		cprint(data[iy]);
		}	
	}*/

	/* initialize with zero */
	for (ib=0; ib<nb;ib++) bb[ib]=sf_cmplx(0.0,0.0);
	

	sf_csolver (fxPEF_lop, sf_ccgstep, nb, 2*(nh-1),
			bb, y, 2*nb,"x0",bb,"verb",verb,"end");
	
	sf_ccgstep_close();
	fxPEF_close();
}

void LPFestim_close()
/*<free allocated memory>*/
{

	free(xx_up);
	free(xx_dw); 	
	free(y); 

}

