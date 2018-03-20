/* Spitz linear filtering for trace interpolation 
   Apply multichannel filter using linear predictor taps*/
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
#include<rsf.h>


#include "mulchanfil.h"
#include "filtMCHF.h"

static int nb; /* PEF order, number of filter taps*/
static int L;
static int Ntot;
static int Nint;
static int nh;
static bool verb;

static sf_complex* bb_for;
static sf_complex* bb_down;
static bool* known; /*mask for known data*/
static sf_complex* in;
static sf_complex* data;
static sf_complex* mchf;
static sf_complex** MCHF; /*multichannel filter Nint x nh     */
static int* ind_unknown;
static int* ind_known;
static sf_complex** eye;
static sf_complex* out_int;

void MCHFfilt_init(int nh1, /*window in space direction*/
					int order,
					int ntraces,
					bool verbosity)
/*<spitz filtering initialization>*/
{
	nh=nh1;	
	nb=order;
	L = ntraces+1;
	verb=verbosity;

	Ntot=L*(nh-1)+1; 
	Nint=Ntot-nh;
	//sf_warning("Ntot = %d, Nint=%d, nh=%d, L=%d",Ntot,Nint,nh,L);
	/* memory ALLOCATIONS*/	
	

	bb_for 	= 	sf_complexalloc( nb + 1 );
	bb_down	= 	sf_complexalloc( nb + 1 );
	known   =   sf_boolalloc( Ntot);
	in		= 	sf_complexalloc( Ntot );	
	data	=	sf_complexalloc( 2 * Ntot );
	MCHF	=	sf_complexalloc2( nh,Nint);	
	mchf	= 	sf_complexalloc( Ntot );
	ind_unknown = sf_intalloc( Nint );
	ind_known = sf_intalloc( nh )  ;
	eye		=	sf_complexalloc2( Ntot,Ntot) ;
	out_int = 	sf_complexalloc(Nint );

	
}

void MCHFfilt_apply(sf_complex *bb /*in*/,
					sf_complex *fft /*in*/, 
					sf_complex *out/*out*/)
/*<MCHF filter construction and application >*/
{
	int lag=1;
	int ih,ib,iN,i,j;

	
	bb_for[0] = sf_cmplx(-1.,0.);
	bb_down[nb] = sf_cmplx(-1.,0.);

	
	for (ib=0; ib<nb;ib++){	
		#ifdef SF_HAS_COMPLEX_H
		bb_down[(nb-1)-ib]=conjf(bb[ib]);
		#else
		bb_down[(nb-1)-ib]=sf_conjf(bb[ib]);
		#endif
		bb_for[ib+1]=bb[ib];
	}

	/*sf_warning("#####PEF filters for");
	for (ib=0; ib<nb+1;ib++) cprint(bb_for[ib]);
		
	sf_warning("#####PEF filters down");
	for (ib=0; ib<nb+1;ib++) cprint(bb_down[ib]);     */

/* now we are going to find the multichannel filter for interpolation */
	i=0;j=0;
	for (iN=0;iN<Ntot;iN++){
		if (! (iN%L) ) {/* 1:L:Ntot */
		known[iN]=true;
		ind_known[j]=iN;
		j++ ;
		} else {
        known[iN]=false;
		ind_unknown[i]=iN;
		i++;
		}	
		in[iN]=sf_cmplx(0.0,0.0);
	}

	
	
	MCHF_init(nb+1,bb_for,bb_down,lag,Ntot);
	
	for (iN=0;iN<Ntot;iN++) {
		mchf[iN]=sf_cmplx(0.0,0.0);
	}

	/* eps=0.00001; dumping parameter */

	//sf_warning("###LOP CG SOLVER for MCHF Ntot=%d",Ntot);
	ih=0;

	for (iN=0;iN<Ntot;iN++){
		
		if (!(iN%L)){

			in[iN]=sf_cmplx(-1.0,0.0);
			/*sf_warning("IN IN IN ");
			for (int iin=0;iin<Ntot;iin++) cprint(in[iin]); */
			
			MCHF_lop(false,false,Ntot,2*Ntot,in,data); 

			/* normal solver  */
			sf_csolver (MCHF_lop, sf_ccgstep,Ntot, 2*Ntot,
			mchf, data, Ntot,"verb",verb,"known",known,"end");
			/* preconditioned solver */
//			sf_csolver_prec (MCHF_lop, sf_ccgstep,sf_ccopy_lop,Ntot, Ntot, 2*Ntot,
//			mchf, data, Ntot,eps,"verb",verb,"known",known,"end");

			sf_ccgstep_close();
			in[iN]=sf_cmplx(0.0,0.0);
			
			for (i=0 ; i<Nint; i++) {			
				MCHF[i][ih]=mchf[ind_unknown[i]];
            }
			/*sf_warning("###MCHF estimated");
			sf_warning("###analized trace # %d",ih);
			for (i=0 ; i<Nint; i++) cprint(MCHF[i][ih] );     */
		    ih++;
		    
	    }
 
		
	//sf_warning("###iteration # %d",iter);
	}
    //sf_warning("total number of traces %d",Ntot);
	//sf_warning("total number of interp traces %d",Nint);
    //for (i=0;i<Nint;i++) sf_warning("ind_unknown %d ind_known %d",ind_unknown[i],ind_known[i]);    
	/*Applying MCHF on the know traces to compute the uknown traces */
	sf_cmatmult_init(MCHF);
	sf_cmatmult_lop (false, false, nh, Nint, fft, out_int );  
	/*for (i=0;i<Nint;i++)
	for (ih=0;ih<nh;ih++) 
     	out_int[i]+=MCHF[i][ih]*fft[ih];	*/

    //sf_warning("After applying MCHF");
	for (i=0;i<Nint;i++) 
		out[ind_unknown[i]]= out_int[i];
	for (ih=0;ih<nh;ih++) 
	    out[ind_known[ih]]= fft[ih];
	
	//for (i=0;i<Ntot;i++) cprint(out[i]);
	MCHF_close();
}

void MCHFfilt_close(void)
/*<free allocated space>*/
{
    free(bb_for);
    free(bb_down);
    free(known);
    free(data);
    free(in);
    free(mchf);
    free(ind_unknown);	
    free(eye[0]);
    free(out_int);
    free(MCHF[0]);
}
