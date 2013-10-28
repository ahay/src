/* Prediction of the LPF components on the fully frequency
band  despite possible limitation on the full bandiwitch
based on Seismic trace interpolation in the F-X domain
S.Spitz Geophysics 56, 785(1991) (Appendix B)
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

#include "cicaf1.h"
#include "predictionP.h"

static int *nq;
static int *k;
static int max_nq;
static int NN;
static int order;
static int ny;
static bool verb;
static sf_complex **Q;
static sf_complex **Qconj;
//static 
static sf_complex *temp_in;
static sf_complex *temp_out;
static sf_complex *diag;
static sf_complex *xx;
static sf_complex *yy;
static sf_complex *q; 


int factorial(int n) {
/*<compute n!>*/
	if ( n <= 1 )
      return 1;
   else
      return n * factorial( n - 1 ); 
}

int max ( int* a, int n) {
/*<compute max val in an int array>*/
        float max_value;
        int i;
		
        max_value = a[0];

        for (i = 1; i < n; i++)
		{
		/*sf_warning("a[i]=%d",(int) a[i]);*/
                if (a[i] > max_value)
				{
                        max_value = a[i];
                }
        }

        return max_value;
		
}

void Pprediction_init(int order1 ,
					  int *k1,/*band frequencies in samples*/
					  int ny1, /*Nyquist frequency sample*/
					  bool verbosity)
{
/*< initialize >*/
	int i1,i2;
	verb=verbosity;
	order=order1;
	k=k1;
	ny=ny1;
	NN = k[1]-k[0]+1;
	
	nq 	= sf_intalloc(order);
	for (i2=0; i2<order;i2++) {
		nq[i2]=(factorial(order) / factorial(order-(i2+1)) / factorial((i2+1)));
		
	} 	

	
	max_nq = (int)max(nq,order);

	Q  		= sf_complexalloc2(max_nq,order);
    Qconj   = sf_complexalloc2(max_nq,order); // Qconj = conjf(Q) 
	for (i1=0;i1<order;i1++)
		for (i2=0;i2<max_nq;i2++) {
		Q[i1][i2]=sf_cmplx(0.0,0.0);
		Qconj[i1][i2]=sf_cmplx(0.0,0.0);
		}
	//initialization
    //sf_warning("###### QUI nq[0]=%d,order=%d,max=%d",nq[0],order,max(  nq, order) );
	q 		= sf_complexalloc(max_nq);
	yy 		= sf_complexalloc(NN-1);
	xx 		= sf_complexalloc(NN-1);
	temp_in  = sf_complexalloc(max_nq);
	temp_out = sf_complexalloc(order);
	diag    = sf_complexalloc(order);
	//lpf     = sf_complexalloc2(order,max_nq);
}

void Pprediction_apply(sf_complex **LPF)
/*<apply >*/
{
    int i1,i2,i3,i4;
	
	sf_complex *up,*dw,*tofree;
	sf_complex **lpf;
	
	/* we apply appendix B rules in Spitz paper to find 
	out of band LP filter components
	First we must create Q prediction filters which have components equal to 
	binomial coefficient(L,j) */

	for (i2=0; i2<order;i2++) {
    
    /* Choosing PL vectors in frequency axis*/
    	for (i1=k[0];i1<k[1];i1++) {
			/*temp[i1-k[0]] = LPF[i1][i2];*/
			xx[i1-k[0]] = LPF[i1][i2];
			yy[i1-k[0]] = LPF[i1+1][i2];
		}
		
		/*for (i1=0; i1<NN-1;i1++){	
			xx[i1]=temp[i1];		
			yy[i1]=(temp[i1+1]);				
   		}*/
	
		/*sf_warning("data vector yy = ###");
		for (i1=0; i1<NN-1;i1++) cprint(yy[i1]);
		sf_warning("### xx vector = ###");
		for (i1=0; i1<NN-1;i1++) cprint(xx[i1]);*/

    	cicaf1_init(NN-1, xx, 1);  
		sf_csolver (cicaf1_lop, sf_ccgstep, nq[i2], NN-1,
			q, yy, NN-1 ,"verb",verb,"end");
		/*sf_warning("#### q1");*/
		sf_ccgstep_close();
		for (i1=0; i1<nq[i2];i1++) {
			/*cprint(q[i1]);*/
			Q[i2][i1]=q[i1];			
		}
       		
		for (i1=0;i1<k[1]-k[0];i1++){	
			/*temp[i1-k[0]] = LPF[i1][i2];*/
			xx[i1]=LPF[k[1]-i1][i2];	
			yy[i1]=LPF[k[1]-i1-1][i2];	
		}

		/*for (i1=0; i1<NN-1;i1++){	
			xx[i1]=temp[NN-i1-1];	
			yy[i1]=temp[NN-i1-2];	
		}*/

		/*sf_warning("###data vector yy = ###");
		for (i1=0; i1<NN-1;i1++) cprint(yy[i1]);
		sf_warning("#### xx vector = ###");
		for (i1=0; i1<NN-1;i1++) cprint(xx[i1]);*/
		

		cicaf1_init(NN-1, xx, 1);  
		sf_csolver (cicaf1_lop, sf_ccgstep, nq[i2], NN-1,
			q, yy, NN-1 ,"verb",verb,"end");
		/*sf_warning("### q2");*/
		sf_ccgstep_close();
		
		cicaf1_close();		

		for (i1=0; i1<nq[i2];i1++) {
			/*cprint(q[i1]);*/
#ifdef SF_HAS_COMPLEX_H
			Q[i2][i1]= .5*( Q[i2][i1]+conjf(q[i1]));
			Qconj[i2][i1]= conjf(Q[i2][i1]);
#else
			Q[i2][i1]= sf_crmul( sf_cadd( Q[i2][i1],conjf(q[i1]) ) ,  .5 );
			Qconj[i2][i1]=  conjf( Q[i2][i1] );
#endif
		}
		
		/*sf_warning("####### Q #########");
		for (i1=0; i1<nq[i2];i1++) 
			cprint(Q[i2][i1]);
		sf_warning("####### Q conj  #########");
		for (i1=0; i1<nq[i2];i1++) 
			cprint(Qconj[i2][i1]);*/
	}
		
	lpf     = sf_complexalloc2(order,max_nq);
	tofree = lpf[0];
	
	memcpy(lpf[0], LPF[k[0]], order*max_nq*sizeof(sf_complex) );
	/*for (i2=0;i2<order;i2++)
		for (i1=0;i1<max_nq;i1++)
			lpf[i1][i2]=LPF[k[0]+i1][i2];*/

	sf_cmatmult_init(Qconj);
	for (i1=0;i1<k[0];i1++) {
		
		for (i2=0;i2<order;i2++) { /* matrix - matrix product */
		/*sf_warning("temp_in vector i2=%d####",i2+1);*/
			for (i3=0;i3<max_nq;i3++) {
				temp_in[i3] = lpf[i3][i2];
				/*cprint(temp_in[i3]);*/
				//computes temp_in vector
			}

			sf_cmatmult_lop (false, false, max_nq,order, 
					 temp_in,temp_out);
			diag[i2] = temp_out[i2];					
		} /* END matrix - matrix product */
		
		/*for ( i3=0;i3<max_nq-1;i3++)
			lpf[i3]=lpf[i3+1];
		lpf[i3]=up; */ /* shift rows down */
		/*sf_warning("pointer");		
		for ( i3=0;i3<max_nq-1;i3++)
			sf_warning("%xx",lpf[i3]);*/

		/*sf_warning("diag %d",i1+1);	
		for ( i3=0;i3<order;i3++)
			cprint(diag[i3]);*/

		/*shift rows up*/		
		up=lpf[0];
		for ( i3=0;i3<max_nq-1;i3++) {
			dw=lpf[i3+1];	
			lpf[i3+1]=up;
			up=dw;
			if (i3==max_nq-2)
				lpf[0]=up;
		}
		for( i4=0;i4<order;i4++) {
		lpf[0][i4]=diag[i4];
		LPF[k[0]-i1-1][i4]=diag[i4];

		}		
	}
	
	free(tofree);
	lpf     = sf_complexalloc2(order,max_nq);
	tofree=lpf[0];
	
	for (i3=0;i3<max_nq;i3++) {    
		memcpy(lpf[i3], LPF[k[1]-i3], order*sizeof(sf_complex) );
		//sf_warning("ny-14 = %d !! k[1]= %d",ny-14,k[1]);
	    //for (i4=0;i4<order;i4++)
			//cprint(LPF[k[1]][i4]);
			//cprint(lpf[1][i4]);
	}
	sf_cmatmult_init(Q);
	//sf_warning("ny=%d,k[1]=%d,ny-k[1]-2=%d",ny,k[1],ny-k[1]-2);
	
	for (i1=0;i1<ny-k[1]-1;i1++) {
		
		for (i2=0;i2<order;i2++) { /* matrix - matrix product */
		//sf_warning("temp_in vector i2=%d####",i2+1);
			for (i3=0;i3<max_nq;i3++) {
				temp_in[i3] = lpf[i3][i2];
				//cprint(temp_in[i3]);
				//computes temp_in vector
			}

			sf_cmatmult_lop (false, false, max_nq,order, 
					 temp_in,temp_out);
			diag[i2] = temp_out[i2];					
		} /* END matrix - matrix product */
		
		/*for ( i3=0;i3<max_nq-1;i3++)
			lpf[i3]=lpf[i3+1];
		lpf[i3]=up; */ /* shift rows down */
		/*sf_warning("pointer");		
		for ( i3=0;i3<max_nq-1;i3++)
			sf_warning("%xx",lpf[i3]);*/

		/*sf_warning("diag %d",i1+1);	
		for ( i3=0;i3<order;i3++)
			cprint(diag[i3]);*/

		/*shift rows up*/		
		up=lpf[0];
		for ( i3=0;i3<max_nq-1;i3++) {
			dw=lpf[i3+1];	
			lpf[i3+1]=up;
			up=dw;
			if (i3==max_nq-2)
				lpf[0]=up;
		}
		//sf_warning("LPF %d = ",k[1]+i1+1);
		for( i4=0;i4<order;i4++) {
		lpf[0][i4]=diag[i4];
		LPF[k[1]+i1+1][i4]=diag[i4];
		//cprint(LPF[k[1]+i1][i4]);
		}		
	}
	//free(tofree);


	//exit(0);


}

void Pprediction_close(){
/*< close >*/
	free(nq);  		
	free(q); 
	free(*Q);free(Q);
	free(*Qconj);free(Qconj);
	free(temp_in);	
	free(temp_out);	
	free(yy);
	free(xx);
	free(diag);

}











