/* adaptive 1d filtering for coherent noise removal using */

/* Copyright (C) 2010 Politecnico di Milanoù

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
#include <rsfgee.h>

#include "filter1d_icaf.h"
#include "filter1d_tcaf.h"

#include "match1d.h"


void pad (float *d, float **D,  int nI, int nJ, int ipad, int jpad)
/*<matrix padding>*/
{
	int i,j,ind;
	for(ind=0,j=0; j < nJ; j++)
	{
		for(i=0; i < nI; i++,ind++ )
		{
			D[i+ipad][j+jpad] = d[ind];
		}
	}

	for(i=0; i < ipad; i++)
	{
		for(j=0; j < nJ+2*jpad; j++)
		{

			D[i][j] = D[2*ipad-(i+1)][j];
			D[nI+2*ipad-1-i][j] = D[nI+i][j];
		}
	}

	for(j=0; j < jpad; j++)
	{
		for(i=0; i < nI+2*ipad; i++)
		{
			D[i][j] = D[i][2*jpad-(j+1)];
			D[i][nJ+2*jpad-1-j] = D[i][nJ+j];
		}
	}
}




void adaptsub_icaf(float *input1, float *input2, float *output,  int ns,  int ntr, int order, int Twin, int Xwin, float mu, char *method, bool verb)
/*<adaptive matched-filter subtraction using internal convolution>*/
{
	float value;
	float  **Data,**Mult;
	float  *d,*m,*x;
	int i,j,ii,jj,ind,ind2,twin,xwin,k,samples;
	int oo=(Twin+1)/2 ;
	twin=(Twin-1)/2;
	xwin=(Xwin-1)/2;
	k=Twin*Xwin;
	samples=ns*ntr;


	Data=sf_floatalloc2(ntr+Xwin-1,ns+Twin-1);
	Mult=sf_floatalloc2(ntr+Xwin-1,ns+Twin-1);

	d = sf_floatalloc(k);
	m = sf_floatalloc(k);
	x = sf_floatalloc(order);

	for (i=0;i<order;i++) x[i]=0.0;

	pad (input1, Data, ns, ntr, twin, xwin);
	pad (input2, Mult, ns, ntr, twin, xwin);



	for(ind=0,j=0; j < ntr; j++)
	{
		for(i=0; i < ns; i++,ind++)
		{
			value=0.0;
			for(ind2=0,jj=0; jj < Xwin; jj++)
			{   // patching 2D
				for(ii=0; ii < Twin; ii++,ind2++)
				{
					d[ind2] =  Data[i+ii][j+jj];
					m[ind2] =  Mult[i+ii][j+jj];
				}
			}



			sf_cdstep_init();

			if ('o' == method[0] || 'O' == method[0]) {

				icaf1_init (k  /* data length */,
							m  /* data [k] */,
							1  /* filter lag (lag=1 is causal) */);

				sf_solver_prec (icaf1_lop, sf_cdstep, sf_copy_lop, order, order, k,
						x, d, 3*order, mu,"end");
			}
			else {

				filter1d_icaf_init(Twin  /* window time length */,
								    Xwin  /* window space length */,
								    m /* signal to be matched */);

				sf_solver_reg (filter1d_icaf_lop, sf_cdstep, sf_copy_lop, order, order, k,
									x, d, 3*order, mu,"end");

				filter1d_icaf_close();
			}

			sf_cdstep_close();

			for(ii=0 ; ii < order ; ii++) {
				value += x[ii] * Mult[i + oo -ii-1][j+xwin];       // matched signal
			}
			output[ind] = input1[ind] - value;       			   // subtraction

		}
		if (verb)
		sf_warning("\r\t\t\t\t\t\t\t\t\t %3.2f%% ",(float)100*(ind+1)/samples);

	}

	free(Data[0]);
	free(Data);
	free(Mult[0]);
	free(Mult);


	free(d);
	free(m);
	free(x);
}

void adaptsub_tcaf(float *input1, float *input2, float *output,  int ns,  int ntr, int order, int Twin, int Xwin, float mu, char *method, bool verb)
/*<adaptive matched-filter subtraction using transient convolution>*/
{
	float value;
	float  **Data,**Mult;
	float  *d,*m,*x;
	int nd,nm;
	int i,j,ii,jj,ind,twin,xwin,samples;

	int o  =  (order-1)/2;
	int oo =  (Twin+1)/2+(order-1)/2 ;
	twin=(Twin-1)/2;
	xwin=(Xwin-1)/2;

	samples=ns*ntr;


	Data=sf_floatalloc2(ntr+Xwin-1,ns+Twin-1);
	Mult=sf_floatalloc2(ntr+Xwin-1,ns+Twin-1);
;

	nm = Twin*Xwin;

	if ('o' == method[0] || 'O' == method[0])
		nd = (Twin*Xwin) + order-1;
	else
		nd = (Twin+order-1)*Xwin;


	d = sf_floatalloc(nd);
	m = sf_floatalloc(nm);
	x = sf_floatalloc(order);

	for (i=0;i<order;i++) x[i]=0.0;

	pad (input1, Data, ns, ntr, twin, xwin);
	pad (input2, Mult, ns, ntr, twin, xwin);



	for(ind=0,j=0; j < ntr; j++)
	{
		for(i=0; i < ns; i++,ind++)
		{
			value=0.0;
			for(jj=0; jj < Xwin; jj++)
			{   // patching 2D
				for(ii=0; ii < Twin; ii++)
				{
					if ('o' == method[0] || 'O' == method[0])
						d[ii+(jj * Twin) + o] =  Data[i+ii][j+jj];
					else
						d[ii+(jj* (Twin+2*o) ) + o] =  Data[i+ii][j+jj];

					m[ii+(jj*Twin)] =  Mult[i+ii][j+jj];
				}
			}

			//sf_cdstep_init();

			if ('o' == method[0] || 'O' == method[0]) {

				tcaf1_init (nm  /* data length */,
							m  /* data [k] */);

				sf_solver_prec (tcaf1_lop, sf_cdstep, sf_copy_lop, order, order, nd,
						x, d, 3*order, mu,"end");
			}
			else {


				filter1d_tcaf_init(Twin  /* window time length */,
								    Xwin  /* window space length */,
								    order /* filter order */,
								    m /* signal to be matched */);


				sf_solver_reg (filter1d_tcaf_lop, sf_cdstep, sf_copy_lop, order, order, nd,
									x, d, 3*order, mu,"end");



				filter1d_tcaf_close();
			}

			sf_cdstep_close();

			for(ii=0 ; ii < order ; ii++) {
				value += x[ii] * Mult[i + oo -(ii+1)][j+xwin];       // matched signal
			}
			output[ind] = input1[ind] - value;       			   // subtraction

		}
		if (verb)
		sf_warning("\r\t\t\t\t\t\t\t\t\t %3.2f%% ",(float)100*(ind+1)/samples);

	}

	free(Data[0]);
	free(Data);
	free(Mult[0]);
	free(Mult);


	free(d);
	free(m);
	free(x);
}

