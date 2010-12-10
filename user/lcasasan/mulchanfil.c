/* this function will creates the linear operator necessary for fx
interpolation with multichannel filtering (backward and forward PEF)*/

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
#include "cicai1.h"
#include "cicai1_dw.h"

#include "mulchanfil.h"
#include "carray.h"

static int nb;
static int lag;
static sf_complex* bb_up;
static sf_complex* bb_dw;
static	sf_complex* dat1;
static	sf_complex* dat2;
static  int nd;

void MCHF_init(int nb1           /* filter size */, 
		         sf_complex* bb_up1  /* filter taps [nb] */,
				 sf_complex* bb_dw1  /* filter taps [nb] */,
			     int lag1  /* filter lag (lag=1 is causal) */,
				 int nd1)
/*< initialize >*/
{

    nb  = nb1;
    lag = lag1;
	bb_up = bb_up1;
	bb_dw = bb_dw1;
	nd=nd1;
	dat1   =  sf_complexalloc(nd);
	dat2   =  sf_complexalloc(nd);	

}

void MCHF_lop(bool adj, bool add, int nx, int ndata, sf_complex* xx, sf_complex* data)
/*<linear operator>*/
{



	int nd1,nd2,id;
		
    nd1=nd;
	nd2=nd;

	

	for (id=0;id<nd1;id++) {
		dat1[id]=data[id];
		dat2[id]=data[id+nd1];
	}

	sf_cadjnull(adj, add, nb, ndata, xx, data);

	cicai1_init(nb, bb_up,lag); 
	cicai1_dw_init(nb, bb_dw,lag); 

	carray( cicai1_lop, cicai1_dw_lop, adj, add, nx, nd1, nd2, xx, dat1, dat2);

	if (!adj){
		for (id=0;id<nd1;id++) {
			data[id]=dat1[id];
			data[id+nd1]=dat2[id];
		}	
	}
}


void MCHF_close()
/*<free allocated space>*/
{
	free(dat1);
	free(dat2);
}
