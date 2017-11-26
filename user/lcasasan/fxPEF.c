/* this function implements the autoregression in the fx domain
   in order to compute the spitz's PEF taps */

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
#include "cicaf1_dw.h"
#include "carray.h"

#include "fxPEF.h"

static int nx;
static int lag;
static sf_complex* xx_up;
static sf_complex* xx_dw;
static sf_complex* dat1;
static sf_complex* dat2;

void fxPEF_init(int lag1 /* lag for internal complex convolution */,
		int nx1 /* data size */,
		sf_complex* xx_up1,
		sf_complex* xx_dw1 )
/*<initialize>*/
{
    lag=lag1;
    nx=nx1;
    xx_up =  xx_up1;
    xx_dw =  xx_dw1;

    dat1   =  sf_complexalloc(nx);
    dat2   =  sf_complexalloc(nx);

}


void fxPEF_close(void)
/*<close>*/
{
    free(dat1);
    free(dat2);
}


void fxPEF_lop(bool adj, bool add, int nb, int ndata, sf_complex* bb, sf_complex* data)
/*<linear operator>*/
{
	


    int nd1,nd2,id;
		
    nd1=ndata/2;
    nd2=ndata/2;

    dat1   =  sf_complexalloc(nd1);
    dat2   =  sf_complexalloc(nd2);	

    for (id=0;id<nd1;id++) {
	dat1[id]=data[id];
	dat2[id]=data[id+nd1];
    }

    sf_cadjnull(adj, add, nb, ndata, bb, data);

    cicaf1_init(nx, xx_up,lag); 
    cicaf1_dw_init(nx, xx_dw,lag); 
    carray(  cicaf1_lop,   cicaf1_dw_lop, adj, add, nb, nd1, nd2, bb, dat1, dat2);
    /* sf_array_simple(adj, add, nb, nd1, nd2, bb, dat1, dat2); */

    if (!adj){
	for (id=0;id<nd1;id++) {
	    data[id]=dat1[id];
	    data[id+nd1]=dat2[id];
	}	
    }
}



