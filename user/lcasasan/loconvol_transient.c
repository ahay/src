/* Convolution in patches with a local filter */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
/*^*/
#include "helix_tcai.h"

#include "loconvol_transient.h"

static sf_filter aa;
static float *yy_pad;
static int nh2;

void loconvol_transient_init(sf_filter aa_in,int ny1)
/*< initialize with the first filter >*/
{
	int i;
    aa = aa_in;
    nh2 = (aa->nh-1)/2;
    sf_warning("\nny1=%d",ny1); fflush(stdout);

    yy_pad = sf_floatalloc(ny1+aa->nh-1);
    for (i=0;i<ny1+aa->nh-1;i++)
    	yy_pad[i]=0.0;
}

void loconvol_transient_lop(bool adj, bool add, int nx, int ny,
		  float *xx, float *yy)
/*< convolve >*/
{
    int i;
	helix_tcai_init(aa);
    aa++;
    sf_warning("\nnx=%d ny=%d nh=%d ntot=%d",nx,ny,aa->nh-1,ny+(aa->nh-1));

    helix_tcai_lop(adj, add, nx, ny+(aa->nh-1), xx, yy_pad);
    //memcpy(yy,yy_pad+nh2, ny);
    for (i=0;i<ny;i++)
    	yy[i] = yy_pad[i+nh2];

}

void loconvol_transient_close()
/*< free >*/
{
	free(yy_pad);
}

