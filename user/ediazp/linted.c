/* Simple 1-D linear interpolation */
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

#include "linted.h"

static float o1, d1;
static const float *coord;

void lint1_init (float o1_in, float d1_in /* regular axis sampling */, 
		 const float *coord_in    /* irregular coordinates */)
/*< initialize >*/
{
    o1 = o1_in;
    d1 = d1_in;
    coord = coord_in;
}

void lint1_lop  (bool adj, bool add, int nm, int nd, float *mm, float *dd)
/*< linear interpolation >*/
{
    int im, id;
    float f, fx, gx;

    sf_adjnull(adj,add,nm,nd,mm,dd);

    for (id=0; id < nd; id++) {
        f = (coord[id]-o1)/d1;     
	im=floorf(f);
        if (0 <= im && im < nm-1) {   
	    fx=f-im;   
	    gx=1.-fx;

	    if(adj) {
		mm[im]   +=  gx * dd[id];
		mm[im+1] +=  fx * dd[id];
	    } else {
		dd[id]   +=  gx * mm[im]  +  fx * mm[im+1];
	    }
	}
    }
}

void lint1_interp (int nin, int nout, float *yin, float *yout, bool extend)
/*< linear interpolation >*/
{
    int k,i,j,jj;
    float *m,xout;

    m=sf_floatalloc(nin-1);
    xout=o1;


    for (j=1; j<nin ;j++){
	    m[j-1]=(yin[j]-yin[j-1])/(coord[j]-coord[j-1]) ; 
    }
    jj=1;
    for (i=0; i<nout; i++){
        k=0;
        for (j=jj; j<nin; j++) {
            if(xout >= coord[j-1] && xout < coord[j]){
               jj=j;
               break; 
            }
        }

        k=jj;	

	    if(xout > coord[nin-1]) { 
    		k=-1;
	    }
        if(xout < coord[0]){
            k=0;
        }	

        // fprintf(stderr,"%d%10.4f%10.4f\n",k,xout,coord[0]);
        yout[i]=-999.0;
    	if(k<0) { 
            if(extend){
        		yout[i]=yin[nin-1];
            }
            
	    }else if (k==0) {
            if(extend){
        		yout[i]=yin[0] ;
            }
        }else {
        //    fprintf(stderr,"%d%10.4f%10.4f%10.4f\n",k,coord[k-1],xout,coord[k]);
		    yout[i]=yin[k-1] + m[k-1]*(xout-coord[k-1]);
        }
        xout+= d1;
    }
}



