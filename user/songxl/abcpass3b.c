/*3D ABC */
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

#include "abcpass3b.h"

static int nx, nz, nh, nbt, nbb, nxl, nxr, nhl, nhr;
static float ct, cb, cxl, cxr, chl, chr;
static float *wt, *wb, *wxl, *wxr, *whl, *whr;

void bd3_init(int n1,  int n2, int n3    /*model size:x*/,
             int nb1, int nb2  /*z:top, bottom*/,
             int nb3, int nb4  /*x:left, right*/,
             int nb5, int nb6  /*h:left, right*/,
             float c1, float c2 /*z:top, bottom*/,
             float c3, float c4 /*x:left, right*/,
             float c5, float c6 /*h:left, right*/)
/*< initialization >*/
{
    nz = n1;  
    nx = n2;  
    nh = n3;  
    nbt = nb1;  
    nbb = nb2;  
    nxl = nb3;  
    nxr = nb4;  
    nhl = nb5;  
    nhr = nb6;  
    ct = c1;  
    cb = c2;  
    cxl = c3;  
    cxr = c4;  
    chl = c5;  
    chr = c6;  
    if(nbt) wt =  sf_floatalloc(nbt);
    if(nbb) wb =  sf_floatalloc(nbb);
    if(nxl) wxl =  sf_floatalloc(nxl);
    if(nxr) wxr =  sf_floatalloc(nxr);
    if(nhl) whl =  sf_floatalloc(nhl);
    if(nhr) whr =  sf_floatalloc(nhr);
    abc_cal(nbt,ct,wt);
    abc_cal(nbb,cb,wb);
    abc_cal(nxl,cxl,wxl);
    abc_cal(nxr,cxr,wxr);
    abc_cal(nhl,chl,whl);
    abc_cal(nhr,chr,whr);
}
   

void bd3_close(void)
/*< free memory allocation>*/
{
    if(nbt) free(wt);
    if(nbb) free(wb);
    if(nxl) free(wxl);
    if(nxr) free(wxr);
    if(nhl) free(whl);
    if(nhr) free(whr);
}

void bd3_decay(float ***a /*3-D matrix*/) 
/*< boundary decay>*/
{
    int iz, ix, ih;
    int nzb, nxb, nhb;
    nzb = nz+nbt+nbb;
    nxb = nx+nxl+nxr;
    nhb = nh+nhl+nhr;
           for (iz=0; iz < nbt; iz++) {  
               for (ix=nxl; ix < nx+nxl; ix++) {
                   for (ih=nhl; ih < nh+nhl; ih++){
                       a[iz][ix][ih] *= wt[iz];
                   }
                   for (ih=0; ih < nhl; ih++){
                       a[iz][ix][ih] *= SF_MIN(wt[iz],whl[ih]);
                   }
                   for (ih=0; ih < nhr; ih++){
                       a[iz][ix][nhb-1-ih] *= SF_MIN(wt[iz],whr[ih]);
                   }
               }
               for (ix=0; ix < nxl; ix++) {
                   for (ih=nhl; ih < nh+nhl; ih++){
                       a[iz][ix][ih] *= SF_MIN(wt[iz],wxl[ix]);
                   }
                   for (ih=0; ih < nhl; ih++){
                       a[iz][ix][ih] *= SF_MIN(SF_MIN(wt[iz],wxl[ix]),whl[ih]);
                   }
                   for (ih=0; ih < nhr; ih++){
                       a[iz][ix][nhb-1-ih] *= SF_MIN(SF_MIN(wt[iz],wxl[ix]),whr[ih]);
                   }
                 
                }
               for (ix=0; ix < nxr; ix++) {
                   for (ih=nhl; ih < nh+nhl; ih++){
                       a[iz][nxb-1-ix][ih] *= SF_MIN(wt[iz],wxr[ix]);
                   }
                   for (ih=0; ih < nhl; ih++){
                       a[iz][nxb-1-ix][ih] *= SF_MIN(SF_MIN(wt[iz],wxr[ix]),whl[ih]);
                   }
                   for (ih=0; ih < nhr; ih++){
                       a[iz][nxb-1-ix][nhb-1-ih] *= SF_MIN(SF_MIN(wt[iz],wxr[ix]),whr[ih]);
                   }
               }
           }
           for (iz=0; iz < nbb; iz++) {  
               for (ix=nxl; ix < nx+nxl; ix++) {
                   for (ih=nhl; ih < nh+nhl; ih++){
                       a[nzb-1-iz][ix][ih] *= wb[iz];
                   }
                   for (ih=0; ih < nhl; ih++){
                       a[nzb-1-iz][ix][ih] *= SF_MIN(wb[iz],whl[ih]);
                   }
                   for (ih=0; ih < nhr; ih++){
                       a[nzb-1-iz][ix][nhb-1-ih] *= SF_MIN(wb[iz],whr[ih]);
                   }
               }
               for (ix=0; ix < nxl; ix++) {
                   for (ih=nhl; ih < nh+nhl; ih++){
                       a[nzb-1-iz][ix][ih] *= SF_MIN(wb[iz],wxl[ix]);
                   }
                   for (ih=0; ih < nhl; ih++){
                       a[nzb-1-iz][ix][ih] *= SF_MIN(SF_MIN(wb[iz],wxl[ix]),whl[ih]);
                   }
                   for (ih=0; ih < nhr; ih++){
                       a[nzb-1-iz][ix][nhb-1-ih] *= SF_MIN(SF_MIN(wb[iz],wxl[ix]),whr[ih]);
                   }
                 
                }
               for (ix=0; ix < nxr; ix++) {
                   for (ih=nhl; ih < nh+nhl; ih++){
                       a[nzb-1-iz][nxb-1-ix][ih] *= SF_MIN(wb[iz],wxr[ix]);
                   }
                   for (ih=0; ih < nhl; ih++){
                       a[nzb-1-iz][nxb-1-ix][ih] *= SF_MIN(SF_MIN(wb[iz],wxr[ix]),whl[ih]);
                   }
                   for (ih=0; ih < nhr; ih++){
                       a[nzb-1-iz][nxb-1-ix][nhb-1-ih] *= SF_MIN(SF_MIN(wb[iz],wxr[ix]),whr[ih]);
                   }
               }
           }
           for (iz=nbt; iz < nz+nbt; iz++) {  
               for (ix=nxl; ix < nxl+nx; ix++) {
                   for (ih=0; ih < nhl; ih++){
                       a[iz][ix][ih] *= whl[ih];
                   }
                   for (ih=0; ih < nhr; ih++){
                       a[iz][ix][nhb-1-ih] *= whr[ih];
                   }
                }
               for (ix=0; ix < nxl; ix++) {
                   for (ih=nhl; ih < nh+nhl; ih++){
                       a[iz][ix][ih] *= wxl[ix];
                   }
                   for (ih=0; ih < nhl; ih++){
                       a[iz][ix][ih] *= SF_MIN(wxl[ix],whl[ih]);
                   }
                   for (ih=0; ih < nhr; ih++){
                       a[iz][ix][nhb-1-ih] *= SF_MIN(wxl[ix],whr[ih]);
                   }
                 
                }
               for (ix=0; ix < nxr; ix++) {
                   for (ih=nhl; ih < nh+nhl; ih++){
                       a[iz][nxb-1-ix][ih] *= wxr[ix];
                   }
                   for (ih=0; ih < nhl; ih++){
                       a[iz][nxb-1-ix][ih] *= SF_MIN(wxr[ix],whl[ih]);
                   }
                   for (ih=0; ih < nhr; ih++){
                       a[iz][nxb-1-ix][nhb-1-ih] *= SF_MIN(wxr[ix],whr[ih]);
                   }
               }
           }
}

void abc_cal(int nb  /* absorbing layer length*/, 
             float c /* decaying parameter*/,
             float* w /* output weight[nb] */)
/*< find absorbing coefficients >*/
{
    int ib;
    if(!nb) return;
    for(ib=0; ib<nb; ib++){
       w[ib]=exp(-c*c*(nb-1-ib)*(nb-1-ib));
    }
}

