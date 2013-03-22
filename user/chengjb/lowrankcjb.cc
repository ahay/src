/* 2-D two-components wavefield modeling using original elastic anisotropic displacement 
  wave equation in VTI media.

   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng and Wei Kang
     
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
#include <assert.h>
#include <rsf.hh>

/* low rank decomposition  */
#include "lowrankcjb.hh"

using namespace std;

void ddlowrankcjb(int (*sample)(vector<int>&, vector<int>&, DblNumMat&), int nx, int nk, double eps, int npk, 
                  int *m2, int *n2, float *ldata, float *fmid, float *rdata)
/*< lowrank symbol approximation: construct left, mid, and right matrixs  >*/
{
   int i, j, k;

   vector<int> md(nx), nd(nk);
   vector<int> lid, rid;

   for (k=0; k < nx; k++)
        md[k] = k;
   for (k=0; k < nk; k++)
        nd[k] = k;

   DblNumMat mid;
   DblNumMat mat;

   /* low rank decomposition for x-component */
   iC( ddlowrank(nx,nk,(*sample),eps,npk,lid,rid,mid) );

   int m, n;
   m=mid.m();
   n=mid.n();

   sf_warning("m=%d n=%d",m, n);

   fmid=sf_floatalloc(m*n);

   k=0;
   for (i=0; i < m; i++)
   for (j=0; j < n; j++)
   {
        fmid[k] = (float)mid(i,j);
        k++;
   }

   iC ( (*sample)(md,lid,mat) );

   ldata=sf_floatalloc(nx*m);
   k=0;
   for (i=0; i < nx; i++)
   for (j=0; j < m; j++)
   {
        ldata[k] = (float)mat(i,j);
        k++;
   }

   iC ( (*sample)(rid,nd,mat) );

   rdata=sf_floatalloc(n*nk);
   k=0;
   for (i=0; i < n; i++)
   for (j=0; j < nk; j++)
   {
        rdata[k] = (float)mat(i,j);
        k++;
   }

   *m2=m;
   *n2=n;
} 

void  reconstruct(float **w, float *ldata, float *fmid, float *rdata, int m, int n, int m2, int n2)
/*< re-construct matrix using the lowrank decomposed matrixs >*/
{
     int  im, in, im2, in2;
     float sum1, sum2;

     for(im=0;im<m;im++)
     for(in=0;in<n;in++)
     {
        sum1=0.0;
        for(im2=0;im2<m2;im2++)
        {

           sum2=0.0;
           for(in2=0;in2<n2;in2++)
              sum2 += fmid[im2*n2+in2]*rdata[in2*n+in];
           sum1 += ldata[im*m2+im2]*sum2;
        }
        w[im][in]=sum1;
     }
}
