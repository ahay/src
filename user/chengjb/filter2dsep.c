/* 2D spatial local and global filtering . */
/*
   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng
     
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


#include "_cjb.h"
#include "_fd.h"

void filter2dsep(float **p3, float **q3, float **p3c, float **q3c, float ****ex, float ****ez,
                 int nx, int nz, int nstep, int hnkx1, int hnkz1)
/*<filter2dsep: 2D spatial filtering for wave-modes separation >*/
{
       int i, j, ii, jj, k, l, kk, ll, ik, jl, ikm, jlm;

       for(i=0; i<nx ; i++)
       {
         ii=i/nstep;
         for(j=0; j<nz ; j++)
         {
          jj=j/nstep;
      	  for(k=-hnkx1; k<=hnkx1; k++)
          {  
             ik=i+k;
             ikm=ik+_m;
	     kk=k+hnkx1;
	     for(l=-hnkz1; l<=hnkz1; l++)
	     {
                jl=j+l;
                jlm=jl+_m;
	      	ll=l+hnkz1;

        	if(ik>=0 && ik<nx && jl>=0 && jl<nz)
	       	{
			p3c[i][j]+=p3[ikm][jlm]*ex[ii][jj][kk][ll];
			q3c[i][j]+=q3[ikm][jlm]*ez[ii][jj][kk][ll];
	      	}
             }
	   }
        }/*j loop */
     }/* i loop */
}

void filter2dsepglobal(float **p3, float **q3, float **p3c, float **q3c, float **xxx, float **zzz, int nx, int nz, int hnkx1, int hnkz1)
/*<filter2dsepglobal: 2D global spatial filtering for wave-modes separation >*/
{
       int i, j, k, l, kk, ll, ik, jl, ikm, jlm;

       for(i=0; i<nx ; i++)
       for(j=0; j<nz ; j++)
       {
      	  for(k=-hnkx1; k<=hnkx1; k++)
          {  
             ik=i+k;
             ikm=ik+_m;
	     kk=k+hnkx1;
	     for(l=-hnkz1; l<=hnkz1; l++)
	     {
                jl=j+l;
                jlm=jl+_m;
	      	ll=l+hnkz1;

        	if(ik>=0 && ik<nx && jl>=0 && jl<nz)
	       	{
			p3c[i][j]+=p3[ikm][jlm]*xxx[kk][ll];
			q3c[i][j]+=q3[ikm][jlm]*zzz[kk][ll];
	      	}
             }
	   }
        }/*j loop */
}
