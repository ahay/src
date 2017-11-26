/* Kirchhoff zero-offset modeling/migration antialiased by parameterization */
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

#include <math.h>

#include <rsf.h>

#include "eigen1.h"


void solveSymmetric22(float a1, float a2, float a3, float *vv00, float *vv01, float *vv10, float *vv11, float *d0, float *d1)
/*< Apply >*/
{
  // Copy matrix to local variables.
  float a00, a01, a11;
  float v00, v10, v01, v11;
  float tiny;
  float c,r,s,t,u,vpr,vqr;
  float dt, vt;

  a00 = a1;
  a01 = a2;
  a11 = a3;

  // Initial eigenvectors. 
  v00 = 1.0f; v01 = 0.0f;
  v10 = 0.0f; v11 = 1.0f;

  // If off-diagonal element is non-zero, zero it with a Jacobi rotation.
  if (a01!=0.0f) {
      
    tiny = 0.1f*sqrt(FLT_EPSILON); // avoid overflow in r*r below

    u = a11-a00;
    
  if (fabs(a01)<tiny*fabs(u)) {
        t = a01/u;

      } else {

        r = 0.5f*u/a01;

        t = (r>=0.0f)?1.0f/(r+sqrt(1.0f+r*r)):1.0f/(r-sqrt(1.0f+r*r));
      }

      c = 1.0f/sqrt(1.0f+t*t);
      s = t*c;
      u = s/(1.0f+c);
      r = t*a01;
      a00 -= r;
      a11 += r;
      //a01 = 0.0f;
      vpr = v00;
      vqr = v10;
      v00 = vpr-s*(vqr+vpr*u);
      v10 = vqr+s*(vpr-vqr*u);
      vpr = v01;
      vqr = v11;
      v01 = vpr-s*(vqr+vpr*u);
      v11 = vqr+s*(vpr-vqr*u);

    }

    // Copy eigenvalues and eigenvectors to output arrays.
    *d0 = a00;
    *d1 = a11;
    *vv00 = v00;  *vv01 = v01;
    *vv10 = v10;  *vv11 = v11;

    // Sort eigenvalues (and eigenvectors) in descending order.
    if (*d0<*d1) {
      dt = *d1;
      *d1 = *d0;
      *d0 = dt;
      vt = *vv10;
      *vv10 = *vv00;
      *vv00 = vt;

      vt = *vv11;
      *vv11 = *vv01;
      *vv01 = vt;

    }

    //sf_warning("inside ux=%f uy=%f vx=%f vy=%f",*vv00,*vv01,*vv10,*vv11);

}
