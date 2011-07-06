/* Local cubic interpolation. */

/*************************************************************************

Copyright Rice University, 2008.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/
#include <trip/base.h>

#include "lci.h"

void lci(int nx,       /* output length                     */
	 int ny,       /* input length                      */
	 float dy,     /* input step,                       */
	 float oy,     /* input origin                      */
	 float * yx,   /* abscissa map (length = nx)        */
	 float * fy,   /* input data, fcn y (length = ny)   */
	 float * fx    /* output data, fcn x (length = nx)  */
	 ) 
/*< local cubic interpolation >*/
{
  
  float denom2 = 1.0/(2.0e+00*dy*dy*dy);
  float denom6 = 1.0/(6.0e+00*dy*dy*dy);
  float maxy = oy+ny*dy;
  
  float y  = 0.0e+00;
  float ym1= 0.0e+00;
  float y0 = 0.0e+00;
  float y1 = 0.0e+00;
  float y2 = 0.0e+00;

  float p0 = 0.0e+00;
  float p1 = 0.0e+00;

  double oops;
  double duh;

  int i;
  int j = 0;
  int k = 0;

  for (i=0; i<nx; i++) {
    
    /* cannot use in c90 - requires c99
    y = fminf(maxy,fmaxf(yx[i],oy));
    */
    oops=oy;
    duh=yx[i];
    duh=fmax(duh,oops);
    oops=maxy;
    y=fmin(oops,duh);

    j = iwave_min(ny-2,iwave_max((int)((y-oy)/dy),0));
    k = iwave_min(iwave_max(0,j-1),ny-4);

    ym1 = y - k*dy;
    y0  = ym1 - dy;
    y1  = y0  - dy;
    y2  = y1  - dy;
    
    p0 = denom6*y0*y1;
    p1 = denom2*ym1*y2;

    fx[i] =
      - p0*y2*fy[k] 
      + p1*y1*fy[k+1] 
      - p1*y0*fy[k+2] 
      + p0*ym1*fy[k+3];
  }
}

