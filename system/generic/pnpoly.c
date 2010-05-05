/* Point in polygon test. */

/* Copyright (c) 1970-2003, Wm. Randolph Franklin

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
   2. Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution.
   3. The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
*/

/* Modified for inclusion in Madagascar */
#include <rsf.h>

int pnpoly(int nvert /* number of vertices */, 
	   float **vert /* [nvert][2] vertices */, 
	   float testx, float testy /* test point */)
/*< test if point is inside a polygon >*/
{
    int i, j, c=0;

    for (i = 0, j = nvert-1; i < nvert; j = i++) {
	if ( ((vert[i][1]>testy) != (vert[j][1]>testy)) &&
	     (testx < (vert[j][0]-vert[i][0]) * (testy-vert[i][1]) / 
	      (vert[j][1]-vert[i][1]) + vert[i][0]) )
	    c = !c;
    }

    return c;
}
