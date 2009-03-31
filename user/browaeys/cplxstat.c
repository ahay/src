/* Statistical mean, variance and correlation for circular data. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include <float.h>
#include <math.h>

#include <rsf.h>

#include "cplxstat.h"

#define SIGN(a) (a > 0 ? 1 : (a < 0 ? -1 : 0))

static const float eps = 1.e-5;

void circ_mean(float *d, int n, float *v, float *t)
/*< circular mean >*/
{
    int i;
    float r,c,s;

    r = SF_PI/180.0;

    c = 0.0;
    s = 0.0;

    for (i = 0; i < n; i++) {
	c += cos(r*d[i]);
	s += sin(r*d[i]);
    }

    c /= n; 
    s /= n;

    /* variance */
    *v = 1.0 - (c*c + s*s);

    /* mean phase */
    *t = atan2(s,c);

    return;
}


void circ_corr(float *d1, float *d2, int n, float *corr, float *shift)
/*< circular correlation >*/
{
    int i;

    float v1,v2,vd,t1,t2,td;
    float m1,m2,md,rm,im;
    float *d;

    circ_mean(d1,n,&v1,&t1);
    m1 = sqrt(1.0 - v1);

    circ_mean(d2,n,&v2,&t2);
    m2 = sqrt(1.0 - v2);

    d = sf_floatalloc(n);

    for (i = 0; i < n; i++) d[i] = d2[i]-d1[i];

    circ_mean(d,n,&vd,&td);
    md = sqrt(1.0 - vd);

    rm = md*cos(td) - m1*m2*cos(t2 - t1);
    im = md*sin(td) - m1*m2*sin(t2 - t1);
            
    /* correlation */
    *corr = sqrt((rm*rm + im*im)/(v1*v2 + eps));

    /* phase shift */
    *shift = atan2(im,rm);

    return;
}
