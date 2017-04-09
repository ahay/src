/* Kirchhoff integral modeling */
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

#ifndef _kirmod_h

typedef struct KTable {
    float t  /* traveltime */;
    float a  /* geometrical spreading */;
    float tx /* traveltime slope (dt/dx) */;
    float ty /* traveltime slope (dt/dy) */;
    float tn /* obliguity (dt/dn) */;
    float an /* angle from the vertical */;
    float ar /* 2.5-D factor (1/r dt/dr) */;
} *ktable;
/*^*/

#endif

void kirmod_table(char type    /* type of velocity distribution */,
		  bool twod    /* 2-D or 2.5-D/3-D */,
		  float z,
		  float x,
		  float y      /* distance between source and receiver */, 
		  float g      /* absolute gradient */,
		  float gx     /* gx+gz*zx */,
		  float gy     /* gy+gz*zy */,
		  float gz     /* gz-gx*zx */,
		  float v1     /* source velocity function */, 
		  float v2     /* receiver velocity function */,
		  float vn     /* "NMO" velocity */,
		  float n      /* "eta" parameter */,
		  float px     /* x+z*zx */,
		  float py     /* y+z*zy */,
		  float pz     /* z-x*zx */,
		  float dz     /* hypotf(1.0,zx) */,
		  ktable table /* [5] output table */)
/*< Compute traveltime attributes >*/
{
    float sigma, rad, v0, a, r, h;

    r = sqrtf(x*x+y*y+z*z)+FLT_EPSILON; /* distance */

    switch (type) {
	case 'a': /* VTI anisotropy */
	    h = z*z/(v1*v1) + x*x/((1.+2*n)*vn*vn); /* hyperbolic part */
	    table->t = sqrtf(((3.+4.*n)*h + sqrtf(h*h + 16.*n*(1.+n)*z*z*x*x/((1.+2*n)*vn*vn*v1*v1)))/(4.*(1.+n))); 
	    if (twod) {
		table->a = sqrtf(r/v1);
	    } else {
		table->a = r;
		table->ar = 1./(r*v1);
	    }
	    px /= (r*v1);
	    py /= (r*v1);
	    pz = fabsf(pz)/(r*dz);
	    table->tn = sqrtf(fabsf(1./(v1*v1)-px*px-py*py));
	    break;		     
	case 'c': /* constant velocity */
	    table->t = r/v1;
	    if (twod) {
		table->a = sqrtf(r/v1);
	    } else {
		table->a = r;
		table->ar = 1./(r*v1);
	    }
	    px /= (r*v1);
	    py /= (r*v1);
	    pz = fabsf(pz)/(r*dz);
	    table->tn = sqrtf(fabsf(1./(v1*v1)-px*px-py*py));
	    break;		    
	case 's': /* linear sloth */
	    v0 = 0.5*(v1+v2);
	    rad = v0*v0-r*r*g*g;           /* => v0^2=1/v^4 */
	    if (rad < 0.) { /* shadow zone */
	      table->t = -FLT_MAX;
	      table->a = 0.;
	      table->ar = 0.;
	      table->tn = 0.;
	      break;
	    }
	    rad = sqrtf(rad);              /* => v0=1/v^2 */
	    sigma = r*sqrtf(2./(v0+rad));  /* => r*v */
	    table->t = sigma*(2*v0+rad)/3.;       /* => r/v */
	    a = sigma*sqrtf(rad);          /* => r */
	    if (twod) {
		table->a = a/sqrtf(sigma);
	    } else {
		table->a = a;
		table->ar = 1./sigma;
	    }
	    px = (px/sigma+0.5*gx*sigma);  /* => p/(rv) */
	    py = (py/sigma+0.5*gy*sigma);  /* => p/(rv) */
	    pz = fabsf(pz+0.5*sigma*sigma*gz)/(sqrtf(v2)*sigma*dz);
	    table->tn = sqrtf(fabsf(v2-px*px-py*py));
	    break;
	case 'v': /* linear velocity */
	    v0 = sqrtf(v1*v2);                     /* => v */
	    table->t = acoshf(1.+0.5*r*r*g*g/(v0*v0))/g;  /* => r/v */
	    a = r*hypotf(1.,0.5*r*g/v0);           /* => r */
	    rad = 1./(a*v0);                       /* => 1./(r*v) */
	    if (twod) {
		table->a = a*sqrtf(rad);
	    } else {
		table->a = a;
		table->ar = rad;
	    }
	    px = (px-0.5*r*r*gx/v2)*rad;           /* => p/(r*v) */
	    py = (py-0.5*r*r*gy/v2)*rad;           /* => p/(r*v) */
	    pz = fabsf(pz*v2-0.5*r*r*gz)/(v0*a*dz);
	    table->tn = sqrtf(fabsf(1./(v2*v2)-px*px-py*py));
	    break;
	default:
	    sf_error("%s: type %c is not implemented",__FILE__,type);
	    break;
    } /* type */
    if (twod) table->ar=0.5;
    table->tx = px;                          /* traveltime slope (dt/dx) */
    table->ty = py;
    table->an = (pz >= 1.0)? 0.: -acosf(pz); /* angle from the vertical */
}

