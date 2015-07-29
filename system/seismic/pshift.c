/* Common phase-shift computation. */ 
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

#include "pshift.h"

static bool depth;
static char rule;

void pshift_init(bool depth1 /* depth (or time) */,
				 char rule1 /* interpolation rule */)
/*< Initialize >*/
{
    depth = depth1;
    rule = rule1;
}

sf_complex pshift(sf_complex w2, float k2, float v1, float v2, float vz, float n)
/*< phase shift for different rules >*/
{
    sf_complex cshift, cshift1, cshift2, y;
    float x;
	
#ifdef SF_HAS_COMPLEX_H
    w2 *= w2;
#else
    w2 = sf_cmul(w2,w2);
#endif
	
    switch (rule) {
		case 'a': /* VTI anisotropy */
#ifdef SF_HAS_COMPLEX_H		       
			if (depth) {
				w2 = w2 * vz * (1. + k2 / (w2 * v1 + 2.*n*k2));
			} else {
				w2 = w2 * (1. + v1 * k2 / (w2  + 2.*n*k2 * v1));
			}
#else
			if (depth) {
				w2 = sf_cmul(sf_crmul(w2,vz),
							 sf_cadd(sf_cmplx(1.,0.),
									 sf_cdiv(sf_cmplx(k2,0.),
											 sf_cadd(sf_crmul(w2,v1),
													 sf_cmplx(2.*n*k2,0.)))));
			} else {
				w2 = sf_cmul(w2,
							 sf_cadd(sf_cmplx(1.,0.),
									 sf_cdiv(sf_cmplx(v1*k2,0.),
											 sf_cadd(w2,
													 sf_cmplx(2.*n*k2 * v1,0.)))));
			}
#endif
			cshift = csqrtf(w2);
			break;
		case 's': /* simple */			
#ifdef SF_HAS_COMPLEX_H
			if (depth) {
				w2 = w2 * v1 + k2;
			} else {
				w2 = w2 + v1 * k2;
			}
#else
			if (depth) {
				w2.r = w2.r * v1 + k2;
				w2.i = w2.i * v1;
			} else {
				w2.r += v1 * k2;
			}
#endif
			cshift = csqrtf(w2);
			break;
		case 'm': /* midpoint */
#ifdef SF_HAS_COMPLEX_H
			if (depth) {
				w2 = 0.5 * w2 * (v1+v2) + k2;
			} else {
				w2 = w2 + 0.5 * (v1+v2) * k2;
			}
#else
			if (depth) {
				w2.r = 0.5 * w2.r * (v1+v2) + k2;
				w2.i = 0.5 * w2.i * (v1+v2);
			} else {
				w2.r += 0.5 * (v1+v2) * k2;
			}
#endif
			cshift = csqrtf(w2);
			break;
		case 'l': /* linear slowth */ 
		default:
#ifdef SF_HAS_COMPLEX_H
			if (depth) {
				cshift1 = csqrtf(w2 * v1 + k2);
				cshift2 = csqrtf(w2 * v2 + k2);
				
				cshift = (cshift1 + cshift2 - 
						  1./(1./cshift1+1./cshift2))/1.5;
			} else {
				cshift1 = csqrtf(w2 + v1 * k2);
				cshift2 = csqrtf(w2 + v2 * k2);
				
				x = 1./(1.+v2/v1);
				cshift = cshift1 + x*(cshift2-cshift1);
				y = x*cshift2/cshift;
				cshift *= (1.-y*(1.-y))/(1.-x*(1.-x));
			}
#else
			if (depth) {
				cshift1 = csqrtf(sf_cadd(sf_crmul(w2,v1),
										 sf_cmplx(k2,0.)));
				cshift2 = csqrtf(sf_cadd(sf_crmul(w2,v2),
										 sf_cmplx(k2,0.)));
				
				cshift = sf_crmul(
								  sf_csub(sf_cadd(cshift1,cshift2), 
										  sf_cdiv(sf_cmplx(1.,0.),
												  sf_cadd(sf_cdiv(sf_cmplx(1.,0.),
																  cshift1),
														  sf_cdiv(sf_cmplx(1.,0.),
																  cshift2)))),1.5);
			} else {
				cshift1 = csqrtf(sf_cadd(w2,sf_cmplx(v1 * k2,0.)));
				cshift2 = csqrtf(sf_cadd(w2,sf_cmplx(v2 * k2,0.)));
				
				x = 1./(1.+v2/v1);
				cshift = sf_cadd(cshift1,
								 sf_crmul(sf_csub(cshift2,cshift1),x));
				y = sf_crmul(sf_cdiv(cshift2,cshift),x);
				cshift = sf_cmul(cshift,
								 sf_crmul(
										  sf_csub(
												  sf_cmplx(1.,0.),
												  sf_cmul(y,
														  sf_csub(sf_cmplx(1.,0.),y))),
										  1./(1.-x*(1.-x))));
			}
#endif
			break;
    }
    return cshift;
}
