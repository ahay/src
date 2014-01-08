/* Generalized p-norm thresholding operator for complex numbers
*/
/*
  Copyright (C) 2013  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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
#include <complex.h>

#include "pthresh.h"

sf_complex pthresholding(sf_complex x, float thr, float p, char* mode)
/*< p-norm thresholding operator for complex numbers >*/
{
    float a=cabsf(x);
    if (strcmp(mode,"hard") == 0)
	return (x)*(a>thr?1.:0.);/* hard thresholding*/
    if (strcmp(mode,"soft") == 0) a=1.0-thr/(a+(a==0));/* soft thresholding */
    if (strcmp(mode,"pthresh") == 0) a=1.0-powf((a+(a==0))/thr, p-2.0);
    /* generalized quasi p-norm thresholding*/
    if (strcmp(mode,"exp") == 0) a=expf(-powf((a+(a==0))/thr, p-2.0));
    /* exponential shrinkage */

    return (x)*(a>0.0?a:0.0);
}
