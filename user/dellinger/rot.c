/* Subroutines for VR programs */
/*
  Copyright (C) 1991 The Board of Trustees of Stanford University
  
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

void rot (double angle, double *x, double *y)
/*< rotation >*/
{
    double          xx, yy;

    angle *= -SF_PI / 180.;

    xx = *x * cos (angle) - *y * sin (angle);
    yy = *x * sin (angle) + *y * cos (angle);

    *x = xx;
    *y = yy;
}

double dot (double *x, double *y)
/*< dot product >*/
{
    double          out;

    out = x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    return out;
}
