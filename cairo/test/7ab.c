/* */
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

#include <cairo.h>

#include <rsf.h>
#include <rsfcairo.h>

static float top, c1, c2;

void cr_init(void) {}

void cr_show(cairo_t *cr)
{
    cr_uorig (-c1,-.5);
    cr_uclip (0.,top-3.,4.,top);
    cr_umove (cr, 0.,top-0.);  
    cr_udraw (cr, 0.,top-4.);
    cr_umove (cr, 0.,top-0.);  
    cr_udraw (cr, 4.,top-0.);
}

int main(int argc, char* argv[])
{
    sf_init(argc,argv);
    
    if (!sf_getfloat("top",&top)) top=5.;
    if (!sf_getfloat("c1",&c1)) c1=0.5;
    if (!sf_getfloat("c2",&c2)) c2=5.;

    cr_main(&argc,&argv);

    exit(0);
}

/* 	$Id: 7ab.c 691 2004-07-04 19:28:08Z fomels $	 */
