/* Double square-root eikonal solver */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

#include "dsreiko.h"

int main(int argc, char* argv[])
{
    float a, b, c, d, e;
    sf_complex root[4];

    sf_init(argc,argv);
    
    if (!sf_getfloat("a",&a)) a=1.;
    if (!sf_getfloat("b",&b)) b=1.;
    if (!sf_getfloat("c",&c)) c=1.;
    if (!sf_getfloat("d",&d)) d=1.;
    if (!sf_getfloat("e",&e)) e=1.;

    ferrari((double)a,(double)b,(double)c,(double)d,(double)e,root);
    
    sf_warning("root 1: (%g,%g)",crealf(root[0]),cimagf(root[0]));
    sf_warning("root 2: (%g,%g)",crealf(root[1]),cimagf(root[1]));
    sf_warning("root 3: (%g,%g)",crealf(root[2]),cimagf(root[2]));
    sf_warning("root 4: (%g,%g)",crealf(root[3]),cimagf(root[3]));

    exit(0);
}
