/* Extract horizons from data */
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

#include "fint1.h"

int main(int argc, char* argv[])
{
  int nt, nh, n2, it, ih, i2;
  float *trace, *coord, t0, dt, t;
  fint1 str;
  sf_file in, out, pick;

  sf_init(argc,argv);
  in = sf_input("in");
  pick = sf_input("pick");
  out = sf_output("out");

  if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
  if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
  if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
  n2 = sf_leftsize(in,1);

  if (!sf_histint(pick,"n1",&nh)) sf_error("No n1= in pick");
  if (n2 != sf_leftsize(pick,1)) sf_error("Wrong number of traces in pick (need %d)",n2);

  sf_putint(out,"n1",nh);

  trace = sf_floatalloc (nt);
  coord = sf_floatalloc (nh);

  str = fint1_init (4,nt,0); 

  for (i2=0; i2 < n2; i2++) {
    sf_floatread (trace,nt,in);
    sf_floatread (coord,nh,pick);
    
    fint1_set(str,trace);
    
    for (ih=0; ih < nh; ih++) {
      t = (coord[ih]-t0)/dt;
      it = floorf(t);
      if (it >=0 && it < nt) {
	coord[ih] = fint1_apply(str,it,t-it,false);
      }
    }

    sf_floatwrite(coord,nh,out);
  }

  exit(0);
}
