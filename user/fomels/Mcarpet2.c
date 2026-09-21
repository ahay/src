/* Carpet flattening. */
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

#include <rsf.h>

int main (int argc, char* argv[])
{
  bool adj;
  int ntr, nr, ie, ne, two, it, nt, i, j, n1, **pairs;
  float eps, **traces, *shift, *sft, *picked;
  sf_map4 pick;
  sf_file shifts, edges, rgts, sfts;

  sf_init(argc, argv);

  if (!sf_getbool("adj",&adj)) adj=false;
  /* adjoint flag */

  shifts = sf_input("shift");
  if (SF_FLOAT != sf_gettype(shifts)) sf_error("Need float data in shift");
  if (!sf_histint(shifts,"n1",&nt)) sf_error("Need n1= in shift");
  if (!sf_histint(shifts,"n2",&nr)) sf_error("Need n2= in shift");
  
  if (adj) {
    sfts = sf_input("in");
    rgts = sf_output("out");
    
    if (!sf_getint("ntr",&ntr)) sf_error("Need ntr=");
    /* number of traces */
    
    sf_putint(rgts,"n2",ntr);
  } else {
    rgts = sf_input("in"); 
    sfts = sf_output("out");

    if (SF_FLOAT != sf_gettype(rgts)) sf_error("Need float data in input");
    if (!sf_histint(rgts,"n1",&n1) || n1 != nt) sf_error("Need n1=%d in input", nt);
    if (!sf_histint(rgts,"n2",&ntr)) sf_error("Need n2= in input");

    sf_putint(sfts,"n2",nr);
  }

  edges = sf_input("edges");

  if (SF_INT != sf_gettype(edges)) sf_error("Need integer data in edges");
  if (!sf_histint(edges,"n1",&two) || two != 2) sf_error("Need n1=2 in edges");
  if (!sf_histint(edges,"n2",&ne) || ne != nr) sf_error("Need n2=%d in edges", nr);

  if (!sf_getfloat("eps",&eps)) eps=0.01;
  /* stretch regularization */
  pick = sf_stretch4_init(nt, 0.0, 1.0, nt, eps);

  pairs = sf_intalloc2(2,ne);
  sf_intread(pairs[0],2*ne,edges);

  sft = sf_floatalloc(nt);
  shift = sf_floatalloc(nt);
  picked = sf_floatalloc(nt);
  traces = sf_floatalloc2(nt,ntr);

  if (adj) {
    for (i=0; i < ntr; i++) {
      for (it=0; it < nt; it++) {
	traces[i][it] = 0.0f;
      }
    }
  } else {
    sf_floatread(traces[0],ntr*nt,rgts);
  }

  for (ie=0; ie < ne; ie++) {
    sf_floatread(shift,nt,shifts);
    for (it=0; it < nt; it++) {
      shift[it] += it;
    }
    sf_stretch4_define(pick, shift, false);
    
    i = pairs[ie][0];
    j = pairs[ie][1];
    
    if (adj) {
      sf_floatread(sft,nt,sfts);
      sf_stretch4_invert_adj(false, pick, sft, picked);
      for (it=0; it < nt; it++) {
	traces[i][it] += sft[it];
	traces[j][it] -= picked[it];
      }
    } else {
      sf_stretch4_invert(false, pick, picked, traces[j]);
      for (it=0; it < nt; it++) {
	sft[it] = traces[i][it] - picked[it];
      }
      sf_floatwrite(sft,nt,sfts);
    }
  }

  if (adj) {
    sf_floatwrite(traces[0],ntr*nt,rgts);
  }

  exit(0);
}
  
