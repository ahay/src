/* Pair-wise shifts */
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

#include "dip3.h"

int main(int argc, char* argv[])
{
  bool verb, drift, ***mm;
  int it, nt, two, ie, ne, pair[2], esize, k0, k1, niter, liter, rect[3], order, nj1;
  off_t pos, mpos;
  float **traces, **mtraces, **dip, p0, eps, pmin, pmax, t0, dt;
  sf_file inp, shift, edge, msk;

  esize = sizeof(float);

  sf_init(argc,argv);
  edge = sf_input("in");
  inp = sf_input("data");
  shift = sf_output("out");

  if (SF_INT != sf_gettype(edge)) sf_error("Need int edge");
  if (!sf_histint(edge,"n1",&two) || two != 2) sf_error("No n1=2 in edge");
  if (!sf_histint(edge,"n2",&ne)) sf_error("No n2= in edge");

  if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
  if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
  if (!sf_histfloat(inp,"d1",&dt)) dt=1.0f;
  if (!sf_histfloat(inp,"o1",&t0)) t0=0.0f;

  sf_putint(shift,"n1",nt);
  sf_putfloat(shift,"d1",dt);
  sf_putfloat(shift,"o1",t0);
  sf_putint(shift,"n2",ne);
  sf_settype(shift, SF_FLOAT);

  traces = sf_floatalloc2(nt,2);
  dip = sf_floatalloc2(nt, 2);

  if (NULL != sf_getstring("mask")) {
    mm = sf_boolalloc3(nt,2,2);
    msk = sf_input("mask");
    mtraces = sf_floatalloc2(nt, 2);
  } else {
    mm = (bool***) sf_alloc(2,sizeof(bool**));
    mm[0] = (bool**) sf_alloc(2,sizeof(bool*));
    mm[0][0] = mm[0][1] = NULL;
    msk = NULL;
    mtraces = NULL;
  }

  if (!sf_getint("niter",&niter)) niter=5;
  /* number of iterations */
  if (!sf_getint("liter",&liter)) liter=20;
  /* number of linear iterations */
  
  if (!sf_getint("rect1",&rect[0])) rect[0]=1;
  /* dip smoothness */
  rect[2]=rect[1]=1;
  
  if (!sf_getfloat("p0",&p0)) p0=0.0f;
  /* initial dip */
  
  if (!sf_getint("order",&order)) order=1;
  /* accuracy order */
  if (!sf_getint("nj1",&nj1)) nj1=1;
  /* antialiasing */
  
  if (!sf_getbool("drift",&drift)) drift=false;
  /* if shift filter */
  
  if (!sf_getbool("verb",&verb)) verb = false;
  /* verbosity flag */
  if (!sf_getfloat("eps",&eps)) eps=0.0f;
  /* regularization */

  if (!sf_getfloat("pmin",&pmin)) pmin = -FLT_MAX;
  /* minimum dip */
  if (!sf_getfloat("pmax",&pmax)) pmax = +FLT_MAX;
  /* maximum dip */

  /* initialize dip estimation */
  dip3_init(nt, 2, 1, rect, liter, eps, false);

  k0 = k1 = -1;
  pos = sf_tell(inp);
  if (NULL != msk) mpos = sf_tell(msk);
  
  for (ie=0; ie < ne; ie++) {
    if (verb) sf_warning("pair %d of %d;", ie+1, ne);
    
    sf_intread(pair,2,edge);

    if (pair[0] != k0) {
      if (NULL != msk) {
	sf_seek(msk,mpos+pair[0]*nt*esize,SEEK_SET);
	sf_floatread(mtraces[0],nt,msk);
      }

      sf_seek(inp,pos+pair[0]*nt*esize,SEEK_SET);
      sf_floatread(traces[0],nt,inp);
      
      k0 = pair[0];
    }

    if (pair[1] != k1) {
      if (NULL != msk) {
	sf_seek(msk,mpos+pair[1]*nt*esize,SEEK_SET);
	sf_floatread(mtraces[1],nt,msk);
      }
      
      sf_seek(inp,pos+pair[1]*nt*esize,SEEK_SET);
      sf_floatread(traces[1],nt,inp);
      
      k1 = pair[1];
    }

    if (NULL != msk) {
      sf_mask32 (false, order, nj1, 1, nt, 2, 1, mtraces[0], mm[0]);
    }

    /* initialize dip */
    for(it=0; it < nt; it++) {
      dip[0][it] = p0;
      dip[1][it] = p0;
    }
	
    /* estimate dip */
    dip3(false, 1, niter, order, nj1, drift, traces[0], dip[0], mm[0][0], pmin, pmax);
	
    /* write t-x dip */
    sf_floatwrite(dip[0],nt,shift);
  }
  if (verb) sf_warning(".");
  
  exit(0);
}
