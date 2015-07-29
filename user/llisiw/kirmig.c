/* Utilities for Kirchhoff migration. */
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

#include <rsf.h>

static int nt;
static float dt; 
static float t0;

void doubint(int nt, float *trace)
/*< causal and anticausal integration >*/
{
    int it;
    float tt;

    tt = trace[0];
    for (it=1; it < nt; it++) {
	tt += trace[it];
	trace[it] = tt;
    }
    tt = trace[nt-1];
    for (it=nt-2; it >=0; it--) {
	tt += trace[it];
	trace[it] = tt;
    }
}

void kirmig_init(int n1, float d1, float o1)
/*< initialize trace parameters >*/
{
    nt=n1;
    dt=d1;
    t0=o1;
}

void pick(bool adj, float ti, float deltat, float *sample, float *trace)
/*< pick a traveltime sample from a trace >*/
{
    int it, itm, itp;
    float ft, tm, tp, ftm, ftp, imp, s;

    ft = (ti-t0)/dt; it = floorf(ft); ft -= it; 
    if ( it < 0 || it >= nt-1) return;
 
    tm = ti-deltat-dt;
    ftm = (tm-t0)/dt; itm = floorf(ftm); ftm -= itm; 
    if (itm < 0) return;
                 
    tp = ti+deltat+dt;
    ftp = (tp-t0)/dt; itp = floorf(ftp); ftp -= itp; 
    if (itp >= nt-1) return;

    imp = dt/(dt+tp-tm);
    imp *= imp;
		
    if (adj) {
	sample[0] += imp*(
	    2.*(1.-ft)*trace[it] + 2.*ft*trace[it+1] -
	    (1.-ftm)*trace[itm] - ftm*trace[itm+1]   - 
	    (1.-ftp)*trace[itp] - ftp*trace[itp+1]);    
    } else {
	s = imp*sample[0];
	
	trace[it]    += 2.*(1.-ft)*s;
	trace[it+1]  += 2.*ft*s;
	trace[itm]   -= (1.-ftm)*s;
	trace[itm+1] -= ftm*s;
	trace[itp]   -= (1.-ftp)*s;
	trace[itp+1] -= ftp*s;
    }	
}
