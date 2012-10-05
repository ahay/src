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

void doubint(int nt, float *trace)
/*< causal and anticausal integratio >*/
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

    
float pick(float ti, float deltat, 
	   const float *trace,
	   int nt, float dt, float t0)
/*< pick a traveltime sample from a trace >*/
{
    int it, itm, itp;
    float ft, tm, tp, ftm, ftp, imp;

    ft = (ti-t0)/dt; it = floorf(ft); ft -= it; 
    if ( it < 0 || it >= nt-1) return 0.0;
 
    tm = ti-deltat-dt;
    ftm = (tm-t0)/dt; itm = floorf(ftm); ftm -= itm; 
    if (itm < 0) return 0.0;
                 
    tp = ti+deltat+dt;
    ftp = (tp-t0)/dt; itp = floorf(ftp); ftp -= itp; 
    if (itp >= nt-1) return 0.0;

    imp = dt/(dt+tp-tm);
    imp *= imp;
		
    return imp*(
	2.*(1.-ft)*trace[it] + 2.*ft*trace[it+1] -
	(1.-ftm)*trace[itm] - ftm*trace[itm+1]    - 
	(1.-ftp)*trace[itp] - ftp*trace[itp+1]);    
}
