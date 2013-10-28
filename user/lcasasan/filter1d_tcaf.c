/* 1D matched filter computed using a 2D window on the data
  Copyright (C) 2010 Politecnico di Milano

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
#include <rsfgee.h>

#include "filter1d_tcaf.h"

static int nt;
static int nw;
static float *M;
static float *m;
static float *d;

void filter1d_tcaf_init(int nt1  /* window time length */,
						 int nw1  /* window space length */,
						 int order /* filter order */,
						 float *M1 /* signal to be matched */)
/*<initialize>*/
{
	nt=nt1;
	nw=nw1;
	M=M1;
	m=sf_floatalloc(nt);
	d=sf_floatalloc(nt+order-1);
}

void filter1d_tcaf_lop(bool adj, bool add, int nm, int nd, float *mm, float *dd)
/*< linear operator >*/
{
	int iw,it;
	//sf_warning("nw=%d nt=%d nt*nw=%d",nw,nt,nt*nw);

	if( nd < (nt+nm-1) * nw ) sf_error("%s: size problem: %d != %d",__FILE__,nd,(nt+nm-1) * nw);
	sf_adjnull (adj, add, nm, nd, mm, dd);

	for (iw=0;iw<nw;iw++) {
		for (it=0;it<nt;it++) {
			m[it]=M[it+iw*nt];
		}

		tcaf1_init (nt  /* data length */,
					m  /* data [k] */);
		if (adj) {

			for (it=0;it<nt+nm-1;it++)
				d[it]=dd[it+iw* (nt+nm-1 ) ];

			tcaf1_lop(adj,true,nm,nt+nm-1,mm,d);
		}
		else {

			tcaf1_lop(adj,false,nm,nt+nm-1,mm,d);
			for (it=0;it<nt+nm-1;it++)
				dd[it+iw*(nt+nm-1 )]=d[it];
		}
	}
}

void filter1d_tcaf_close()
/*<de-allocations>*/
{
	free(m);
	free(d);
}
