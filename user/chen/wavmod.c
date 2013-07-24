/* 3D wave modeling */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>
#include "fd3.h"

#include "obsys.h"
/*^*/

static HVel vel;
static int nxyz, nt, st, ntj, jt, jtm;
static float *wav, *owv, *nwv, *ud, **data;
static int *pg, ng;
static bool verb;

void wavmod_init(HVel pvel,
	float dt, int nt0, int st0, int jt0, int jtm0,
	int *pg0, int ng0, bool vb)
/*< initialize >*/
{
	int i1, n1, n2, n3;
	double dt2;
	float d1, d2, d3;

	vel = pvel; verb = vb;
	nt = nt0;	st = st0;	jt = jt0;	jtm= jtm0;
	pg = pg0;	ng = ng0;

	n1 = sf_n(vel->z); n2=1; n3=1;
	d1 = sf_d(vel->z); d2=0.0; d3=0.0;
	if(vel->nd >= 2) 
	{
		n2 = sf_n(vel->x);
		d2 = sf_d(vel->x);
	}
	if(vel->nd >= 3) 
	{
		n3 = sf_n(vel->y);
		d3 = sf_d(vel->y);
	}

	nxyz = n1*n2*n3;

	owv = sf_floatalloc(nxyz);
	wav = sf_floatalloc(nxyz);
	nwv = sf_floatalloc(nxyz);
	ud  = sf_floatalloc(nxyz);
	ntj = (nt-st+1)/jt;
	data  = sf_floatalloc2(ntj, ng);

	fd3_init(d1, d2, d3);
	dt2 = dt*dt;
	for(i1=0; i1<nxyz; i1++)
		vel->p[i1] = dt2*vel->p[i1]*vel->p[i1];
}

void wavmod_shot(sf_file dat, sf_file wfl, int ns, int *ps, float **ws)
/*< shot modelling >*/
{
	int it, ig, is;
	int i1, n1, n2, n3;
	float *p;

	n1 = sf_n(vel->z); n2=1; n3=1;
	if(vel->nd >= 2) n2 = sf_n(vel->x);
	if(vel->nd >= 3) n3 = sf_n(vel->y);

	memset(owv, 0, nxyz*sizeof(float));
	memset(wav, 0, nxyz*sizeof(float));
	memset(nwv, 0, nxyz*sizeof(float));
	memset(ud, 0, nxyz*sizeof(float));

	for (it=0; it < nt; it++)
	{
		fd3_laplacian(n1, n2, n3,  wav, ud);

		for (is=0; is<ns; is++)
		ud[ps[is]] += ws[is][it];

#ifdef _OPENMP
#pragma omp parallel for	 \
	schedule(dynamic,n1)	  \
	private(i1)		  
#endif
		for(i1=0; i1<nxyz; i1++)
			nwv[i1] = 2.0*wav[i1] - owv[i1] + ud[i1]*vel->p[i1];
		p = owv;	
		owv = wav;
		wav = nwv;
		nwv = p;

		if(it >= st && (it-st)%jt==0)
		{
#ifdef _OPENMP
#pragma omp parallel for	 \
	schedule(dynamic,8)	  \
	private(ig)		  
#endif
			for(ig=0; ig<ng; ig++)
			data[ig][(it-st)/jt] = wav[pg[ig]];
		}
		if(wfl!=NULL && it>=st && (it-st)%jtm == 0) /* wave */
		sf_floatwrite(wav, nxyz, wfl);
		if(verb) sf_warning("%d of %d;", it, nt);
	}
	/* output seismic data */
	sf_floatwrite(data[0], ntj*ng, dat);
}

void wavmod_close()
/*< release memory >*/
{
	free(data[0]);
	free(data);
	free(owv);
	free(wav);
	free(nwv);
	free(vel);
	free(ud);
}


