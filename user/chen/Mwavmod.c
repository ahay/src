/* 1-2-3D finite difference modeling */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
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
#include "wavmod.h"

int main(int argc, char* argv[])
{
    int jt, jtm, is, st;			/* index */
    int nt, n1, n2, ns, ng, m[SF_MAX_DIM], nm;				 	/* dimensions */
    float ot, o1, d1;				/* original */
    sf_file in, dat, wfl, vel, sgrid, ggrid; /* I/O files */
    float **wvlt;
    sf_axis ax;
    HVel hv;
    HCoord hs, hg;
    bool verb;

    sf_init(argc,argv);

    in = sf_input ("in");   
/* source wavelet \n
   \t nt X 1: regular shot gather
   \t nt X ns: simultaneous source shotting
*/
    vel  = sf_input ("vel");  /* velocity field */
    sgrid  = sf_input ("sgrid");  /* source grid */
    ggrid  = sf_input ("ggrid");  /* geophone grid */
    dat = sf_output("out");  /* seismic data */
 
    if(sf_getstring("wfl")!=NULL) 
	wfl = sf_output("wfl"); /* wavefield movie file */
    else wfl=NULL;

    if (!sf_getint("jt",&jt)) jt=1; 
    /* time interval in observation system */
    if (!sf_getint("jtm",&jtm)) jtm=100; 
    /* time interval of wave movie */
    if (!sf_getfloat("ot", &ot)) ot = 0.0; 
    /* time delay */
    if (!sf_getbool("verb", &verb)) verb = false; 
    /* verbosity */

    /* velocity and observation system */
    hv = obs_vel(vel);
    m[0] = sf_n(hv->z);	nm = 1;
    if(hv->nd >= 2)	{m[1] = sf_n(hv->x); nm=2;}
    if(hv->nd >= 3)	{m[2] = sf_n(hv->y); nm=3;}
    hs = obs_coord(sgrid, m, nm);
    hg = obs_coord(ggrid, m, nm);
    ns = sf_n(hs->a2);
    ng = sf_n(hg->a2);

    /* waveform */
    ax = sf_iaxa(in, 1);
    n1 = sf_n(ax);	o1 = sf_o(ax);	d1 = sf_d(ax);
    if(!sf_histint(in, "n2", &n2) || n2 != ns) n2=1;
    wvlt = sf_floatalloc2(n1, n2);
    sf_floatread(wvlt[0], n1*n2, in);


    if(ot<o1) ot=o1;
    st = (ot-o1)/d1;
    nt = (n1-st+1)/jt;
    sf_setn(ax, nt);
    sf_setd(ax, d1*jt);
    sf_seto(ax, ot);
    sf_oaxa(dat, ax, 1);
    sf_oaxa(dat, hg->a2, 2);
    if(n2==1) sf_oaxa(dat, hs->a2, 3);

    if(wfl!=NULL)
    {
	sf_oaxa(wfl, hv->z, 1);
	if(hv->nd >= 2) sf_oaxa(wfl, hv->x, 2);
	if(hv->nd >= 3) sf_oaxa(wfl, hv->y, 3);
	sf_setn(ax, (n1-st+1)/jtm);
	sf_setd(ax, d1*jtm);
	sf_seto(ax, ot);
	sf_oaxa(wfl, ax, hv->nd+1);
	if(n2==1) sf_oaxa(wfl, hs->a2, hv->nd+2);
    }

    wavmod_init(hv, d1, n1, st, jt, jtm, hg->p, ng, verb);

    if(n2==1)
	for (is=0; is < ns; is++)
	{
	    wavmod_shot(dat, wfl, 1, hs->p+is, wvlt);
	    if(verb) sf_warning("shot %d of %d", is, ns/n2);
	}
    else 
	wavmod_shot(dat, wfl, ns, hs->p, wvlt);

    wavmod_close();

    free(wvlt[0]);
    free(wvlt);
    return (0);
}

