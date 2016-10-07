/* Picking local maxima on the first axis with evenly spaced windows. */
/*
  Copyright (C) 2016 University of Texas at Austin

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

int nspace;

void normalize(float *dat, int n)
/* normalization with maximum absolute number */
{
    int i;
    float max=0.;

    for (i=0; i<n; i++){
        if(fabs(dat[i])>max) max=fabs(dat[i]);
    }

    for (i=0; i<n; i++){
        dat[i] /= max;
    }
}

static int pick_compare (const void *p1, const void *p2)
{
    float f1 = cimagf(* (sf_complex*) p1);
    float f2 = cimagf(* (sf_complex*) p2);
    return (f1 < f2)? 1: (f1 > f2)? -1: 0;
}

static int sort_compare (const void *p1, const void *p2)
{
    float f1 = * (float*) p1;
    float f2 = * (float*) p2;
    return (f1 < f2)? 1: (f1 > f2)? -1: 0;
}

int selection(int *itau, float *psemb, int np0)
/* eliminate very close events */
{
	int np, ip, jp;

	np=np0;
	for(ip=0; ip<np0-1; ip++){
		if(ip==np-1) break;

		if(itau[ip]-itau[ip+1]<=nspace){
			if(psemb[ip+1]<psemb[ip]) jp=ip+1;
			else jp=ip;
			for( ;jp<np0-1; jp++){
				itau[jp]=itau[jp+1];
				psemb[jp]=psemb[jp+1];
			}
			ip--;
			np--;
		}
	}

	return np;
}

int main(int argc, char* argv[])
{
    bool parab, removal;
    int i1, n1, i2, n2, ip, np, iw, nw, wsize, wmin, wmax, *itau, *npick;
    float o1, d1, t0, t1, t2, t, a, space, *trace=NULL;
    float min, max, x, *tau, *semb, *psemb;
    sf_complex *pick=NULL;
    sf_file in, out, semblance, npicks;

    sf_init(argc, argv);
    in = sf_input("in");
    out = sf_output("out");
	semblance = sf_input("semblance");
    npicks = sf_output("npicks");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;

    if (!sf_getfloat("min",&min)) min=o1;
    /* minimum value of time */

    if (!sf_getfloat("max",&max)) max=o1+(n1-1)*d1; 
    /* maximum value of time */ 

    if (!sf_getint("np",&np)) np=n1;
    /* maximum number of picks */

    if (!sf_getint("nw",&nw)) nw=1;
    /* number of windows */
    if(nw>np) nw=np;

    if (!sf_getfloat("space",&space)) space=100.; 
    /* minimum distance bewteen picked events */ 
	nspace=space/d1+0.5;

    if (!sf_getbool("parab",&parab)) parab=false;
    /* if y, parabolic interpolation */

    if (!sf_getbool("removal",&removal)) removal=true;
    /* if y, remove adjacent events based on semblance */

    wsize=((max-min)/d1+1+nw/2)/nw;
    wmin=min/d1+0.5;
    wmax=max/d1+0.5;

    sf_putint(out,"n1",np);
    sf_putint(npicks,"n1",n2);
    sf_putint(npicks,"n2",1);
    sf_settype(npicks,SF_INT);

    trace = sf_floatalloc(n1);
    pick = sf_complexalloc(n1);
	tau = sf_floatalloc(np);
	itau = sf_intalloc(np);
	semb = sf_floatalloc(n1);
	psemb = sf_floatalloc(np);
	npick = sf_intalloc(n2);

	for (i2=0; i2 < n2; i2++) {
		sf_floatread(trace,n1,in);
		sf_floatread(semb,n1,semblance);

		/* normalization in equally-spaced windows */
		normalize(trace, n1);

		for (iw=0; iw<nw-1; iw++) {
			normalize(trace+wmin+iw*wsize, wsize);
		}

		normalize(trace+wmin+(nw-1)*wsize, wmax-wmin-(nw-1)*wsize); 

		/* find all maxima */
		t0 = trace[0];
		t1 = trace[1];
		ip = 0;
		for (i1=2; i1 < n1; i1++) {
			t2 = trace[i1];

			if (t1 > t0 && t1 > t2) {

				if(parab){
					/* parabolic approximation */
					t = 0.5*(t2-t0)/(2*t1-t0-t2);
					a = t1+0.25*(t2-t0)*t;

					if (t < -1.) {
						t=-1;
						a=t0;
					} else if (t > 1.) {
						t=1.;
						a=t2;
					} 

					x = o1+(i1-1+t)*d1;
				}else{
					x = o1+(i1-1)*d1;
					a = t1;
				}

				if (x >= min && x <= max) {
					pick[ip] = sf_cmplx(x,a);	
					ip++;
				}
			}

			t0 = t1;
			t1 = t2;
		}

		if (0==ip) {
			pick[0] = sf_cmplx(o1-d1,0.);
			ip++;
		}

		/* sort in amplitude */
		qsort(pick,ip,sizeof(sf_complex),pick_compare);
	
		if (ip>np) ip=np;

		for (i1=0; i1<ip; i1++){
			tau[i1]=crealf(pick[i1]);
		}

		/* sort in coordinate */
		qsort(tau,ip,sizeof(float),sort_compare);

		/* select events with large semblance */
		if(removal){
			for (i1=0; i1<ip; i1++){
				itau[i1] = tau[i1]/d1+0.5;
				psemb[i1] = semb[itau[i1]];
			}
			ip = selection(itau, psemb, ip);
		}

		/* record number of events */
		for(i1=ip; i1<np; i1++){
			itau[i1]=0;
		}
		npick[i2]=ip;

		/* convert to float data type */
		for (i1=0; i1<np; i1++){
			tau[i1]=d1*itau[i1];
		}

		sf_floatwrite(tau,np,out);
	}
	sf_intwrite(npick,n2,npicks);

	exit(0);
} 
