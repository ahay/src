/* Illustrating non-hyperbolic approximations */
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
#include <float.h>

#include <rsf.h>

int main(int argc, char* argv[]) 
{
    char *dist, *appr;
    int np, nq, ip, iq;
    float dp, dq, p, q, q2, q1, t=0., ta=0., *trace, eps;
    sf_file out;

    sf_init(argc,argv);
    out = sf_output("out");

    sf_setformat(out,"native_float");

    if (!sf_getint("np",&np)) np=300;
    if (!sf_getint("nq",&nq)) nq=300;
    
    if (!sf_getfloat("dp",&dp)) dp = 1./(np-1);
    if (!sf_getfloat("dq",&dq)) dq = 4./(nq-1);

    eps=0.001;

    sf_putint(out,"n1",np);
    sf_putint(out,"n2",nq);
    sf_putfloat(out,"d1",dp);
    sf_putfloat(out,"d2",dq);
    sf_putfloat(out,"o1",0.);
    sf_putfloat(out,"o2",1.+eps);

    dist = sf_getstring("dist");
    /* distribution type */
    if (NULL == dist) sf_error("Need dist=");

    appr = sf_getstring("appr");
    /* approximation type */
    if (NULL == appr) sf_error("Need appr=");

    trace = sf_floatalloc(np);

    for (iq=0; iq < nq; iq++) {
	q = 1. + eps + iq*dq;
	q2 = q*q;
	for (ip=0; ip < np; ip++) {
	    p = ip*dp;
	    p *= p;
	    switch (dist[0]) {
		case 'l': /* linear */
		    q1 = logf(q);
		    t = logf((1 + q2 + p*(q2-1) + 
			      sqrtf(1. + q2*(q2-2.) + 
				    p*(2.*q2*q2-2. + p*(1 + q2*(q2-2.)))))/
			     (2.*q))/q1;
		    switch (appr[0]) {
			case 'h': /* hyperbolic */
			    ta = sqrtf(1. + 2.*p/q1);
			    break;
			case 'm': /* Malovichko */
			    ta = 1. + (q2-1.)*(
				sqrtf((q2-1. + 2*p*(1 + q2))/
				      (q2-1.))-1.)/((1 + q2)*q1);
			    break;
			case 'f': case 's': case 'n': /* Fomel-Stovas */
			    ta = 1 + 2*p/q1 +
				(2*p*p*(q2-1. - q1*(q2+1.)))/
				(q1*q1*(q2-1.)*
				 (1 + (p*(1 - q2 + q1*(q2+1.)))/(q1*(q2-1.)) + 
				  sqrtf(powf(-(p*(q2-1.)) + 
					     q1*(q2-1. + p*(q2+1.)),2)/
					powf(q2-1.,2) - 
					p*p*(2*sqrtf(3)*
					     powf((q1+1. - q2 + q1*q2)/
						  (q2-1.),1.5)
					     - (3*(1.-q2 + q1*(q2+1.)))/
					     (q2-1.)))/q1));
			    ta = sqrtf(ta);
			    break;
			default:
			    sf_error("Unknown approximation %s",appr);
			    break;
		    }	    
		    break;
		default:
		    sf_error("Unknown distribution %s",dist);
		    break;
	    }
	    trace[ip] = 100.*fabsf(t-ta)/t;
	}
	sf_floatwrite(trace,np,out);
    }

    exit(0);
}
