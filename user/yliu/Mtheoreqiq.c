/* Output theoretical interval Q value and equivalent Q value. */
/*
  Copyright (C) 2022 Jilin University
  
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
#include <math.h>

int main(int argc, char *argv[])
{
    bool  verb;
    int   n1, nc, *q=NULL, *at=NULL, *iq=NULL, i1, i, ic;
    float d1, o1, tt, zz;
    float *eq = NULL, *T;

    sf_file out ,teq;

    sf_init(argc, argv);
    
    out = sf_output("out");
    sf_settype(out,SF_INT);
    
    /* verbosity */
    if (!sf_getbool("verb",&verb)) verb = false;

    /* basic parameters */
    if (!sf_getint("n1",&n1)) n1 = 1000;
    /* size of time axis */

    if (!sf_getfloat("d1",&d1)) d1 = 0.001;
    /* sampling on time axis */
    
    if (!sf_getfloat("o1",&o1)) o1 = 0.;
    /* origin on time axis */
    
    if (!sf_getint("nc",&nc)) nc = 1;
    /* number of layer */
    if (nc < 1) sf_error("Need nq >= 1");

    if (nc >= 1) {
	at = sf_intalloc(nc);
	q = sf_intalloc(nc);

	/* snprintf(key,6,"at"); */
	if (!sf_getints("at",at,nc)) {
	    /* (at=[at1,at2,...] layer quality factor) */
	    sf_error("need at");
	}
	
	/* snprintf(key,6,"q"); */
	if (!sf_getints("q",q,nc)) {
	    /* (q=[q1,q2,...] layer quality factor) */
	    for (i = 0; i < nc; i++) {
		q[i] = 100;
	    }
	} 
    }

    sf_putint(out,"n1",n1);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"o1",o1);

    if (NULL != sf_getstring("teq")) {
	teq = sf_output("teq");
	sf_settype(teq,SF_FLOAT);
	eq = sf_floatalloc(n1);
	
	sf_putint(teq,"n1",n1);
	sf_putfloat(teq,"d1",d1);
	sf_putfloat(teq,"o1",o1);
    } else {
	teq = NULL;
    }
    
    iq = sf_intalloc(n1);
    T = sf_floatalloc(nc);

    for (i = 0; i < nc; i++) {
	T[i] = at[i] * d1;	
	for (i1 = 0; i1 < n1; i1++) {
	    tt = o1 + i1 * d1;
	    if (i == 0) {
		if (tt <= T[0]) {
		    iq[i1] = q[0];
		}
	    }
	    if (i > 0 && i < nc) {
		if (T[i-1] < tt && tt <= T[i]) {
		    iq[i1] = q[i];
		}
	    }
	    if (i == nc-1) {
		if (tt >T[i]) {
		    iq[i1] = q[i];
		}
	    }
	}
    }
    
    sf_intwrite(iq,n1,out);
    
    /* theoretic equivalent Q value */
    if (NULL != teq) {
	for (i = 0; i < nc; i++) {
	    for (i1 = 0; i1 < n1; i1++) {
		if (i == 0) {
		   if (i1 <= at[0]) {
		       eq[i1] = q[0];
		   }
		}

		if (i == 1) {
		   if (i1 > at[0] && i1 <=at[1]) {
		       eq[i1] = q[1];
		   }
		}

		zz = 0.;	
		if (i > 1) {
		    for (ic = 1; ic < i; ic++) {
			zz += (at[ic]-at[ic-1])*d1/q[ic];
		    }
		    if (i1 > at[i-1] && i1 <= at[i]) {
			eq[i1] = (i1-at[0])*d1/(zz+(i1-at[i-1])*d1/q[i]);
		    }
		}
		if (i == nc-1) {
		    if (i1 > at[nc-1]) {
			eq[i1] = (i1-at[0])*d1/(zz+(i1-at[nc-2])*d1/q[nc-1]);
		    }		    
		}
	    }
	}
    }
    if (NULL != teq) sf_floatwrite(eq,n1,teq);
    
}
