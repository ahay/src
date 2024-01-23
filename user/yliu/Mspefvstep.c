/* Streaming prediction-error filter with variable step. */
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

void gspef_loop(float *inp, float *res);

float gspef_estimate_constraint1(float data, float *tmp, int i1);

float gspef_estimate_constraint2(float data, float *tmp, int i1);

void gspef_pick(float *inp, float *tmp, int i2, int i1, int nw);

void init_flt();

void sweep_flt(int flt_case);

bool verb;
int n1, n2, n3, na;
float lambda, lambda1, lambda2;
float *ftmp1, *ftmp2;

int main(int argc, char *argv[]) {

    float *dinp, *dout;

    sf_file inp, out;

    sf_init(argc, argv);

    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("verb", &verb))
	verb = false;

    if (!sf_histint(inp, "n1", &n1))
	sf_error("No n1= in input");
    if (!sf_histint(inp, "n2", &n2))
	n2 = 1;
    if (!sf_histint(inp, "n3", &n3))
	n3 = 1;
    n2 = n2 * n3;
    if (verb)
	sf_warning("n1=%d n23=%d", n1, n2);
    /* figure out data size */

    if (!sf_getint("na", &na)) sf_error("Need na=");
    /* filter size */

    if (!sf_getfloat("lambda1", &lambda1)) sf_error("Need lambda1=");
    /* constant scale parameter in time axis */
 
   if (!sf_getfloat("lambda2", &lambda2)) lambda2 = 0.0f;
    /* constant scale parameter in space axis */

    if (verb) sf_warning("lambda1=%f lambda2=%f", lambda1, lambda2);
    
    lambda1 *= lambda1;
    lambda2 *= lambda2;
    lambda = lambda1 + lambda2;
    /* set up lambda parameters */

    dinp = sf_floatalloc(n1 * n2);
    dout = sf_floatalloc(n1 * n2);

    sf_floatread(dinp, n1 * n2, inp);
    /* read input data  */

    gspef_loop(dinp, dout);
    /* deconvolution with SPEF with variable step*/

    sf_floatwrite(dout, n1 * n2, out);
    /* write output data */
}

/* arrange looping order  */
void gspef_loop(float *inp, float *res) 
{
    int ia1, i1, i2, pt, nf;
    int *nw, *ngp;
    float *tmp, *infilt;
    sf_file lag, infil;

    if (NULL != sf_getstring("infil")) { /* initial filter */
	infil = sf_input("infil");
	if(!sf_histint(infil, "n1", &nf)) nf = na;
	if (nf!=na) sf_error("nf need= na");
	infilt = sf_floatalloc(nf);
    } else {
	infil = NULL;
	infilt = NULL;
    }
    
    if (NULL != sf_getstring("lag")) { /* variable step */
	lag = sf_input("lag");
	ngp = sf_intalloc(n1);
    } else {
	lag = NULL ;
	ngp = NULL ;
    }

    nw = sf_intalloc(n1);
    tmp = sf_floatalloc(na);
    /* restore picked data */

    init_flt();
    /* initialize filters */

    for (i2 = 0; i2 < n2; i2++) {

	sweep_flt(1);
	/* set ftmp1 to be zero */
	if (NULL != infil) {
	    sf_floatread(infilt, nf, infil);
	}
	if (NULL != infil) {
	    for (ia1=0; ia1 < na; ia1++) {
		ftmp1[ia1] = infilt[ia1];
	    } /* initialize ftmp1 with global filter */
	} else {
	    for (ia1=0; ia1 < na; ia1++) {
		ftmp1[ia1] = 0.0f;
	    } /* initialize ftmp1 with zero */
	}

	if (NULL != lag) {
	    sf_intread(ngp,n1,lag);
	}
	for (i1 = 0; i1 < n1; i1++) {
	    nw[i1] = 0;
	}
	for (i1 = 0; i1 < n1; i1++) {
	    nw[i1] = na + ngp[i1];
	}

	for (i1 = 0; i1 < n1; i1++) {

	    if (verb) sf_warning("loop over %3d %3d;", i2+1, i1+1);
	    pt = i2 * n1 + i1;

	    gspef_pick(inp, tmp, i2, i1, nw[i1]);
	    /* SPEF for calculating data points */

	    if (i2 == 0) {
		res[pt] = gspef_estimate_constraint1(inp[pt], tmp, i1);
	    } else {
		res[pt] = gspef_estimate_constraint2(inp[pt], tmp, i1);
	    }
	}
    }

    if (verb) sf_warning(".");

    free(ftmp1);
    free(ftmp2);
}

float gspef_estimate_constraint1(float data, float *tmp, int i1) 
{
    int i;
    float dTd, dTa, rn, res, A, B, C, D;

    dTd = 0.0f;
    for (i = 0; i < na; i++) {
	dTd += tmp[i] * tmp[i];
    } /* dT*d */

    dTa = 0.0f;
    for (i = 0; i < na; i++) {
	A = lambda1 * ftmp1[i];
	B = tmp[i] * A / lambda1;
	dTa += B;
    } /* dT*a */

    rn = (data + dTa) / (lambda1 + dTd);
    res = lambda1 * rn;

    for (i = 0; i < na; i++) {
	C = lambda1 * ftmp1[i];
	D = C / lambda1;

	ftmp1[i] = D - (rn * tmp[i]);
	ftmp2[i1 * na + i] = ftmp1[i];
    } /* update ftmp1, ftmp2 */

    return res;
}

float gspef_estimate_constraint2(float data, float *tmp, int i1) 
{
    int i;
    float dTd, dTa, rn, res, A, B, C, D;

    dTd = 0.0f;
    for (i = 0; i < na; i++) {
	dTd += tmp[i] * tmp[i];
    } /* dT*d */

    dTa = 0.0f;
    for (i = 0; i < na; i++) {
	A = lambda1 * ftmp1[i] + lambda2 * ftmp2[i1 * na + i];
	B = tmp[i] * A / lambda;
	dTa += B;
    } /* dT*a */

    rn = (data + dTa) / (lambda + dTd);
    res = lambda * rn;

    for (i = 0; i < na; i++) {
	C = lambda1 * ftmp1[i] + lambda2 * ftmp2[i1 * na + i];
	D = C / lambda;

	ftmp1[i] = D - (rn * tmp[i]);
	ftmp2[i1 * na + i] = ftmp1[i];
    } /* update ftmp1, ftmp2 */

    return res;
}

void gspef_pick(float *inp, float *tmp, int i2, int i1, int nw) {
    
    int pi, pt;
    for (pi = 0; pi < na; pi++) {

	pt = i1 - nw + pi;

	if (pt < 0) { /* prevent dT cross the boundary */
	    tmp[pi] = 0.0f;
	} else { 
	    tmp[pi] = inp[i2 * n1 + pt];
	}
    }
}

/* initialize filter */
void init_flt()
{
    ftmp1  = sf_floatalloc(na);
    ftmp2  = sf_floatalloc(n1*na);
    /* memory allocation */

    sweep_flt(1);
    sweep_flt(2);
    /* sweep filter with zero */

}

/* sweeping filter */
void sweep_flt(int flt_case)
{
    int ia;
    
    if (flt_case == 1) {	
	for (ia=0; ia < na; ia++) {
	    ftmp1[ia] = 0.0f;
	} /* initialize ftmp1 with zero */
    } else if (flt_case ==2) {
	for (ia=0; ia < na*n1; ia++) {
	    ftmp2[ia] = 0.0f;
	}/* initialize ftmp2 */

    } else {
	sf_error("Pick wrong case!");
    }

}

