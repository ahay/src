/* 2D missing data interpolation using f-x streaming prediction filter. */
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

/* control calculation loop */
void loop_control(kiss_fft_cpx *inp, kiss_fft_cpx *ext_data,
                  kiss_fft_cpx *conj_ext_data);

/* filter memory allocate and initialize */
void filter_init();

/* filter sweeping */
void filter_sweep(int FILTER_CASE);

/* filter estimate and data interpolation */
kiss_fft_cpx filter_estimete(kiss_fft_cpx data, kiss_fft_cpx *stream,
                             kiss_fft_cpx *conj_stream, int CONSTRAINT_TYPE,
                             int MASK_TYPE, int i2);

/* pick data stream for calculation */
void data_stream_pick(kiss_fft_cpx *inp, kiss_fft_cpx *stream, int i2, int i1);

/* add boundary for data calculation */
void bound_add(kiss_fft_cpx *data, kiss_fft_cpx *ext_data,
               kiss_fft_cpx *conj_ext_data, int bn1, int bn2, int bns);

/* forward FFT */
void forw_fft(float *dd, kiss_fft_cpx *ff, int fnt, int fnw);

/* backward FFT */
void back_fft(kiss_fft_cpx *ff, float *dd, int bnt, int bnw);

bool VERB;
int FUNC_TYPE, FILTER_TYPE, BOUND_TYPE, BOUND_SIZE;

int n1, nw, n2, n3, n23, na, naa;
float o1, d1;
int *mask;
kiss_fft_cpx lambda, lambda1, lambda2;
kiss_fft_cpx *ftmp1, *ftmp2;
kiss_fft_cpx zero;

int main(int argc, char *argv[])
{
    int nt;
    float *dinp, *dout;
    kiss_fft_cpx *dfft, *temp, *ctmp;
    sf_file inp, out, msk = NULL;
    
    sf_init(argc, argv);
    
    inp = sf_input("in");
    out = sf_output("out");
    
    if (!sf_getbool("verb", &VERB)) VERB = false;
    /* default=false, verbosity flag */
    
    FUNC_TYPE = 2;
    /* default = 2,
       1 : denoise
       2 : interpolation*/
    
    if (!sf_getint("ftype", &FILTER_TYPE))
	FILTER_TYPE = 1;
    /* default = 1,
       1 : causal filter design
       2 : noncausal filter design */
    
    
    zero.r = 0.0f;
    zero.i = 0.0f;
    
    if (!sf_histfloat(inp, "o1", &o1))
	o1 = 0.;
    if (!sf_histfloat(inp, "d1", &d1))
	d1 = 1.;
    if (!sf_histint(inp, "n1", &n1))
	sf_error("No n1= in input");
    if (!sf_histint(inp, "n2", &n2))
	n2 = 1;
    if (!sf_histint(inp, "n3", &n3))
	n3 = 1;
    n23 = n2 * n3;
    
    /* initialize FFT setting */
    nt = 2 * kiss_fft_next_fast_size((n1 + 1) / 2);
    if (nt % 2)
	nt++;
    nw = nt / 2 + 1;
    /* dw = 1. / (nt * d1); */
    
    lambda = lambda1 = lambda2 = zero;
    if (!sf_getfloat("lambda1", &lambda1.r))
	sf_error("Need lambda1=");
    /* lambda1.r */
    if (!sf_getfloat("lambda2", &lambda2.r))
	lambda2 = zero;
    /* lambda2.r */
    lambda1 = sf_cmul(lambda1, lambda1); /* lambda1^2 */
    lambda2 = sf_cmul(lambda2, lambda2); /* lambda2^2 */
    lambda = sf_cadd(lambda1, lambda2);  /* lambda */
    
    if (!sf_getint("na", &na))
	sf_error("Need na=");
    if (FILTER_TYPE == 2) {
	na = na * 2;
    }
    if (VERB)
	sf_warning("Filter size %d", na);
    BOUND_SIZE = na;
    
    /* memory allocate */
    dinp = sf_floatalloc(n1 * n23);
    dout = sf_floatalloc(n1 * n23);
    dfft = (kiss_fft_cpx *)sf_complexalloc(nw * n23);
    
    temp = (kiss_fft_cpx *)sf_complexalloc(nw * (n2 + BOUND_SIZE));
    ctmp = (kiss_fft_cpx *)sf_complexalloc(nw * (n2 + BOUND_SIZE));
    if (VERB)
	sf_warning("Allocating memory done.");
    
    if (FUNC_TYPE == 2) {
	msk = sf_input("mask");
	if (SF_INT != sf_gettype(msk))
	    sf_error("Need int type in mask");
	/* mask */
	
	mask = sf_intalloc(n23);
	/* memory allocate */
	
	sf_intread(mask, n23, msk);
	if (VERB)
	    sf_warning("Read mask data done");
	/* read mask data */
    }
    
    /* read data */
    sf_floatread(dinp, n1 * n23, inp);
    if (VERB)
	sf_warning("Reading data done.");
    
    /* fft with nft points */
    forw_fft(dinp, dfft, nt, nw);
    if (VERB)
	sf_warning("FFT transform with nft done.");
    
    /* add boundary for dfft */
    bound_add(dfft, temp, ctmp, nw, n2, BOUND_SIZE);
    if (VERB)
	sf_warning("Adding boundary for dfft done.");
    
    /* denoise/ interpolation */
    loop_control(dfft, temp, ctmp);
    /* estima_flt(dfft, temp, ctmp); */
    
    /* fft-1 */
    back_fft(dfft, dout, nt, nw);
    if (VERB)
	sf_warning("Inverse FFT transform done.");
    
    /* output denoise/recover data */
    sf_floatwrite(dout, n1 * n23, out);
    
    return 0;
}

/* control calculation loop */
void loop_control(kiss_fft_cpx *inp, kiss_fft_cpx *ext_data,
                  kiss_fft_cpx *conj_ext_data)
{
    int i1, i2, pt;
    kiss_fft_cpx *stream, *conj_stream;
    
    stream = (kiss_fft_cpx *)sf_complexalloc(na);
    conj_stream = (kiss_fft_cpx *)sf_complexalloc(na);
    
    filter_init();
    
    for (i1 = 0; i1 < nw; i1++) { /* loop over frequency axis */
	
	filter_sweep(1);
	
	for (i2 = 0; i2 < n2; i2++) { /* loop over x axis */
	    
	    if (VERB)
		sf_warning("loop over [%d/%d] [%d/%d];",
			   i1 + 1, nw, i2 + 1, n2);
	    pt = i2 * nw + i1;
	    
	    data_stream_pick(ext_data, stream, i2, i1);
	    data_stream_pick(conj_ext_data, conj_stream, i2, i1);
	    
	    if (lambda2.r <= 0.) {
		inp[pt] =
		    filter_estimete(inp[pt], stream, conj_stream, 
				    1, mask[i2], i2);
	    }
	    else {
		if (i1 == 0) {
		    inp[pt] = filter_estimete(inp[pt], stream, 
					      conj_stream, 12, mask[i2], i2);
		}
		else {
		    inp[pt] = filter_estimete(inp[pt], stream, 
					      conj_stream, 12, mask[i2], i2);
		}
	    }
	    
	    if (mask[i2] == 0) {
		ext_data[(i2 + BOUND_SIZE) * nw + i1] = inp[pt];
		conj_ext_data[(i2 + BOUND_SIZE) * nw + i1] = sf_conjf(inp[pt]);
	    } /* update temp of d */
	}
    }
    
    if (VERB)
	sf_warning(".");
    
    free(ftmp1);
    free(ftmp2);
}

/* filter estimate and data interpolation */
kiss_fft_cpx filter_estimete(kiss_fft_cpx data, kiss_fft_cpx *stream,
                             kiss_fft_cpx *conj_stream, int CONSTRAINT_TYPE,
                             int MASK_TYPE, int i2)
{
    int ia;
    kiss_fft_cpx dTd, dTa, dn, rn, res, A, B, C, D;
    
    if (CONSTRAINT_TYPE == 12) {
	lambda = sf_cadd(lambda1, lambda2);
    } else if (CONSTRAINT_TYPE == 1) {
	lambda = lambda1;
    }
    
    dTd = zero;
    for (ia = 0; ia < na; ia++) {
	dTd = sf_cadd(dTd, sf_cmul(stream[ia], conj_stream[ia]));
    } /* update dd */
    
    dTa = zero;
    for (ia = 0; ia < na; ia++) {
	if (CONSTRAINT_TYPE == 12) {
	    A = sf_cadd(sf_cmul(lambda1, ftmp1[ia]),
			sf_cmul(lambda2, ftmp2[i2 * na + ia]));
	} else if (CONSTRAINT_TYPE == 1) {
	    A = sf_cmul(lambda, ftmp1[ia]);
	} else {
	    sf_error("Error of constraint type!");
	}
	
	B = sf_cmul(stream[ia], sf_cdiv(A, lambda));
	dTa = sf_cadd(dTa, B);
    } /* update da */
    
    if (FUNC_TYPE == 1) { /* denoise */
	
	rn = sf_cdiv(sf_cadd(data, dTa), sf_cadd(lambda, dTd));
	res = sf_cmul(lambda, rn);
	dn = sf_csub(data, res);
    } else if (FUNC_TYPE == 2) { /* interpolation */
	
	if (MASK_TYPE) { /* calculate filter */
	    dn = data;
	    rn = sf_cdiv(sf_cadd(data, dTa), sf_cadd(lambda, dTd));
	    res = sf_cmul(lambda, rn);
	} else { /* interpolate gap */
	    dn = sf_csub(zero, dTa);
	    rn = zero;
	}
    } else {
	sf_error("Wrong type been chosed!");
    }
    
    for (ia = 0; ia < na; ia++) {
	
	if (CONSTRAINT_TYPE == 12) {
	    
	    C = sf_cadd(sf_cmul(lambda1, ftmp1[ia]),
			sf_cmul(lambda2, ftmp2[i2 * na + ia]));
	} else if (CONSTRAINT_TYPE == 1) {
	    C = sf_cmul(lambda1, ftmp1[ia]);
	}
	
	D = sf_cdiv(C, lambda);
	ftmp1[ia] = sf_csub(D, sf_cmul(rn, conj_stream[ia]));
	
	if (lambda2.r > 0.)
	    ftmp2[i2 * na + ia] = ftmp1[ia];
    } /* update ftmp1, ftmp2, filter */
    
    return dn;
}

/* filter memory allocate and initialize */
void filter_init()
{
    /* memory allocation */
    ftmp1 = (kiss_fft_cpx *)sf_complexalloc(na);
    ftmp2 = (kiss_fft_cpx *)sf_complexalloc(n2 * na);
    
    /* sweep filter with zero */
    filter_sweep(1);
    filter_sweep(2);
}

/* filter sweeping */
void filter_sweep(int FILTER_CASE)
{
    int ia;
    if (FILTER_CASE == 1) {
	
	for (ia = 0; ia < na; ia++) {
	    ftmp1[ia] = zero;
	} /* initialize ftmp1 */
    } else if (FILTER_CASE == 2) {
	
	for (ia = 0; ia < n2 * na; ia++) {
	    ftmp2[ia] = zero;
	} /* initialize ftmp2 */
    } else {
	sf_error("Pick wrong case!");
    }
}

/* pick data stream for calculation */
void data_stream_pick(kiss_fft_cpx *inp, kiss_fft_cpx *stream, int i2, int i1)
{
    int ia, pt;
    
    for (ia = 0; ia < na; ia++) {
	pt = ((BOUND_SIZE + i2) - (1 + ia)) * nw + i1;
	stream[ia] = inp[pt];
    }
}

/* forward fft transform */
void forw_fft(float *dd, kiss_fft_cpx *ff, int fnt, int fnw)
{
    int i1, i2;
    float *fft_dtmp, shift, fdw;
    kiss_fft_cpx *ctrace, ce;
    kiss_fftr_cfg cfg;
    
    fdw = 1. / (fnt * d1);
    fft_dtmp = sf_floatalloc(fnt);
    ctrace = (kiss_fft_cpx *)sf_complexalloc(fnw);
    cfg = kiss_fftr_alloc(fnt, 0, NULL, NULL);
    
    /* fft1 t - >f */
    for (i2 = 0; i2 < n23; i2++) {
	
	for (i1 = 0; i1 < fnt; i1++) {
	    if (i1 < n1) {
		fft_dtmp[i1] = dd[i2 * n1 + i1];
	    } else {
		fft_dtmp[i1] = 0.; /* padding with zero */
	    }
	}
	
	kiss_fftr(cfg, fft_dtmp, ctrace);
	
	for (i1 = 0; i1 < fnw; i1++) {
	    ff[i2 * fnw + i1] = ctrace[i1]; /* centerting  */
	}
	
	if (0. != o1) {
	    for (i1 = 0; i1 < fnw; i1++) {
		shift = -2.0 * SF_PI * i1 * fdw * o1;
		ce.r = cosf(shift);
		ce.i = sinf(shift);
		ff[i2 * fnw + i1] = sf_cmul(ff[i2 * fnw + i1], ce);
	    }
	}
    }
}

/* backward fft transform */
void back_fft(kiss_fft_cpx *ff, float *dd, int bnt, int bnw)
{
    int i1, i2;
    float *ifft_dtmp, shift, bdw, wght;
    kiss_fft_cpx *ictrace, ce;
    kiss_fftr_cfg icfg;
    
    bdw = 1. / (bnt * d1);
    wght = 1.0f / bnt;
    ifft_dtmp = sf_floatalloc(bnt);
    ictrace = (kiss_fft_cpx *)sf_complexalloc(bnw);
    icfg = kiss_fftr_alloc(bnt, 1, NULL, NULL);
    
    for (i2 = 0; i2 < n23; i2++) {
	if (0. != o1) {
	    for (i1 = 0; i1 < bnw; i1++) {
		shift = +2.0 * SF_PI * i1 * bdw * o1;
		ce.r = cosf(shift);
		ce.i = sinf(shift);
		ff[i2 * bnw + i1] = sf_cmul(ff[i2 * bnw + i1], ce);
	    }
	}
	
	for (i1 = 0; i1 < bnw; i1++) {
	    ictrace[i1] = ff[i2 * bnw + i1]; /* centerting  */
	}
	
	kiss_fftri(icfg, ictrace, ifft_dtmp);
	
	for (i1 = 0; i1 < n1; i1++) {
	    dd[i2 * n1 + i1] = ifft_dtmp[i1] * wght;
	}
    }
}

/* add boundary and output conjugate */
void bound_add(kiss_fft_cpx *data, kiss_fft_cpx *ext_data,
               kiss_fft_cpx *conj_ext_data, int bn1, int bn2, int bns)
{
    int i, i1 = 0, i2 = 0;
    int bn2ns = bn2 + bns;
    if (VERB)
	sf_warning("origin data size n1=%d,n2=%d", bn1, bn2);
    if (na > n1)
	sf_error("Crossing the line in n1 axis!");
    
    for (i = 0; i < bn2ns; i++) {
	ext_data[i] = zero;
	conj_ext_data[i] = zero;
    } /* initialize tmp and ctmp */
    
    for (i1 = 0; i1 < bn1; i1++) {
	
	for (i2 = bns; i2 < bn2ns; i2++) {
	    ext_data[i2 * bn1 + i1] = data[(i2 - bns) * bn1 + i1];
	} /* cp origin data */
	
	for (i2 = 0; i2 < bns; i2++) {
	    ext_data[i2 * bn1 + i1] = data[(bns - i2) * bn1 + i1];
	} /* extend boundary with na points */
    }
    
    for (i2 = 0; i2 < bn2ns; i2++) {
	for (i1 = 0; i1 < bn1; i1++) {
	    conj_ext_data[i2 * bn1 + i1] = sf_conjf(ext_data[i2 * bn1 + i1]);
	}
    } /* calculate ctmp */
}
