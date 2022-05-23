/* 3D missing data interpolation using f-x-y streaming prediction filter. */
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

/* add boundary and output temporal data */
void bound_add(kiss_fft_cpx *data, kiss_fft_cpx *data_tmp,
               kiss_fft_cpx *data_ctmp, int bn1, int bn2, int bn3, int bns1,
               int bns2);
/* forward FFT */
void forw_fft(float *dd, kiss_fft_cpx *ff, int fnt, int fnw);

/* backward FFT */
void back_fft(kiss_fft_cpx *ff, float *dd, int bnt, int bnw);

/* filter memory allocate and initialize */
void filter_init();

/* filter sweeping */
void filter_sweep(int flt_case);

/* assign value to filter */
void filter_assign(kiss_fft_cpx *flt_out, kiss_fft_cpx *flt_in);

/* judge line state */
int line_check(int NUM);

/* update temp data when interpolation */
void temp_update(kiss_fft_cpx *temp, kiss_fft_cpx *ctmp, kiss_fft_cpx *data,
                 int i1, int i2, int i3);

/* pick data stream for calculation */
void data_stream_pick(kiss_fft_cpx *inp, kiss_fft_cpx *stream, int i1, int i2,
                      int i3);

/* control calculation loop */
void loop_control(kiss_fft_cpx *inp, kiss_fft_cpx *ext_data,
                  kiss_fft_cpx *conj_ext_data);

/* filter estimate and data interpolation */
kiss_fft_cpx filter_estimete(kiss_fft_cpx data, kiss_fft_cpx *stream,
                             kiss_fft_cpx *conj_stream, int CONSTRAINT_TYPE,
                             int MASK_TYPE, int i2, int i3);

bool VERB; /* verbosity flag */
int FUNC_TYPE, FILTER_TYPE;
int BOUND_SIZE1, BOUND_SIZE2;

float o1, d1;
int n1, nw, n2, n3, n23;
int na1, na2, naa1, naa2, naa12;

int *mask;
kiss_fft_cpx *ftmp1, *ftmp2, *ftmp3, *filter;
kiss_fft_cpx lambda, lambdax, lambday, lambdaf;

kiss_fft_cpx zero;

int main(int argc, char *argv[])
{
    
    int nt;
    float *dinp, *dout;
    kiss_fft_cpx *dfft, *temp, *ctmp;
    
/*  sf_axis a1, a2, a3; */
    
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
    /* define 0+0i for convenient */
    
    if (!sf_histint(inp, "n1", &n1))
	sf_error("No n1= in input");
    if (!sf_histint(inp, "n2", &n2))
	sf_error("No n2= in input");
    if (!sf_histint(inp, "n3", &n3))
	n3 = 1;
    n23 = n2 * n3;
    if (VERB)
	sf_warning("n1=%d n23=%d", n1, n23);
    
    if (!sf_histfloat(inp, "d1", &d1))
	d1 = 1.;
    if (!sf_histfloat(inp, "o1", &o1))
	o1 = 0.;
    
    /* initialize FFT setting */
    nt = 2 * kiss_fft_next_fast_size((n1 + 1) / 2);
    if (nt % 2)
	nt++;
    nw = nt / 2 + 1;
    
    lambda = lambdax = lambday = lambdaf = zero;
    if (!sf_getfloat("lambdax", &lambdax.r))
	sf_error("Need lambdax=");
    /* lambdax.r */
    if (!sf_getfloat("lambday", &lambday.r))
	lambday = zero;
    /* lambday.r */
    if (!sf_getfloat("lambdaf", &lambdaf.r))
	lambdaf = zero;
    /* lambdaf.r */
    lambdax = sf_cmul(lambdax, lambdax); /* lambdax^2 */
    lambday = sf_cmul(lambday, lambday); /* lambday^2 */
    lambdaf = sf_cmul(lambdaf, lambdaf); /* lambdaf^2 */
    
    if (VERB) {
	sf_warning("lambdax^2 %0.10f %0.10fi", lambdax.r, lambdax.i);
	sf_warning("lambday^2 %0.10f %0.10fi", lambday.r, lambday.i);
	sf_warning("lambdaf^2 %0.10f %0.10fi", lambdaf.r, lambdaf.i);
    }
    
    if (!sf_getint("na1", &na1))
	sf_error("Need na1=");
    if (!sf_getint("na2", &na2))
	sf_error("Need na2=");
    
    if (FUNC_TYPE == 1) {
	naa1 = 2 * na1 + 1;
	naa2 = 2 * na2 + 1;
    } else if (FUNC_TYPE == 2) {
	naa1 = 2 * na1 + 1; /* size along x axis */
	naa2 = na2;         /* size along y axis */
    } else {
	sf_error("Wrong Function Type!");
    }
    
    if (VERB)
	sf_warning("Filter size %d %d", naa1, naa2);
    /* PEF filter size */

    BOUND_SIZE1 = na1;
    BOUND_SIZE2 = na2;
    naa12 = naa1 * naa2;
    
    /* memory allocate */
    dinp = sf_floatalloc(n1 * n23);
    dout = sf_floatalloc(n1 * n23);
    
    dfft = (kiss_fft_cpx *)sf_complexalloc(nw * n23);
    temp = (kiss_fft_cpx *)sf_complexalloc(nw * (n2 + BOUND_SIZE1 * 2) *
					   (n3 + BOUND_SIZE2 * 2));
    ctmp = (kiss_fft_cpx *)sf_complexalloc(nw * (n2 + BOUND_SIZE1 * 2) *
					   (n3 + BOUND_SIZE2 * 2));
    if (VERB)
	sf_warning("Allocating memory done.");
    /* memory allocate */
    
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
    free(dinp);
    if (VERB)
	sf_warning("FFT transform with nft done.");
    
    /* add boundary for dfft */
    bound_add(dfft, temp, ctmp, nw, n2, n3, BOUND_SIZE1, BOUND_SIZE2);
    if (VERB)
	sf_warning("Adding boundary for dfft done.");
    
    /* denoise / interpolation */
    loop_control(dfft, temp, ctmp);
    free(temp);
    free(ctmp);
    
    /* fft-1 */
    back_fft(dfft, dout, nt, nw);
    if (VERB)
	sf_warning("Inverse FFT transform done.");
    
    sf_floatwrite(dout, n1 * n23, out);
    /* output denoise/recover data */
    
    free(dfft);
    free(dout);
}

/* control calculation loop */
void loop_control(kiss_fft_cpx *inp, kiss_fft_cpx *ext_data,
                  kiss_fft_cpx *conj_ext_data)
{
    int li1, li2, li3, LINE_STATE = 0;
    int data_pt, mask_pt;
    kiss_fft_cpx *stream, *conj_stream;
    
    stream = (kiss_fft_cpx *)sf_complexalloc(naa12);
    conj_stream = (kiss_fft_cpx *)sf_complexalloc(naa12);
    
    filter_init();
    /* initialize filters */
    
    for (li1 = 0; li1 < nw; li1++) { /* loop over f axis */
	
	filter_sweep(2);
	filter_sweep(1);
	filter_assign(ftmp1, ftmp3);
	/* initialize filter at begining of each frequency slice */
	
	for (li3 = 0; li3 < n3; li3++) { /* loop over y axis */
	    
	    LINE_STATE = line_check(li3);
	    
	    if (LINE_STATE == 0) {
		
		for (li2 = 0; li2 < n2; li2++) { /* forward order */
		    
		    data_pt = li3 * n2 * nw + li2 * nw + li1;
		    mask_pt = li3 * n2 + li2;
		    
		    if (VERB) {
			sf_warning("i1=%3d, i2=%3d, i3=%3d, mask=%d;", 
				   li1, li2, li3, mask[mask_pt]);
		    }
		    
		    data_stream_pick(ext_data, stream, li1, li2, li3);
		    data_stream_pick(conj_ext_data, conj_stream, li1, li2, li3);
		    
		    if (li1 == 0 && li3 == 0) {
			inp[data_pt] = filter_estimete(inp[data_pt],
						       stream, conj_stream,
						       123, mask[mask_pt],
						       li2, li3);
		    } else if (li1 == 0 && li3 != 0) {
			inp[data_pt] = filter_estimete(inp[data_pt], 
						       stream, conj_stream,
						       123, mask[mask_pt], 
						       li2, li3);
		    } else if (li1 != 0 && li3 == 0) {
			inp[data_pt] = filter_estimete(inp[data_pt], 
						       stream, conj_stream,
						       123, mask[mask_pt], 
						       li2, li3);
		    } else if (li1 != 0 && li3 != 0) {
			inp[data_pt] = filter_estimete(inp[data_pt], stream, 
						       conj_stream, 123, 
						       mask[mask_pt], li2, li3);
		    } else {
			sf_error("Need constrains to compute!");
		    }
		    
		    if (mask[mask_pt] == 0)
			temp_update(ext_data, conj_ext_data, inp, 
				    li1, li2, li3);
		}
	    } else {

		for (li2 = n2 - 1; li2 > -1; li2--) { /* backward order */
		    
		    data_pt = li3 * n2 * nw + li2 * nw + li1;
		    mask_pt = li3 * n2 + li2;
		    
		    if (VERB) {
			sf_warning("i1=%3d, i2=%3d, i3=%3d, mask=%d;", 
				   li1, li2, li3,
				   mask[mask_pt]);
		    }
		    
		    data_stream_pick(ext_data, stream, li1, li2, li3);
		    data_stream_pick(conj_ext_data, conj_stream, li1, li2, li3);
		    
		    if (li1 == 0) {
			inp[data_pt] = filter_estimete(inp[data_pt], stream,
						       conj_stream, 123, 
						       mask[mask_pt], li2, li3);
		    } else {
			inp[data_pt] = filter_estimete(inp[data_pt], stream, 
						       conj_stream, 123, 
						       mask[mask_pt], li2, li3);
		    }
		    
		    if (mask[mask_pt] == 0)
			temp_update(ext_data, conj_ext_data, inp, 
				    li1, li2, li3);
		}
	    }
	    
	} /* over y */
    }   /* loop f */
    if (VERB)
	sf_warning(".");
    
    free(ftmp1);
    free(ftmp2);
    free(ftmp3);
    free(mask);
}

/* filter estimate and data interpolation */
kiss_fft_cpx filter_estimete(kiss_fft_cpx data, kiss_fft_cpx *stream,
                             kiss_fft_cpx *conj_stream, int CONSTRAINT_TYPE,
                             int MASK_TYPE, int i2, int i3)
{
    int ia;
    kiss_fft_cpx sTs, sTa, sn, rn, res, A, B, C, D;
    
    if (CONSTRAINT_TYPE == 1) {
	lambda = lambdax;
    } else if (CONSTRAINT_TYPE == 12) {
	lambda = sf_cadd(lambdax, lambday);
    } else if (CONSTRAINT_TYPE == 13) {
	lambda = sf_cadd(lambdax, lambdaf);
    } else if (CONSTRAINT_TYPE == 123) {
	lambda = sf_cadd(lambdax, sf_cadd(lambday, lambdaf));
    }
    
    sTs = zero;
    for (ia = 0; ia < naa12; ia++) {
	sTs = sf_cadd(sTs, sf_cmul(stream[ia], conj_stream[ia]));
    } /* update sTs */
    
    sTa = zero;
    for (ia = 0; ia < naa12; ia++) {
	if (CONSTRAINT_TYPE == 1) {
	    A = sf_cmul(lambdax, ftmp1[ia]);
	} else if (CONSTRAINT_TYPE == 12) {
	    A = sf_cadd(sf_cmul(lambdax, ftmp1[ia]),
			sf_cmul(lambday, ftmp2[i2 * naa12 + ia]));
	} else if (CONSTRAINT_TYPE == 13) {
	    A = sf_cadd(sf_cmul(lambdax, ftmp1[ia]),
			sf_cmul(lambdaf, ftmp3[(i3 * n2 + i2) * naa12 + ia]));
	} else if (CONSTRAINT_TYPE == 123) {
	    A = sf_cadd(
		sf_cmul(lambdax, ftmp1[ia]),
		sf_cadd(sf_cmul(lambday, ftmp2[i2 * naa12 + ia]),
			sf_cmul(lambdaf, ftmp3[(i3 * n2 + i2) * naa12 + ia])));
	} else {
	    sf_error(("Error of constraint type!"));
	}
	
	B = sf_cmul(stream[ia], sf_cdiv(A, lambda));
	sTa = sf_cadd(sTa, B);
    } /* update sTa */
    
    if (FUNC_TYPE == 1) { /* denoise */
	
	rn = sf_cdiv(sf_cadd(data, sTa), sf_cadd(lambda, sTs));
	res = sf_cmul(lambda, rn);
	sn = sf_csub(data, res);
    } else if (FUNC_TYPE == 2) { /* interpolation */
	
	if (MASK_TYPE) { /* calculate filter */
	    sn = data;
	    rn = sf_cdiv(sf_cadd(data, sTa), sf_cadd(lambda, sTs));
	} else { /* interpolate gap */
	    sn = sf_csub(zero, sTa);
	    rn = zero;
	} /* calculate residual */
    } else {
	sf_error("Wrong type been chosed!");
    }
    
    for (ia = 0; ia < naa12; ia++) {
	
	if (CONSTRAINT_TYPE == 1) {
	    
	    C = sf_cmul(lambdax, ftmp1[ia]);
	} else if (CONSTRAINT_TYPE == 12) {
	    C = sf_cadd(sf_cmul(lambdax, ftmp1[ia]),
			sf_cmul(lambday, ftmp2[i2 * naa12 + ia]));
	} else if (CONSTRAINT_TYPE == 13) {
	    C = sf_cadd(sf_cmul(lambdax, ftmp1[ia]),
			sf_cmul(lambdaf, ftmp3[(i3 * n2 + i2) * naa12 + ia]));
	} else if (CONSTRAINT_TYPE == 123) {
	    C = sf_cadd(
		sf_cmul(lambdax, ftmp1[ia]),
		sf_cadd(sf_cmul(lambday, ftmp2[i2 * naa12 + ia]),
			sf_cmul(lambdaf, ftmp3[(i3 * n2 + i2) * naa12 + ia])));
	}
	
	D = sf_cdiv(C, lambda);
	ftmp1[ia] = sf_csub(D, sf_cmul(rn, conj_stream[ia]));
	
	/* if (lambday.r > 0.) */
	ftmp2[i2 * naa12 + ia] = ftmp1[ia];
	/* if (lambdaf.r > 0.) */
	ftmp3[(i3 * n2 + i2) * naa12 + ia] = ftmp1[ia];
    } /* update ftmp1, ftmp2, filter */
    
    return sn;
}

/* initialize filter */
void filter_init()
{
    ftmp1 = (kiss_fft_cpx *)sf_complexalloc(naa12);
    ftmp2 = (kiss_fft_cpx *)sf_complexalloc(naa12 * n2);
    ftmp3 = (kiss_fft_cpx *)sf_complexalloc(naa12 * n2 * n3);
    /* memory allocation */
    
    filter_sweep(1);
    filter_sweep(2);
    filter_sweep(3);
    /* sweep filter with zero */
}

/* sweeping filter */
void filter_sweep(int flt_case)
{
    int ia;
    
    if (flt_case == 1) {

	for (ia = 0; ia < naa12; ia++) {
	    ftmp1[ia] = zero;
	} /* initialize ftmp1 */
    } else if (flt_case == 2) {
	for (ia = 0; ia < naa12 * n2; ia++) {
	    ftmp2[ia] = zero;
	} /* initialize ftmp2 */
    } else if (flt_case == 3) {
	for (ia = 0; ia < naa12 * n2 * n3; ia++) {
	    ftmp3[ia] = zero;
	} /* initialize ftmp3 */
    } else {
	sf_error("Pick wrong case!");
    }
}

/* pick data stream for calculation */
void data_stream_pick(kiss_fft_cpx *inp, kiss_fft_cpx *stream, int pi1, int pi2,
                      int pi3)
{
    
    int i2, i3, pn1, pn2, pt;
    /* base = 1; */
    pn2 = n2 + BOUND_SIZE1 * 2;
    pn1 = nw;
    
    for (i3 = 0; i3 < naa2; i3++) { /* y axis */
	for (i2 = 0; i2 < naa1; i2++) { /* x axis */
	    
	    pt = (i3 + pi3) * (pn2 * pn1) + (i2 + pi2) * pn1 + pi1;
	    stream[i3 * naa1 + i2] = inp[pt];
	    if (i3 == na2 && i2 == na1)
		stream[i3 * naa1 + i2] = zero;
	}
    }
}

/* add boundary and output temporal data */
void bound_add(kiss_fft_cpx *data, kiss_fft_cpx *data_tmp,
               kiss_fft_cpx *data_ctmp, int bn1, int bn2, int bn3, int bns1,
               int bns2)
{
    int i, i1 = 0, i2 = 0, i3 = 0, ipt, ipc, ipd;
    int bn2ns1, bn3ns2, nbond2, nbond3;
    
    bn2ns1 = bn2 + bns1;
    bn3ns2 = bn3 + bns2;
    
    nbond2 = bn2 + bns1 * 2;
    nbond3 = bn3 + bns2 * 2;
    
    if (VERB)
	sf_warning("origin data size n1=%d,n2=%d,n2=%d", bn1, bn2, bn3);
    if (VERB)
	sf_warning("temp data size n1=%d,n2=%d,n2=%d", bn1, nbond2, nbond3);
    if (na1 + 1 > bn2)
	sf_error("Crossing the line in n1 axis!");
    if (na2 + 1 > bn3)
	sf_error("Crossing the line in n2 axis!");
    
    for (i = 0; i < bn1 * nbond2 * nbond3; i++) {
	data_tmp[i] = zero;
	data_ctmp[i] = zero;
    } /* initialize tmp and ctmp */
    
    for (i1 = 0; i1 < bn1; i1++) { /* loop over frequency axis */
	
	for (i3 = bns2; i3 < bn3ns2; i3++) { /* loop over y axis */
	    for (i2 = bns1; i2 < bn2ns1; i2++) { /* loop over x axis */
		ipt = i3 * (bn1 * nbond2) + i2 * bn1 + i1;
		ipd = (i3 - bns2) * (bn1 * bn2) + (i2 - bns1) * bn1 + i1;
		data_tmp[ipt] = data[ipd];
	    }
	} /* cp origin data */
	
	for (i3 = bns2; i3 < bn3ns2; i3++) { /* loop over y axis */
	    for (i2 = 0; i2 < bns1; i2++) { /* loop over x axis */
		ipt = i3 * (bn1 * nbond2) + i2 * bn1 + i1;
		ipd = i3 * (bn1 * nbond2) + (2 * bns1 - i2) * bn1 + i1;
		data_tmp[ipt] = data_tmp[ipd];
	    }
	    for (i2 = 0; i2 < bns1; i2++) {
		ipt = i3 * (bn1 * nbond2) + (bn2ns1 + i2) * bn1 + i1;
		ipd = i3 * (bn1 * nbond2) + (bn2ns1 - i2 - 1) * bn1 + i1;
		data_tmp[ipt] = data_tmp[ipd];
	    }
	} /* extend boundary with naa1 points at x axis */
	
	for (i3 = 0; i3 < bns2; i3++) { /* loop over y axis */
	    for (i2 = 0; i2 < bn2ns1; i2++) { /* loop over x axis */
		ipt = i3 * (bn1 * nbond2) + i2 * bn1 + i1;
		ipd = (2 * bns2 - i3) * (bn1 * nbond2) + i2 * bn1 + i1;
		data_tmp[ipt] = data_tmp[ipd];
	    }
	    for (i2 = 0; i2 < bn2ns1; i2++) { /* loop over x axis */
		ipt = (bn3ns2 + i3) * (bn1 * nbond2) + i2 * bn1 + i1;
		ipd = (bn3ns2 - i3 - 1) * (bn1 * nbond2) + i2 * bn1 + i1;
		data_tmp[ipt] = data_tmp[ipd];
	    }
	} /* extend boundary with naa2 points at y axis */
	
	for (i3 = 0; i3 < nbond3; i3++) {
	    for (i2 = 0; i2 < nbond2; i2++) {
		ipt = i3 * (bn1 * nbond2) + i2 * bn1 + i1;
		ipc = ipt;
		data_ctmp[ipc] = sf_conjf(data_tmp[ipt]);
	    }
	} /* calculate ctmp */
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

/* judge parity */
int line_check(int NUM)
{
    if (NUM % 2 == 0) {
	return 0;
    } else {
	return 1;
    }
}

/* update temp data when interpolation */
void temp_update(kiss_fft_cpx *temp, kiss_fft_cpx *ctmp, kiss_fft_cpx *data,
                 int i1, int i2, int i3)
{
    int dpos, tpos, pn2, pn1;
    
    pn2 = n2 + BOUND_SIZE1 * 2;
    pn1 = nw;
    
    dpos = i3 * n2 * nw + i2 * nw + i1;
    
    tpos = (i3 + BOUND_SIZE2) * pn2 * pn1 + (i2 + BOUND_SIZE1) * pn1 + i1;
    
    temp[tpos] = data[dpos];
    ctmp[tpos] = sf_conjf(temp[tpos]);
}

/* assign value to filter */
void filter_assign(kiss_fft_cpx *flt_out, kiss_fft_cpx *flt_in)
{
    int ia;
    
    for (ia = 0; ia < naa12; ia++) {
	flt_out[ia] = flt_in[ia];
    }
}
