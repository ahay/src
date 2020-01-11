/* f-x-y domian interpolation/denoise using 3D streaming PEF */
/*
  Copyright (C) 2017 Jilin University

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

void boundary(kiss_fft_cpx *data, kiss_fft_cpx *data_tmp,
              kiss_fft_cpx *data_ctmp,int bn1,int bn2,int bn3,
              int bns1, int bns2);
void forw_fft(float *dd,kiss_fft_cpx *ff, int fnt, int fnw);
void back_fft(kiss_fft_cpx *ff, float *dd, int bnt, int bnw);
void pick_cef(kiss_fft_cpx *inp, kiss_fft_cpx *out,
              int pn1, int pn2, int pi1, int pi2, int pi3);
void init_flt();
void sweep_flt(int flt_case);
void loop_control(kiss_fft_cpx *data, kiss_fft_cpx *data_tmp,
                  kiss_fft_cpx *data_ctp);
void est_point(kiss_fft_cpx *data, kiss_fft_cpx *data_tmp,
               kiss_fft_cpx *data_ctmp, int i1, int i2, int i3);

bool verb;      /* verbosity flag */
int type;       /* 1: denoise, 2 : interpolation */

/* input */
float *dinp, *dout;
float o1, d1, o2, d2, o3, d3;
int n1, n2, n3, n23;

int *mask;
kiss_fft_cpx  *dfft, *res;
float dw;
int nt, nw;

kiss_fft_cpx *temp, *ctemp;    /* temp store */
int bdsz2, bdsz1;

kiss_fft_cpx *ftmp1, *ftmp2, *ftmp3, *filter;
kiss_fft_cpx lambda, lambda1, lambda2, lambda3;
int na1, na2, naa1, naa2, naa12;

kiss_fft_cpx zero;

int main(int argc, char *argv[])
{
  sf_file inp, out, msk=NULL;

  sf_init(argc, argv);

  inp  = sf_input ("in");
  out  = sf_output("out");

  if (!sf_getbool("verb", &verb))  verb  = false;
  /* default=false, verbosity flag */

  if (!sf_getint("type", &type))   type  = 1;
  /* default = 1,
     1 : denoise with noncausul filter desin
     2 : interpolation with causul filter design*/

  zero.r = 0.0f;
  zero.i = 0.0f;
  /* define 0+0i for convenient */

  if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
  if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
  if (!sf_histint(inp,"n3",&n3)) n3 = 1;
  n23 = n2 * n3;
  if (verb) sf_warning("n1=%d n23=%d", n1, n23);

  if (!sf_histfloat(inp,"d1",&d1)) d1 = 1. ;
  if (!sf_histfloat(inp,"o1",&o1)) o1 = 0. ;

  if (!sf_histfloat(inp,"d2",&d2)) d2 = 1. ;
  if (!sf_histfloat(inp,"o2",&o2)) o2 = 0. ;

  if (!sf_histfloat(inp,"d3",&d3)) d3 = 1. ;
  if (!sf_histfloat(inp,"o3",&o3)) o3 = 0. ;

  nt  = 2*kiss_fft_next_fast_size((n1+1)/2);
  if (nt%2) nt++;
  nw  = nt/2+1;
  dw  = 1./(nt*d1);
  /* initialize FFT setting */

  lambda = lambda1  = lambda2 = zero;
  if (!sf_getfloat("lambda1", &lambda1.r)) sf_error("Need lambda1=");
  /* lambda1 */

  if (!sf_getfloat("lambda2", &lambda2.r)) sf_error("Need lambda2=");
  /* lambda2 */

  if (!sf_getfloat("lambda3", &lambda3.r)) sf_error("Need lambda3=");
  /* lambda3 */

  lambda1 = sf_cmul(lambda1, lambda1);    /* lambda1^2 */
  lambda2 = sf_cmul(lambda2, lambda2);    /* lambda2^2 */
  lambda3 = sf_cmul(lambda3, lambda3);    /* lambda3^2 */

  lambda = sf_cadd(lambda1, sf_cadd(lambda2, lambda3));      /* lambda */
  if (verb) {
    sf_warning("lambda %0.10f %0.10fi", lambda.r, lambda.i);
    sf_warning("lambda1^2 %0.10f %0.10fi", lambda1.r, lambda1.i);
    if (lambda2.r > 0.0f) {
      sf_warning("lambda2^2 %0.10f %0.10fi", lambda2.r, lambda2.i);
      if (lambda3.r > 0.0f)
        sf_warning("lambda3^2 %0.10f %0.10fi", lambda3.r, lambda3.i);
    } else if ( lambda3.r > 0.0f) {
      sf_error("Need lambda2=");
    }
  }

  if (!sf_getint("na1", &na1)) sf_error("Need na1=");
  /* fiter size along x axis */
  if (!sf_getint("na2", &na2)) sf_error("Need na2=");
  /* fiter size along y axis */
   if (type == 1) { /* fiter size for denoise */
    naa1 = 2*na1 + 1; /* size along x axis */
    naa2 = 2*na2 + 1; /* size along y axis */
  } else if (type == 2){ /* fiter size for interpolation */
    naa1 = 2*na1 + 1;
    naa2 = na2;
  }
  if (verb) sf_warning("Filter size %d %d", naa1, naa2);
  /* PEF filter size */
  bdsz1 = na1;
  bdsz2 = na2;
  naa12 = naa1 * naa2;
  /* boundary size */

  dinp = sf_floatalloc(n1*n23);
  dout = sf_floatalloc(n1*n23);
  dfft  = (kiss_fft_cpx*)sf_complexalloc(nw*n23);
  res   = (kiss_fft_cpx*)sf_complexalloc(nw*n23);
  /* memory allocate */

  temp  = (kiss_fft_cpx*)sf_complexalloc(nw*(n2+bdsz1*2)*(n3+bdsz2*2));
  ctemp = (kiss_fft_cpx*)sf_complexalloc(nw*(n2+bdsz1*2)*(n3+bdsz2*2));
  if (verb) sf_warning("Allocating memory done.");
  /* memory allocate */

  if (type == 2) {
    msk = sf_input ("mask");
    if (SF_INT!=sf_gettype(msk)) sf_error("Need int type in mask");
    /* if type=2, need mask for interpolation */

    mask = sf_intalloc(n23);
    /* memory allocate */

    sf_intread(mask, n23, msk);
    if (verb) sf_warning("Read mask data done");
    /* read mask data */
  }

  sf_floatread(dinp, n1*n23, inp);
  /* read data */
  if (verb) sf_warning("Reading data done.");

  forw_fft(dinp, dfft, nt, nw);
  /* fft with nft points */
  if (verb) sf_warning("FFT transform with nft done.");

  boundary(dfft, temp, ctemp, nw, n2, n3, bdsz1, bdsz2);
  /* add boundary for dfft */
  if (verb) sf_warning("Adding boundary for dfft done.");

  loop_control(dfft, temp, ctemp);
  /* denoise/interpolation */

  back_fft(dfft, dout, nt, nw);
  /* fft-1 */
  if (verb) sf_warning("Inverse FFT transform done.");

  sf_floatwrite(dout, n1*n23, out);
  /* output denoise/recover data */

}


/* forward fft transform */
void forw_fft(float *dd, kiss_fft_cpx *ff, int fnt, int fnw)
{
  int i1, i2;
  float *fft_dtmp, shift, fdw;
  kiss_fft_cpx *ctrace, ce;
  kiss_fftr_cfg cfg;

  fdw = 1./(fnt*d1);
  fft_dtmp = sf_floatalloc(fnt);
  ctrace = (kiss_fft_cpx*)sf_complexalloc(fnw);
  cfg  = kiss_fftr_alloc(fnt, 0, NULL, NULL);

  /* fft1 t - >f */
  for (i2=0; i2 < n23; i2++) {

    for (i1=0; i1 < fnt; i1++) {
      if (i1 < n1) {
        fft_dtmp[i1] = dd[i2*n1+i1];
      } else {
        fft_dtmp[i1] = 0. ;        /* padding with zero */
      }
    }

    kiss_fftr(cfg, fft_dtmp, ctrace);

    for (i1=0; i1 < fnw; i1++)	{
      ff[i2*fnw+i1] = ctrace[i1];  /* centerting  */
      /* if (verb) sf_warning("ff=%f", ctrace[i1]); */
    }

    if (0. != o1) {
      for (i1=0; i1 < fnw; i1++) {
        shift = -2.0 * SF_PI * i1 * fdw * o1;
        ce.r = cosf(shift);
        ce.i = sinf(shift);
        ff[i2*fnw+i1]=sf_cmul(ff[i2*fnw+i1],ce);
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

  bdw = 1./(bnt*d1);
  wght = 1.0f/bnt;
  ifft_dtmp = sf_floatalloc(bnt);
  ictrace = (kiss_fft_cpx*)sf_complexalloc(bnw);
  icfg  = kiss_fftr_alloc(bnt, 1, NULL, NULL);

  for (i2=0; i2 < n23; i2++) {
    if (0. != o1) {
      for (i1=0; i1 < bnw; i1++) {
        shift = +2.0 * SF_PI * i1 * bdw * o1;
        ce.r = cosf(shift);
        ce.i = sinf(shift);
        ff[i2*bnw+i1] = sf_cmul(ff[i2*bnw+i1], ce);
      }
    }

    for (i1=0; i1 < bnw; i1++)	{
     ictrace[i1] = ff[i2*bnw+i1];  /* centerting  */
    }

    kiss_fftri(icfg, ictrace, ifft_dtmp);

    for (i1=0; i1<n1; i1++) {
      dd[i2*n1+i1] = ifft_dtmp[i1]*wght;
    }
  }

}

/* add boundary and output conjugate */
void boundary(kiss_fft_cpx *data, kiss_fft_cpx *data_tmp,
              kiss_fft_cpx *data_ctmp,int bn1,int bn2,int bn3,
              int bns1, int bns2)
{
  int i, i1=0, i2=0, i3=0, ipt, ipc, ipd;
  int bn2ns1, bn3ns2, nbond2, nbond3;

  bn2ns1 = bn2 + bns1;
  bn3ns2 = bn3 + bns2;

  nbond2 = bn2 + bns1*2;
  nbond3 = bn3 + bns2*2;

  if (verb) sf_warning("origin data size n1=%d,n2=%d,n2=%d",
                       bn1, bn2, bn3);
  if (verb) sf_warning("temp data size n1=%d,n2=%d,n2=%d",
                       bn1, nbond2, nbond3);
  if (na1+1 > bn2) sf_error("Crossing the line in n1 axis!");
  if (na2+1 > bn3) sf_error("Crossing the line in n2 axis!");

  for(i=0; i < bn1*nbond2*nbond3; i++) {
    data_tmp[i]  = zero;
    data_ctmp[i] = zero;
  }/* initialize tmp and ctmp */

  for (i1=0; i1 < bn1; i1++) {  /* loop over frequency axis */

    for (i3=bns2; i3 < bn3ns2; i3++) {    /* loop over y axis */
      for (i2=bns1; i2 < bn2ns1; i2++) {  /* loop over x axis */
        ipt =     i3   *(bn1*nbond2) +     i2   *bn1 + i1;
        ipd = (i3-bns2)*  (bn1*bn2) + (i2-bns1)*bn1 + i1;
        data_tmp[ipt] = data[ipd];
      }
    }/* cp origin data */

    for (i3=bns2; i3 < bn3ns2; i3++) {   /* loop over y axis */
      for (i2=0; i2 < bns1; i2++) {      /* loop over x axis */
        ipt = i3*(bn1*nbond2) +     i2   *bn1 + i1;
        ipd = i3*(bn1*nbond2) + (2*bns1-i2)*bn1 + i1;
        data_tmp[ipt] = data_tmp[ipd];
      }
      for (i2=0; i2 < bns1; i2++) {
        ipt = i3*(bn1*nbond2) + (bn2ns1+i2)  *bn1 + i1;
        ipd = i3*(bn1*nbond2) + (bn2ns1-i2-1)*bn1 + i1;
        data_tmp[ipt] = data_tmp[ipd];
      }
    }/* extend boundary with naa1 points at x axis */

    for (i3=0; i3 < bns2; i3++) {         /* loop over y axis */
      for (i2=0; i2 < bn2ns1; i2++) {     /* loop over x axis */
        ipt =      i3    *(bn1*nbond2) + i2*bn1 + i1;
        ipd = (2*bns2-i3)*(bn1*nbond2) + i2*bn1 + i1;
        data_tmp[ipt] = data_tmp[ipd];
      }
      for (i2=0; i2 < bn2ns1; i2++) {     /* loop over x axis */
        ipt = (bn3ns2+i3)  *(bn1*nbond2) + i2*bn1 + i1;
        ipd = (bn3ns2-i3-1)*(bn1*nbond2) + i2*bn1 + i1;
        data_tmp[ipt] = data_tmp[ipd];
      }
    }/* extend boundary with naa2 points at y axis */

    for (i3=0; i3 < nbond3; i3++) {
      for (i2=0; i2 < nbond2; i2++) {
        ipt = i3*(bn1*nbond2) + i2*bn1 + i1;
        ipc = ipt;
        data_ctmp[ipc] = sf_conjf(data_tmp[ipt]);
        /* sf_warning("ctmp[%d]=%f", i1, ctmp[i2*nw2+i1].r); */
      }
    }/* calculate ctmp */
  }

}

/* initialize filter */
void init_flt()
{
  ftmp1  = (kiss_fft_cpx*)sf_complexalloc(naa12);
  ftmp2  = (kiss_fft_cpx*)sf_complexalloc(naa12*n2);
  ftmp3  = (kiss_fft_cpx*)sf_complexalloc(naa12*n2*n3);
  /* memory allocation */

  sweep_flt(1);
  sweep_flt(2);
  sweep_flt(3);
  /* sweep filter with zero */

}

/* sweeping filter */
void sweep_flt(int flt_case)
{
  int ia;

  if (flt_case == 1) {

    for (ia=0; ia < naa12; ia++) {
      ftmp1[ia] = zero;
    } /* initialize ftmp1 */

  } else if (flt_case ==2) {
    for (ia=0; ia < naa12*n2; ia++) {
      ftmp2[ia] = zero;
    }/* initialize ftmp2 */

  } else if (flt_case ==3) {
    for (ia=0; ia < naa12*n2*n3; ia++) {
      ftmp3[ia] = zero;
    }/* initialize ftmp3 */

  } else {
    sf_error("Pick wrong case!");
  }

}

/* arrange snaky looping processing order  */
void loop_control(kiss_fft_cpx *data, kiss_fft_cpx *data_tmp,
                  kiss_fft_cpx *data_ctp)
{
  int ci1, ci2, ci3;

  init_flt();
  /* initialize filters */

  for (ci1=0; ci1 < nw; ci1++) { /* loop over f axis */

    sweep_flt(1);
    sweep_flt(2);
    /* initialize filter at begining of each frequency slice */

    for (ci3=0; ci3 < n3; ci3++) {  /* loop over y axis */

      if (ci3%2 == 0) {                 /* loop over x axis */
        for (ci2=0; ci2 < n2; ci2++) {     /* forward order */

          if(verb) sf_warning("i1=%3d, i2=%3d, i3=%3d;", ci1, ci2, ci3);
          est_point(data, data_tmp, data_ctp, ci1, ci2, ci3);
          /* calculate point */
        }
      } else {
        for (ci2=n2-1; ci2 > -1; ci2-- ) { /* backward order */

          if(verb) sf_warning("i1=%3d, i2=%3d, i3=%3d;", ci1, ci2, ci3);
          est_point(data, data_tmp, data_ctp, ci1, ci2, ci3);
          /* calcluate point */
        }
      }

    }
  }
  if(verb) sf_warning(".");
}

/* streaming prediction filter */
void est_point(kiss_fft_cpx *data, kiss_fft_cpx *data_tmp,
               kiss_fft_cpx *data_ctmp, int i1, int i2, int i3)
{
  int pt, pt_sli, pt_tmp, ia;
  kiss_fft_cpx sTs, sTa, sn, rn, A, B, C, D;
  kiss_fft_cpx *tmp, *ctp;

  tmp    = (kiss_fft_cpx*)sf_complexalloc(naa12);
  ctp   = (kiss_fft_cpx*)sf_complexalloc(naa12);

  pt     = i3        *(n2*nw)           + i2        *nw + i1;
  pt_tmp = (i3+bdsz2)*((n2+bdsz1*2)*nw) + (i2+bdsz1)*nw + i1;
  pt_sli = i3        *n2                + i2;

  pick_cef( data_tmp, tmp, nw, n2+bdsz1*2, i1, i2, i3);
  pick_cef(data_ctmp, ctp, nw, n2+bdsz1*2, i1, i2, i3);

  sTs = zero;
  for (ia=0; ia < naa12; ia++) {
    sTs = sf_cadd(sTs, sf_cmul(tmp[ia], ctp[ia]));
  }
  /* update dd */

  sTa = zero;
  for (ia=0; ia < naa12; ia++) {
    A = sf_cadd(sf_cmul(lambda1, ftmp1[ia]),
                sf_cadd(sf_cmul(lambda2, ftmp2[i2*naa12+ia]),
                        sf_cmul(lambda3, ftmp3[pt_sli*naa12+ia])));

    B = sf_cmul(tmp[ia], sf_cdiv(A, lambda));
    sTa = sf_cadd(sTa, B);
  }/* update da */

  if (type == 1) { /* denoise */
    sn = data[pt];
    rn = sf_cdiv(sf_cadd(sn, sTa), sf_cadd(lambda, sTs));
    res[pt] = sf_cmul(lambda, rn);
    data[pt] = sf_csub(sn, res[pt]);
    /* calculate rn */

  } else if ( type == 2) { /* interpolation */
    if (mask[i3*n2+i2]) {
      sn = data[pt];
      rn = sf_cdiv(sf_cadd(sn, sTa), sf_cadd(lambda, sTs));
      res[pt] = sf_cmul(lambda, rn);
      /* calculate rn */
    } else {
      sn = sf_csub(zero, sTa);
      rn = zero;
      data[pt] = sn;
      /* interpolation missing point */

      data_tmp[pt_tmp]  = sn;
      data_ctmp[pt_tmp] = sf_conjf(sn);
      /* update temp of d */
    }
  } else {
    sf_error("Wrong type been chosed!");
  }

  for (ia=0; ia < naa12; ia++) {

    C = sf_cadd(sf_cmul(lambda1, ftmp1[ia]),
                sf_cadd(sf_cmul(lambda2, ftmp2[i2*naa12+ia]),
                        sf_cmul(lambda3, ftmp3[pt_sli*naa12+ia])));
    D = sf_cdiv(C, lambda);
    ftmp1[ia] = sf_csub(D, sf_cmul(rn, ctp[ia]));

    ftmp2[i2*naa12+ia]     = ftmp1[ia];
    ftmp3[pt_sli*naa12+ia] = ftmp1[ia];
    // filter[pt*naa12+ia]    = ftmp1[ia];
  }/* update ftmp1, ftmp2, filter */

  free(tmp);
  free(ctp);
}


/* pick data to esimate filter */
void pick_cef(kiss_fft_cpx *inp, kiss_fft_cpx *out,
              int pn1, int pn2, int pi1, int pi2, int pi3)
{
  int i2=0, i3, point;

  for (i3=0; i3 < naa2; i3++) {   /* y axis */
    for (i2=0; i2 < naa1; i2++) {   /* x axis */
      point = (i3+pi3) * (pn2*pn1) + (i2+pi2) * pn1 + pi1;
      out[i3*naa1+i2] = inp[point];
      if (i3==na2 && i2==na1) out[i3*naa1+i2] = zero;
    }
  }

}
