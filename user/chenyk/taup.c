/* Tau-p tranform unitlity */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
#include "taup.h"
#define LTABLE 8
#define NTABLE 513

float findmax(int ny, int nx, kiss_fft_cpx **c) 
/*< find the max value >*/

{
  int iy,ix; float smag,cmax=0;
  for (iy=0; iy<ny; iy++)
  for (ix=0; ix<nx; ix++)
    if (cmax<(smag=sf_cabsf(c[iy][ix]))) cmax = smag;
  return cmax;
}

void slant(float dt,int nt, int nx,float *gx, float pxmin, float dpx, int npx, float **data, float **taup, float **semb, float *sembp, float *sembpMax, float *fshift, float *fden) 
/*< taup transform 
Direct mapping of seismic data to the domain of intercept
time and ray parameter -A plane-wave decomposition, Paul L. Stoffa, Geophysics, 1981 >*/
{
  int i,it,ix,ipx; float px,pxy,*fnum,snum,sden;
  for (ipx=0; ipx<npx; ipx++) {
    fnum = taup[ipx];
    memset(&fnum[0],0,nt*sizeof(float)); memset(&fden[0],0,nt*sizeof(float));
    px = pxmin+ipx*dpx;

    for (ix=0; ix<nx; ix++) {
      pxy = px*gx[ix];

      shfs8r(dt,nt,0.0,data[ix],0.0,0.0,nt,pxy,fshift);

      for (it=0; it<nt; it++) {fnum[it] += fshift[it]; fden[it] += fshift[it]*fshift[it]; }// slant stack Stoffa Buhl Geophysics 81 (9)
    }
    for (it=0;it<nt;it++) fden[it] *= nx; // semblance Stoffa Buhl 1981 (10)
    for (it=5;it<nt-5;it++) {
      snum = sden = 0;
      for (i=it-5; i<it+5; i++) snum += fnum[i]*fnum[i], sden += fden[i];
      if (sden) semb[ipx][it] = snum/sden; // semblance Stoffa Buhl 1981 (10) 
      else      semb[ipx][it] = 0.; 
    }

    for (it=0; it<nt; it++) taup[ipx][it] = fnum[it];
  }

// max semblance along taup 
  for (sembpMax[0]=0,ipx=0; ipx<npx; ipx++) {
    sembp[ipx]=semb[ipx][0];
    for (it=0; it<nt; it++) 
      if (sembp[ipx]<semb[ipx][it]) sembp[ipx]=semb[ipx][it];
    sembpMax[0] = fmaxf(sembpMax[0],sembp[ipx]);
  }
}

void shfs8r (float dx, int nxin, float fxin, float yin[], 
	float yinl, float yinr, int nxout, float fxout, float yout[])
/*< shift a uniformly-sampled real function via a table of 8-coeff. sinc approximations. >*/
{
	static float table[NTABLE][LTABLE];
	static int tabled=0;
	int jtable,ishift,ktable,ixout,itable,jshift,ilo,ihi;
	float frac,shift,dshift,tablei,youti;

	/* tabulate sinc interpolation coefficients if not already tabulated */
	if (!tabled) {
		for (jtable=1; jtable<NTABLE-1; jtable++) {
			frac = (float)jtable/(float)(NTABLE-1);
			mksinc(frac,LTABLE,&table[jtable][0]);
		}
		for (jtable=0; jtable<LTABLE; jtable++) {
			table[0][jtable] = 0.0;
			table[NTABLE-1][jtable] = 0.0;
		}
		table[0][LTABLE/2-1] = 1.0;
		table[NTABLE-1][LTABLE/2] = 1.0;
		tabled = 1;
	}

	/* determine most appropriate set of tabulated coefficients */
	shift = (fxout-fxin)/dx;
	ishift = (int)shift;
	if (shift<0.0) ishift--;
	dshift = shift-(float)ishift;
	ktable = round(dshift*(NTABLE-1));

	/* convolve sinc approximation with input function */
	for (ixout=0; ixout<nxout; ixout++)
		yout[ixout] = 0.0;
	for (itable=0; itable<LTABLE; itable++) {
		tablei = table[ktable][itable];
		jshift = ishift+itable-LTABLE/2+1;
		for (ixout=0,youti=yinl*tablei; ixout<-jshift; ixout++)
			yout[ixout] += youti;
		ilo = fmaxf(0,-jshift);
		ihi = fminf(nxout-1,nxin-jshift-1);
		for (ixout=ilo; ixout<=ihi; ixout++)
			yout[ixout] += yin[ixout+jshift]*tablei;
		for (ixout=nxin-jshift,youti=yinr*tablei; ixout<nxout; ixout++)
			yout[ixout] += youti;
	}
}

void mksinc (float d, int lsinc, float sinc[])
/*< Compute least-squares optimal sinc interpolation coefficients. >*/
{
	int j;
	double s[20],a[20],c[20],work[20],fmax;

	/* compute auto-correlation and cross-correlation arrays */
	fmax = 0.066+0.265*log((double)lsinc);
	fmax = (fmax<1.0)?fmax:1.0;
	for (j=0; j<lsinc; j++) {
		a[j] = dsinc(fmax*j);
		c[j] = dsinc(fmax*(lsinc/2-j-1+d));
	}

	/* solve symmetric Toeplitz system for the sinc approximation */
	stoepd(lsinc,a,c,s,work);
	for (j=0; j<lsinc; j++)
		sinc[j] = s[j];
}

double dsinc (double x)
/*< Return sinc(x) = sin(PI*x)/(PI*x) (double version) >*/
{
	double pix;

	if (x==0.0) {
		return 1.0;
	} else {
		pix = SF_PI*x;
		return sin(pix)/pix;
	}
}

void stoepd (int n, double r[], double g[], double f[], double a[])
/*< Solve a symmetric Toeplitz linear system of equations Rf=g for f >*/
{
	int i,j;
	double v,e,c,w,bot;

	if (r[0] == 0.0) return;

	a[0] = 1.0;
	v = r[0];
	f[0] = g[0]/r[0];

	for (j=1; j<n; j++) {
		
		/* solve Ra=v as in Claerbout, FGDP, p. 57 */
		a[j] = 0.0;
		f[j] = 0.0;
		for (i=0,e=0.0; i<j; i++)
			e += a[i]*r[j-i];
		c = e/v;
		v -= c*e;
		for (i=0; i<=j/2; i++) {
			bot = a[j-i]-c*a[i];
			a[i] -= c*a[j-i];
			a[j-i] = bot;
		}

		/* use a and v above to get f[i], i = 0,1,2,...,j */
		for (i=0,w=0.0; i<j; i++)
			w += f[i]*r[j-i];
		c = (w-g[j])/v;
		for (i=0; i<=j; i++)
			f[i] -= c*a[j-i];
	}
}


