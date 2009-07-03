/* Rice HPCSS seismic modeling and migration. */

/*************************************************************************

Copyright Rice University, 2009.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/

/* Modified for distribution with Madagascar */

#include <rsf.h>

#include "wavefun.h"

#ifndef _sf_wavefun_h

typedef struct {
    int nz;          /* number of gridpoints in z */
    int nx;          /* number of gridpoints in x */
    int nt;          /* number of time steps */
    int nm;          /* number of time steps to skip between movie frames
			(<=0 for no movie) */
    float dz;        /* step in z */
    float dx;        /* step in x */
    float dt;        /* step in t */
    float freq;      /* center frequency of Ricker wavelet */
    int isz;         /* source depth, in units of dz */
    int isxbeg;      /* far left source x coord in units of dx */
    int isxend;      /* far right source x coord in units of dx */
    int iskip;       /* interval between sources in units of dx */
    int igz;         /* recvr depth, in units of dz */
    int igxbeg;      /* far left receiver x coord in units of dx */
    int igxend;      /* far right receiver x coord in units of dx */
    int imbeg;       /* midpoint begin */
    int imend;       /* midpoint end */
    int imskip;      /* midpoint skip */
    int ihmax;       /* offset radius, units of dx */
    sf_file vfile;    /* velocity file - defines space grid */
    sf_file tfile;    /* trace output file */
    sf_file mfile;    /* movie file */
    sf_file rmfile;   /* receiver movie file */
    sf_file imfile;   /* image file */
} WINFO;
/*^*/

#endif

void getinputs(WINFO * wi) {
    /*< get input parameters >*/

    if (NULL != sf_getstring("velocity")) {
	/* velocity file - defines space grid */
	wi->vfile = sf_input("velocity");
	if (!sf_histint(wi->vfile,"n1",&(wi->nz)) ||
	    !sf_histint(wi->vfile,"n2",&(wi->nx)) ||
	    !sf_histfloat(wi->vfile,"d1",&(wi->dz))||
	    !sf_histfloat(wi->vfile,"d2",&(wi->dx)))
	    sf_error("Need to specify n1=, n2=, d1=, d2= in velocity");
    } else {
	wi->vfile = NULL;
	if (!sf_getint("nz",&(wi->nz))) sf_error("Need nz=");
	/* number of gridpoints in z */
	if (!sf_getint("nx",&(wi->nx))) sf_error("Need nx=");
	/* number of gridpoints in x */
	if (!sf_getfloat("dz",&(wi->dz))) sf_error("Need dz=");
	/* step in z */
	if (!sf_getfloat("dx",&(wi->dx))) sf_error("Need dx=");
	/* step in x */
    }
    if (NULL != sf_getstring("source")) {
	/* source movie file */
	wi->mfile = sf_output("source");
    } else {
	wi->mfile = NULL;
    }
    if (NULL != sf_getstring("trace")) {
	/* trace output file */
	wi->tfile = sf_output("trace");
    } else {
	wi->tfile = NULL;
    }
    if (NULL != sf_getstring("receiver")) {
	/* receiver movie file */
	wi->tfile = sf_output("receiver");
    } else {
	wi->tfile = NULL;
    }
    if (NULL != sf_getstring("image")) {
	wi->imfile= sf_output("image");
	/* image file */
    } else {
	wi->imfile= NULL;
    }

    if (!sf_getint("nt",&(wi->nt))) sf_error("Need nt=");
    /* number of time steps */
    if (!sf_getfloat("dt",&(wi->dt))) sf_error("Need dt=");
    /* step in t */

  if (!sf_getint("nm",&(wi->nm))) wi->nm=0;
  /* number of time steps to skip between movie frames
     (<=0 for no movie) */

  if (!sf_getint("isz",&(wi->isz))) wi->isz=1;
  /* source depth, in units of dz */

  if (!sf_getint("isxbeg",&(wi->isxbeg))) wi->isxbeg=(wi->nx)/2;
  /* far left source x coord in units of dx */

  if (!sf_getint("isxend",&(wi->isxend))) wi->isxend=(wi->nx)/2;
  /* far right source x coord in units of dx */

  if (!sf_getint("iskip",&(wi->iskip))) wi->iskip=1;
  /* interval between sources in units of dx */

  if (!sf_getint("igz",&(wi->igz))) wi->igz=1;
  /* recvr depth, in units of dz */

  if (!sf_getint("igxbeg",&(wi->igxbeg))) wi->igxbeg=1;
  /* far left receiver x coord in units of dx */

  if (!sf_getint("igxend",&(wi->igxend))) wi->igxend=0;
  /* far right receiver x coord in units of dx */

  if (!sf_getfloat("fpeak",&(wi->freq))) wi->freq=0.01;
  /* center frequency of Ricker wavelet */
  
  /* default is zero offset */

  if (!sf_getint("ihmax",&(wi->ihmax))) wi->ihmax=0;
  /* offset radius, units of dx */

  if (!sf_getint("imbeg",&(wi->imbeg))) wi->imbeg=wi->ihmax;
  /* midpoint begin */

  if (!sf_getint("imend",&(wi->imend))) wi->imend=wi->nx-wi->ihmax-1;
  /* midpoint end */

  if (!sf_getint("imskip",&(wi->imskip))) wi->imskip=1;
  /* midpoint skip */

  /* sanity-check */
  if (wi->iskip<1 || wi->imskip<1) 
      sf_error("either source or midpoint skip is "
	       "nonpositive - ABORT");

  if (wi->imbeg<wi->ihmax) 
      sf_error("first midpoint located within "
	       "offset radius of domain boundary - ABORT");

  if (wi->imend>wi->nx-wi->ihmax-1) 
      sf_error("first midpoint located within "
	       "offset radius of domain boundary - ABORT");

  if (wi->nx<1||wi->nz<1) 
      sf_error("number of spatial samples input is "
	       "nonpositive - ABORT");

  if (wi->isxbeg<1 || wi->isxend>wi->nx-2 ||
      wi->isz<1 || wi->isz>wi->nz-2) {
      sf_warning("source indices isz=%d, isxbeg=%d isxend=%d place source",
		 wi->isz,wi->isxbeg,wi->isxend);
      sf_error("outside of wavefield update region [1,%d] x [1,%d]",
	       wi->nz-1,wi->nx-1);
  }

  if (wi->igz<1 || wi->igz>wi->nz-2 ||
      wi->igxbeg<1 || wi->igxbeg > wi->nx-2 ||
      wi->igxend<1 || wi->igxend > wi->nx-2) {
      sf_warning("receiver depth index igz=%d or cable endpoint",wi->igz);
      sf_warning("indices igxbeg=%d, igxend=%d place receivers outside of",
		 wi->igxbeg,wi->igxend);
      sf_error("wavefield update region [1,%d] x [1,%d]",
	       wi->nz-1,wi->nx-1);
  }
}

void fassign(float * a, float c, int n) 
/*< assign constant to array - presumes memory already allocated >*/
{
  int i;
  for (i=0;i<n;i++) a[i]=c;
}

void fzeros(float * a, int n) 
/*< assign zero to array - presumes memory already allocated >*/
{
  int i;
  for (i=0;i<n;i++) a[i]=0.0;
}

void fsquare(float * a, int n) 
/*< square array elements - presumes memory already allocated
  and initialized >*/
{
  int i;
  for (i=0;i<n;i++) a[i]=a[i]*a[i];
}

float fgetmax(float * a, int n) 
/*< get max val from array >*/
{
  int i;
  float m=-SF_HUGE;
  for (i=0;i<n;i++) m=fmaxf(m,a[i]);
  return m;
}

float fgetmin(float * a, int n) 
/*< get min val from array >*/
{
  int i;
  float m=SF_HUGE;
  for (i=0;i<n;i++) m=fminf(m,a[i]);
  return m;
}

float fgetrick(float t, float fpeak) 
/*< return value at t of causal Ricker wavelet with peak frequency f >*/    
{
    float st=SF_PI*fpeak*(t-(1.2/fpeak));
    st*=st;
    return (1-2*st)*exp(-st);
}
