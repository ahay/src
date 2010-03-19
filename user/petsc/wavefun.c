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

#ifndef _wavefun_h

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

void getinputs(bool mod,  /* modeling or migration */ 
	       WINFO * wi /* parameters */) 
/*< get input parameters >*/
{
    wi->vfile = sf_input("in");
    /* velocity file - defines space grid */

    if (!sf_histint(wi->vfile,"n1",&(wi->nz)) ||
	!sf_histint(wi->vfile,"n2",&(wi->nx)) ||
	!sf_histfloat(wi->vfile,"d1",&(wi->dz))||
	!sf_histfloat(wi->vfile,"d2",&(wi->dx)))
	sf_error("Need to specify n1=, n2=, d1=, d2= in velocity");

    if (NULL != sf_getstring("source")) {
	/* source movie file */
	wi->mfile = sf_output("source");
    } else {
	wi->mfile = NULL;
    }

    if (mod) {
	wi->tfile = sf_output("out");
	if (!sf_getint("nt",&(wi->nt))) sf_error("Need nt=");
	/* number of time steps */
	if (!sf_getfloat("dt",&(wi->dt))) sf_error("Need dt=");
	/* step in t */	
    } else { 
	wi->tfile = sf_input("trace");
	if (!sf_histint(wi->tfile,"n1",&(wi->nt))) 
	    sf_error("Need n1= in trace");
	if (!sf_histfloat(wi->tfile,"d1",&(wi->dt))) 
	    sf_error("Need d1= in trace");
    }

    if (NULL != sf_getstring("receiver")) {
	/* receiver movie file */
	wi->rmfile = sf_output("receiver");
    } else {
	wi->rmfile = NULL;
    }

    if (mod) {
	wi->imfile= NULL;
    } else {
	wi->imfile= sf_output("out");
	/* image file */
    } 


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

    if (mod) {
	sf_putint(wi->tfile,"n1",wi->nt);
	sf_putfloat(wi->tfile,"d1",wi->dt);
	sf_putfloat(wi->tfile,"o1",0.);
	sf_putstring(wi->tfile,"label1","Time");
	sf_putstring(wi->tfile,"unit1","s");
	sf_putint(wi->tfile,"n2",wi->igxend-wi->igxbeg+1);
	sf_putfloat(wi->tfile,"d2",wi->dx);
	sf_putfloat(wi->tfile,"o2",wi->igxbeg*wi->dx);
	sf_putstring(wi->tfile,"label2","Receiver");
	sf_putint(wi->tfile,"n3",(wi->isxend-wi->isxbeg)/wi->iskip+1);
	sf_putfloat(wi->tfile,"d3",wi->iskip*wi->dx);
	sf_putfloat(wi->tfile,"o3",wi->isxbeg*wi->dx);
	sf_putstring(wi->tfile,"label3","Source");
    } 

    if (NULL != wi->mfile && wi->nm) {
	sf_putint(wi->mfile,"n3",wi->nt/wi->nm);
	sf_putfloat(wi->mfile,"o3",0.);
	sf_putfloat(wi->mfile,"d3",wi->nm*wi->dt);
	sf_putstring(wi->mfile,"label3","Time");
	sf_putstring(wi->mfile,"unit3","s");
    }

    if (NULL != wi->rmfile) {
	sf_putint(wi->rmfile,"n3",wi->nt);
	sf_putfloat(wi->rmfile,"o3",0.);
	sf_putfloat(wi->rmfile,"d3",wi->dt);
	sf_putstring(wi->rmfile,"label3","Time");
	sf_putstring(wi->rmfile,"unit3","s");
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
