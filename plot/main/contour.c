#include <math.h>

#include <rsf.h>
#include <rsfplot.h>

static void contour (float **z, int n1, int n2, float c);

int main (int argc, char* argv[])
{
    int n1, n2, n3, i3, nc0, nc, ic, n12, i1, i2;
    float **z, zi, dc, c0, zmin, zmax, *c;
    float o1, o2, d1, d2, min1, min2, max1, max2;
    bool hasc, hasdc, hasc0, transp;
    sf_file in;

    sf_init(argc,argv);
    in = sf_input("in");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;
    
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;

    if (!sf_getfloat("min1",&min1)) min1=o1;
    if (!sf_getfloat("min2",&min2)) min2=o2;
    if (!sf_getfloat("max1",&max1)) max1=o1+(n1-1)*d1;
    if (!sf_getfloat("max2",&max2)) max2=o2+(n2-1)*d2;

    if (!sf_getint("nc",&nc0)) nc=nc0=5;

    c = sf_floatalloc(nc);
    vp_plot_init(nc);

    hasc = sf_getfloats("c",c,nc);
    hasdc = sf_getfloat("dc",&dc);
    hasc0 = sf_getfloat("c0",&c0);
    
    z = sf_floatalloc2(n1,n2);

    if (!sf_getbool ("transp",&transp)) transp=false;

    for (i3=0; i3 < n3; i3++) {
	sf_read(z[0],sizeof(float),n12,in);
	
	if (!hasc) {
	    if (!hasdc || !hasc0) {
		zmin = z[0][0];
		zmax = z[0][0];
		for (i2=0; i2 < n2; i2++) {
		    for (i1=0; i1 < n1; i1++) {
			zi= z[i2][i1];
			if      (zi < zmin) zmin=zi;
			else if (zi > zmax) zmax=zi;
		    }
		}
		if (hasdc) {
		    for (c0 = floorf(zmin/dc) * dc - dc; c0 < zmin; c0 += dc) ;
		} else if (hasc0) {		
		    nc = vp_optimal_scale(nc0, zmin-c0, zmax-c0, &zi, &dc);
		} else {
		    nc = vp_optimal_scale(nc0, zmin, zmax, &c0, &dc);
		}
	    }
	    for (ic=0; ic < nc; ic++) {
		c[ic] = c0 + dc*ic;
	    }
	}

	vp_stdplot_init (min1, max1, min2, max2,
			 transp,false,false,false);
	vp_frame_init(in,"tlb");

	for (ic = 0; ic < nc; ic++) {
	    vp_plot_set (ic);
	    contour (z,n1,n2,c[ic]);
	} 

    } /* i3 */

    exit(0);
}

/*
            north (0)
  (ix,iy+1) --------- (ix+1,iy+1)
            | cell  |
   west (3) | ix,iy | east (1)
            |       |
  (ix,iy)   --------- (ix+1,iy)
            south (2)
*/
static void contour (float **z, int n1, int n2, float c)
{

}
#ifdef IUYFIJF

register int    ix, iy, non;
int             jx, jy, sset (), wset ();
register float  zxymc, zemc, znmc;
float           x, y,  *pzxy;

    /*
     * find all the intersections 
     */
    non = 0;      /* clear intersection counter */
    for (iy = 0; iy < n2 - 1; iy++)
    {
  for (ix = 0; ix < n1 - 1; ix++)
  {

    /******* mask *****
    if (ix == iy || ix + 1 == iy || ix == iy + 1) {
      clrw (pzxy);
      clrs (pzxy);
      continue;
      } */

      pzxy = &z[iy][ix];
      zxymc = (*pzxy) - c;/* z(x,y) - c */
      zemc = z[iy][ix + 1] - c;  /* (z to the east) - c */
      znmc = z[iy + 1][ix] - c;  /* (z to the north) - c */
#define OPPSIGN(A,B) (1==(((A)<0.)+((B)<0.)))
      if (OPPSIGN (zxymc, znmc))  /* if west edge intersected */
      {
    setw (pzxy);  /* set the west bit */
    non++;    /* and increment counter */
      }
      else    /* else */
    clrw (pzxy);  /* clear the west bit */
      if (OPPSIGN (zxymc, zemc))  /* if south edge intersected */
      {
    sets (pzxy);  /* set the south bit */
    non++;    /* and increment counter */
      }
      else    /* else */
    clrs (pzxy);  /* clear the south bit */
  }
    }
    for (ix = 0, iy = n2 - 1; ix < n1 - 1; ix++)  /* northern boundary */
    {

      
    /******* mask *****
      if (ix == iy || ix + 1 == iy) {
	clrs (pzxy);
	clrw (pzxy);
	continue;
	} */

  pzxy = &z[iy][ix];
  zxymc = (*pzxy) - c;  /* z(x,y) - c */
  zemc = z[iy][ix + 1] - c;  /* (z to the east) - c */
  if (OPPSIGN (zxymc, zemc))  /* if south edge intersected */
  {
      sets (pzxy);  /* set the south bit */
      non++;    /* and increment counter */
  }
  else      /* else */
      clrs (pzxy);  /* clear the south bit */
  clrw (pzxy);    /* clear the west bit */
    }
    for (iy = 0, ix = n1 - 1; iy < n2 - 1; iy++)  /* eastern boundary */
    {

      /******* mask *****
    if (ix == iy || ix == iy + 1) {
      clrw (pzxy);
      clrs (pzxy);
      continue;
      } */

  pzxy = &z[iy][ix];
  zxymc = (*pzxy) - c;  /* z(x,y) - c */
  znmc = z[iy + 1][ix] - c;  /* (z to the north) - c */
  if (OPPSIGN (zxymc, znmc))  /* if west edge intersected */
  {
      setw (pzxy);  /* set the west bit */
      non++;    /* and increment counter */
  }
  else      /* else */
      clrw (pzxy);  /* clear the west bit */
  clrs (pzxy);    /* clear the south bit */
    }

    /*
     * draw contours intersecting a boundary 
     */
    for (ix = 0, iy = n2 - 1; ix < n1 - 1 && non > 0; ix++)
  /* north boundary */
    {
  if (sset (&z[iy][ix]))
  {
      x = ix + delta (c, z[iy][ix], z[iy][ix + 1]);
      y = iy;
      moveu (x, y, coordinatec);
      clrs (&z[iy][ix]);
      non--;
      jx = ix;
      jy = iy - 1;
      while (cconnect (z, n1, n2, c, &jx, &jy))
    non--;
  }
    }
    for (ix = n1 - 1, iy = 0; iy < n2 - 1 && non > 0; iy++)  /* east boundary */
    {
  if (wset (&z[iy][ix]))
  {
      x = ix;
      y = iy + delta (c, z[iy][ix], z[iy + 1][ix]);
      moveu (x, y, coordinatec);
      clrw (&z[iy][ix]);
      non--;
      jx = ix - 1;
      jy = iy;
      while (cconnect (z, n1, n2, c, &jx, &jy))
    non--;
  }
    }
    for (ix = 0, iy = 0; ix < n1 - 1 && non > 0; ix++)  /* south boundary */
    {
  if (sset (&z[iy][ix]))
  {
      x = ix + delta (c, z[iy][ix], z[iy][ix + 1]);
      y = iy;
      moveu (x, y, coordinatec);
      clrs (&z[iy][ix]);
      non--;
      jx = ix;
      jy = iy;
      while (cconnect (z, n1, n2, c, &jx, &jy))
    non--;
  }
    }
    for (ix = 0, iy = 0; iy < n2 - 1 && non > 0; iy++)  /* west boundary */
    {
  if (wset (&z[iy][ix]))
  {
      x = ix;
      y = iy + delta (c, z[iy][ix], z[iy + 1][ix]);
      moveu (x, y, coordinatec);
      clrw (&z[iy][ix]);
      non--;
      jx = ix;
      jy = iy;
      while (cconnect (z, n1, n2, c, &jx, &jy))
    non--;
  }
    }

    /*
     * draw interior contours 
     */
    for (iy = 0; iy < n2 - 1 && non > 0; iy++)
    {
  for (ix = 0; ix < n1 - 1 && non > 0; ix++)
  {
      if (sset (&z[iy][ix]))  /* check south edge of cell */
      {
    x = ix + delta (c, z[iy][ix], z[iy][ix + 1]);
    y = iy;
    moveu (x, y, coordinatec);
    clrs (&z[iy][ix]);
    non--;    /* clear start */
    jx = ix;
    jy = iy;
    if (cconnect (z, n1, n2, c, &jx, &jy))
        sets (&z[iy][ix]);  /* finish = start */
    while (cconnect (z, n1, n2, c, &jx, &jy))
        non--;
      }
  }
    }
  return(0);
}

/* 
cconnect draws a line from one intersection of the cell (ix,iy)
   to another intersection of the cell, provided the latter intersection exists,
   and then clears the latter intersection and updates ix and iy.
   cconnect returns 0 if the latter intersection does not exist or if the 
   latter intersection is a grid boundary; otherwise returns 1.
*/

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int cconnect (float **z, int n1, int n2,float  c, int *ix, int *iy)
_XFUNCPROTOEND
#else
int cconnect (z, n1, n2, c, ix, iy)
    int             n1, n2, *ix, *iy;
    float         **z, c;
#endif
{
  int mask;

register int    jx, jy;
float           x, y;

    jx = (*ix);
    jy = (*iy);

    mask = (z[jy][jx] >= 0. && z[jy][jx+1] >= 0. && 
	    z[jy+1][jx] >= 0. && z[jy+1][jx+1] >= 0.);

    if (sset (&z[jy + 1][jx]))  /* if exiting north */
    {
  jy++;

  x = jx + delta (c, z[jy][jx], z[jy][jx + 1]);
  y = jy;

  if (mask) 
    drawu (x, y, coordinatec);
  else
    moveu (x, y, coordinatec);

  clrs (&z[jy][jx]);
  if (++(*iy) >= n2 - 1)
      return (0);
    }
    else
    if (wset (&z[jy][jx + 1]))  /* if exiting east */
    {
  jx++;

  x = jx;
  y = jy + delta (c, z[jy][jx], z[jy + 1][jx]);

  if (mask)  
    drawu (x, y, coordinatec);
  else
    moveu (x, y, coordinatec);

  clrw (&z[jy][jx]);
  if (++(*ix) >= n1 - 1)
      return (0);
    }
    else
    if (sset (&z[jy][jx]))  /* if exiting south */
    {

  x = jx + delta (c, z[jy][jx], z[jy][jx + 1]);
  y = jy;

  if (mask) 
    drawu (x, y, coordinatec);
  else
    moveu (x, y, coordinatec);

  clrs (&z[jy][jx]);
  if (--(*iy) < 0)
      return (0);
    }
    else
    if (wset (&z[jy][jx]))  /* if exiting west */
    {
      
  x = jx;
  y = jy + delta (c, z[jy][jx], z[jy + 1][jx]);

  if (mask) 
    drawu (x, y, coordinatec);
  else
    moveu (x, y, coordinatec);

  clrw (&z[jy][jx]);
  if (--(*ix) < 0)
      return (0);
    }
    else
  return (0);    /* no exit found */
    return (1);
}

/* 
subroutines to set, clear, and check status of bits */
#define SOUTH 0x00000001
#define WEST 0x00000002

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void sets (register int *i)
_XFUNCPROTOBEGIN
#else
void sets (i)
    register int   *i;
#endif
{
    *i |= SOUTH;
}




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void  clrs (register int *i)
_XFUNCPROTOEND
#else
void  clrs (i)
    register int   *i;
#endif
{
    *i &= ~SOUTH;
}
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int  sset (register int *i)
_XFUNCPROTOEND
#else
int  sset (i)
    register int   *i;
#endif
{
    return ((*i) & SOUTH);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void  setw (register int *i)
_XFUNCPROTOEND
#else
void  setw (i)
    register int   *i;
#endif
{
    *i |= WEST;
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void  clrw (register int *i)
_XFUNCPROTOEND
#else
void  clrw (i)
    register int   *i;
#endif
{
    *i &= ~WEST;
}
#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int  wset (register int *i)
_XFUNCPROTOEND
#else
int  wset (i)
    register int   *i;
#endif
{
    return ((*i) & WEST);
}

/* 
subroutine to compute (a-b)/(c-b) for use in linear interpolation */

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
float delta (float a, float b, float  c)
_XFUNCPROTOEND
#else
float delta (a, b, c)
    float           a, b, c;
#endif
{
float           t;

    t = c - b;      /* avoids pathological comparison */
    if (t != 0.0)
  return ((a - b) / t);
    else
  return ((float) 0.5);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int transp (float *x, float *y)
_XFUNCPROTOBEGIN
#else
int transp (x, y)
    float          *x, *y;
#endif

{
float           xyexch;

    xyexch = *x;
    *x = *y;
    *y = xyexch;
  return(0);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int rotate (float x, float *x1, float max, float min)
_XFUNCPROTOEND
#else
int rotate (x, x1, max, min)
    float           x, *x1, max, min;
#endif
{
float           temp, temp2;

    temp = x;
    temp2 = (min + max) - temp;
    *x1 = temp2;
  return(0);
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int draw (float x0, float y0, struct coordinfo coord)
_XFUNCPROTOBEGIN
#else
int draw (x0, y0, coord)
    float           x0, y0;
    struct coordinfo coord;
#endif
{
float           x1, y1;
float x, y; x=x0;y=y0;

    if (coord.transp)
  transp (&x, &y);
    if (coord.xreverse)
    {
  rotate (x, &x2, coord.min1, coord.max1);
  x = x2;
    }
    if (coord.yreverse)
    {
  rotate (y, &x2, coord.min2, coord.max2);
  y = x2;
    }
    gl_draw (x, y);
  return(0);
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int drawu (float x0, float y0, struct coordinfo coord)
_XFUNCPROTOBEGIN
#else
int drawu (x0, y0, coord)
    float           x0, y0;
    struct coordinfo coord;
#endif

{
float           x2;
float x, y; x=x0;y=y0;

    if (coord.transp)
  transp (&x, &y);
    if (coord.xreverse)
    {
  rotate (x, &x2, coord.min1, coord.max1);
  x = x2;
    }
    if (coord.yreverse)
    {
  rotate (y, &x2, coord.min2, coord.max2);
  y = x2;
    }
    gl_udraw (x, y);
  return(0);
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int moveu (float x0, float y0, struct coordinfo coord)
#else
int moveu (x0, y0, coord)
    float           x0, y0;
    struct coordinfo coord;
#endif

{
float           x2;
float x, y; x=x0;y=y0;

    if (coord.transp)
  transp (&x, &y);
    if (coord.xreverse)
    {
  rotate (x, &x2, coord.min1, coord.max1);
  x = x2;
    }
    if (coord.yreverse)
    {
  rotate (y, &x2, coord.min2, coord.max2);
  y = x2;
    }
    gl_umove (x, y);
  return(0);
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int  move (float x0, float y0,  struct coordinfo coord)
_XFUNCPROTOEND
#else
int  move (x0, y0, coord)
    float           x0, y0;
    struct coordinfo coord;
#endif

{
float           x2;
float x, y; x=x0;y=y0;

    if (coord.transp)
  transp (&x, &y);
    if (coord.xreverse)
    {
  rotate (x, &x2, coord.min1, coord.max1);
  x = x2;
    }
    if (coord.yreverse)
    {
  rotate (y, &x2, coord.min2, coord.max2);
  y = x2;
    }
    gl_umove (x, y);
  return(0);
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int where (float *x,float *y,  struct coordinfo coord)
_XFUNCPROTOEND
#else
int where (x, y, coord)
    float          *x, *y;
    struct coordinfo coord;
#endif

{
float           x1, y1;

    gl_where (&x1, &y1);
    if (coord.transp)
  transp (&x1, &y1);
    if (coord.xreverse)
    {
  rotate (x1, &x2, coord.min1, coord.max1);
  x1 = x2;
    }
    if (coord.yreverse)
    {
  rotate (y1, &x2, coord.min2, coord.max2);
  y1 = x2;
    }
    *x = x1;
    *y = y1;
  return(0);
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int markup (int npts,int  mtype, int  msize, float  x0, float  y0, struct coordinfo coord)
_XFUNCPROTOEND
#else
int markup (npts, mtype, msize, x0, y0, coord)
    int             npts, mtype, msize;
    float           x0, y0;
    struct coordinfo coord;
#endif

{
float x,y; x=x0;y=y0;
    if (coord.transp)
  transp (&x, &y);
    if (coord.xreverse)
    {
  rotate (x, &x2, coord.min1, coord.max1);
  x = x2;
    }
    if (coord.yreverse)
    {
  rotate (y, &x2, coord.min2, coord.max2);
  y = x2;
    }
    gl_upmark (npts, mtype, msize, x, y);
  return(0);
}

#endif
