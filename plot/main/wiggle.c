/* Plot data with wiggly traces.

Takes < in.rsf > out.vpl
*/ 

#include <rsf.h>
#include <rsfplot.h>

MAIN ()
{
   vp_filep(outstream); /* tell plotting routines where to stick it */

    initial1 ();
    xpostion ();
    initial2 ();
    memoryallocation ();
    update ();
    counter = 0;
    framecounter = 0;
    graphinitial ();
    for (; n3c > 0; n3c--)
    {
  datainput ();
  findzdata ();
  plotwiggle ();
    }

  return(0);
}

/* initial1 initializes and fetches n1, n2, n3 o1, d1
*/


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void initial1 (void )
_XFUNCPROTOEND
#else
void initial1 ()
#endif
{
int             ii, n1num, d1num, o1num, dummy;

/* initializing data space */
    fastplt = 0;
    data.esize= 4;
    data.n1[0] = 0;
    data.d1[0] = 1;
    data.o1[0] = data.d1[1];



    n1num = fetch ("n1", "d", data.n1);
    switch (n1num)
    {
    case 0:
  seperr ("n1 is needed!\n");
  break;
    case 1:
  break;
    default:
  if (n1num > NPMAX)
      seperr ("entered too many values for n1  exceeded NPMAX\n");
  fprintf (stderr, "************** Warning *******************\n");
  fprintf (stderr, "Too many values for n1 were entered,\n");
  fprintf (stderr, "using only first one (%d) \n ", data.n1[0]);
  fprintf (stderr, "******************************************\n");
  for (ii = 1; ii < NPMAX; ii++)
      data.n1[ii] = data.n1[0];
  break;
    }

    if (fetch ("n3", "d", &data.n3) == 0)
  data.n3 = 1;
    if (data.n3 == 0)
  seperr ("n3 = 0 assuming there is no data\n");
    ii = ssize ("in") / (data.n1[0] * data.n3 * sizeof (float));
/* setting n2 to the number of n1's */

    if (fetch ("n2", "d", &data.n2) == 0)
    {
  data.n2 = ii;
    }
    n3c = data.n3;
    d1num = fetch ("d1", "f", data.d1);
    switch (d1num)
    {
    case 0:
  break;
    case 1:
  break;
    default:
  if (d1num > NPMAX)
      seperr ("entered too many values for d1  exceeded NPMAX\n");
  fprintf (stderr, "************** Warning *******************\n");
  fprintf (stderr, "Too many values for d1 were entered,\n");
  fprintf (stderr, "using only first one (%f) \n ", data.d1[0]);
  fprintf (stderr, "******************************************\n");
  for (ii = 1; ii < NPMAX; ii++)
      data.d1[ii] = data.d1[0];
  break;
    }
    dth = data.d1[0] / 2.;
    o1num = fetch ("o1", "f", data.o1);
    switch (o1num)
    {
    case 0:
  break;
    case 1:
  break;
    default:
  if (o1num > NPMAX)
      seperr ("entered too many values for o1  exceeded NPMAX\n");
  fprintf (stderr, "************** Warning *******************\n");
  fprintf (stderr, "Too many values for o1 were entered,\n");
  fprintf (stderr, "using only first one (%f) \n ", data.o1[0]);
  fprintf (stderr, "******************************************\n");
  for (ii = 1; ii < NPMAX; ii++)
      data.o1[ii] = data.o1[0];
  break;
    }
    dummy = fetch ("o3", "f", &data.o3);
    dummy = fetch ("d3 ", "f", &data.d3);
    dummy = fetch ("esize ", "d", &data.esize);
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void xpostion (void)
_XFUNCPROTOEND
#else
void xpostion ()
#endif
{
int             dummy;

    if (auxin ("xpos") != NULL)
       xposfd = fileno(auxin ("xpos"));
    else
       xposfd = -1;
    xpos = (float *) calloc (data.n2, sizeof (float));
    if (xposfd != -1)
    {
  if (amountread = sreed ("xpos", xpos, data.n2 * data.esize) != data.n2 * data.esize)
  {
      fprintf (stderr, "******************** WARNING ***********************");
      fprintf (stderr, "Amount of data read in xpos was not amount specified\n");
      fprintf (stderr, "check n2 \n");
      fprintf (stderr, "****************************************************");
  }
    }
    if (xposfd == -1)
    {
  if (!fetch ("d2", "f", &data.d2))
      data.d2 = 1;
  if (!fetch ("o2", "f", &data.o2))
      data.o2 = data.d2;
  coordinate.min2 = data.o2;
  coordinate.max2 = data.o2 + (data.n2 - 1) * data.d2;
  ix0 = (coordinate.min2 - data.o2) / data.d2 - .5;
  if (ix0 < 0)
      ix0 = 0;
  ixmax = (coordinate.max2 - data.o2) / data.d2 + .5;
  if (ixmax < 0)
      ixmax = data.n2 - 1;

  for (ix = 0; ix < data.n2; ix++)
  {
      xpos[ix] = data.o2 + ix * data.d2;
  }
    }
    else
    {
  coordinate.min2 = xpos[0];
  coordinate.max2 = xpos[data.n2 - 1];
  if (data.n2 == 1)
      data.d2 = (coordinate.max2 - coordinate.min2);
  else
      data.d2 = (coordinate.max2 - coordinate.min2) / (data.n2 - 1);
  dummy = fetch ("d2", "f", &data.d2);
  for (ix = 0; ix < data.n2; ix++)
  {
      if (xpos[ix] >= coordinate.min2)
    break;
  }
  ix0 = ix - 1;
  if (ix0 < 0)
      ix0 = 0;
  for (ix = data.n2 - 1; ix >= 0; ix--)
  {
      if (xpos[ix] <= coordinate.max2)
    break;
  }
  ixmax = ix + 1;
  if (ixmax >= data.n2)
      ixmax = data.n2 - 1;
    }
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void initial2 (void)
_XFUNCPROTOEND
#else
void initial2 ()
#endif
{
    if (!getch ("overplot", "1", &overplot))
    {
  overplot = 0;
    }
    /* polygon stuff */
    ipoly = 0;
    if (!getch ("poly", "1", &poly))
  poly = 0;
    if (poly)
    {
  ipoly = 1;
  if (!getch ("fatp", "d", &fatp))
      fatp = 1;
  if (!getch ("xmask", "d", &xmask))
      xmask = 1;
  if (!getch ("ymask", "d", &ymask))
      ymask = 1;
    }
    if (!getch ("ntile", "f", &ntile))
  if (!getch ("pclip", "f", &ntile))
      ntile = 98;
    if (!getch ("zplot", "f", &zplot))
  zplot = .75;
    if (!getch ("clip", "f", &zdata))
  zdata = 0;
    zplot *= data.d2;
    coordinate.min1 = data.o1[0];
    coordinate.max1 = data.o1[0] + (data.n1[0] - 1) * data.d1[0];
    it0 = (coordinate.min1 - data.o1[0]) / data.d1[0] - .5;
    if (it0 < 0)
  it0 = 0;
    itmax = (coordinate.max1 - data.o1[0]) / data.d1[0] + .5;
    if (itmax > data.n1[0])
  itmax = data.n1[0] - 1;
    coordinate.min2 = coordinate.min2 - zplot;
    coordinate.max2 = coordinate.max2 + zplot;
    if (!getch ("seemean", "1", &seemean))
  seemean = 0;
    if (!getch ("preder", "1", &npreder))
  npreder = 0;
}

#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void memoryallocation (void)
_XFUNCPROTOEND
#else
void memoryallocation ()
#endif
{
    ptile = (float *) calloc ((data.n1[0] + data.n2), sizeof(float));
    pdata = (float *) calloc ((npreder + 1) * data.n1[0] * data.n2, sizeof(float));
    position1 = (float *) calloc (data.n1[0], sizeof(float));
    position2 = (float *) calloc (data.n1[0], sizeof(float));
    p = (float *) calloc (data.n1[0], sizeof(float));
    ip = (int *) (ptile + data.n2);
    px = (float *) calloc (2 + data.n1[0], sizeof (float));
    py = (float *) calloc (2 + data.n1[0], sizeof (float));
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void datainput (void)
_XFUNCPROTOEND
#else
void datainput ()
#endif
{
int             temp;

    temp = data.n1[0] * data.n2 * data.esize;
    if (amountread = sreed ("in", pdata, temp) != temp)
    {
  fprintf (stderr,"\n******************** WARNING ***********************");
  fprintf (stderr,"\nAmount of data read in is not amount specified\n");
  fprintf (stderr,"\nAmount read=%d, Amount expected=%d.",
     amountread,temp);
  fprintf (stderr,"\nCheck n1, n2 values.");
  fprintf (stderr,"\n****************************************************");
    }

}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void update (void)
_XFUNCPROTOEND
#else
void update ()
#endif
{
    minusone = -1;
    one = 1;
    putch ("n1", "d", &minusone);
    putch ("n2", "d", &one);
    putch ("n3", "d", &one);
    set_output_data_format("vplot");
    hclose ();
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void graphinitial (void)
_XFUNCPROTOEND
#else
void graphinitial ()
#endif
{
float           tempminmax;

    if (!getch ("wantframe", "1", &wantframe))
  wantframe = 1;
    if (!getch ("wantframenum", "1", &wantframenum))
  wantframenum = 1;
    strcpy (axis1.wherelabel, "b");
    strcpy (axis2.wherelabel, "l");
    coordinate.transp = 0;
    coordinate.xreverse = 0;
    coordinate.yreverse = 0;
    strcpy (title.wheretitle, "t");
    gl_coordint (&position, &coordinate, &axis1, &axis2);
/*    if (coordinate.transp)
    {
    tempminmax = coordinate.max1;
    coordinate.max1 = coordinate.max2;
    coordinate.max2 = tempminmax;
    tempminmax = coordinate.min1;
    coordinate.min1 = coordinate.min2;
    coordinate.min2 = tempminmax;
    } 
    getminmax(&coordinate);
*/
    gl_minmax (&coordinate);
    coordinate.pad = 0;
    coordinate.npad = getch ("pad", "1", &coordinate.pad);
    gl_padint (&coordinate);
    gl_axisint (&axis1, &axis2, &coordinate, &position);

    gl_gridint(&grid, &coordinate, &axis1, &axis2);
    gl_titleint (&title);
    gl_colorint (&color);
    plotint ();
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void getminmax (struct coordinfo *coordinate)
_XFUNCPROTOEND
#else
void getminmax (coordinate)
    struct coordinfo *coordinate;
#endif
{
int             dummy;

/*
* This routine will fetch the min,and max values  and pad. 
*/
    dummy = getch ("min1 tmin", "f", &coordinate->min1);
    dummy = getch ("max1 tmax", "f", &coordinate->max1);
    dummy = getch ("min2 xmin", "f", &coordinate->min2);
    dummy = getch ("max2 xmax", "f", &coordinate->max2);
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void plotint (void)
_XFUNCPROTOEND
#else
void plotint ()
#endif

{
int i;
    if (!getch ("tlinecol", "d", timeline.col))
  for( i=0;i<n3c; timeline.col[i++] = 2);
    if (!getch ("tlinefat", "d", &timeline.fat))
  timeline.fat = 0;
    if (!getch ("plotcol tracecol", "d", plot.col))
  for( i=0;i<n3c; plot.col[i++] = 6);
    if (!getch ("plotfat tracefat", "d", plot.fat))
  for( i=0;i<n3c; plot.fat[i++] = 0);
    if (coordinate.transp)
    {
  timeline.grid1 = 0;
  if (grid.grid2) 
           timeline.grid2 = 0;
        else
           timeline.grid2 = 1;
  timeline.g2num = axis2.dnum;
  timeline.g1num = axis2.dnum;
        
    }
    else
    {
  timeline.grid2 = 0;
  if(grid.grid1)
  timeline.grid1 = 0;
  else
  timeline.grid1 = 1;
  timeline.g1num = axis1.dnum;
    }
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void findzdata (void)
_XFUNCPROTOEND
#else
void findzdata ()
#endif
{
    itile = it0 + (itmax - it0) * ntile / 100.;
    if (zdata <= 0.)
    {
  for (ix = ix0; ix <= ixmax; ix++)
  {
      q = pdata + data.n1[0] * ix;
      qe = q + data.n1[0] * data.n2;
      for (it = it0; it <= itmax; it++)
      {
    ip[it] = it;
    p[it] = fabs (q[it]);
      }
      pfind (p, ip, it0, itmax, itile);
      ptile[ix] = p[ip[itile]];
  }
  for (ix = ix0; ix <= ixmax; ix++)
      ip[ix] = ix;
  itile = ix0 + (ixmax - ix0) * ntile / 100;
  pfind (ptile, ip, ix0, ixmax, itile);
  zdata = ptile[ip[itile]];
    }
    if (zdata <= 0.)
  /*ntile percent of data are zero.  Define: zdata = L1 (norm) = sum (
          abs(non-zeros) ) / # of non-zeros */
    {
  scale = 0;
  counts = 0;
  for (ix = 0; ix < data.n2; ix++)
  {
      for (it = 0; it < data.n1[0]; it++)
      {
    if (pdata[ix * data.n1[0] + it] != 0.)
    {
        scale += fabs (pdata[ix * data.n1[0] + it]);
        counts += 1;
    }
      }
  }
  if (scale > 0.)
      scale = counts * zplot / scale;

    }
    else
  scale = zplot / zdata;
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void pfind (float *x,int *p, int m, int n, int k)
_XFUNCPROTOEND
#else
void pfind (x, p, m, n, k)
    int            *p, m, n, k;
    float          *x;
#endif
{
int             i, j;

    if (m < n)
    {
  ppart (x, p, m, n, &i, &j);
  if (k <= j)
      pfind (x, p, m, j, k);
  else
  if (i <= k)
      pfind (x, p, i, n, k);
    }
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void ppart (float *x,int *p,int m,int  n, int *pk,int *pj)
_XFUNCPROTOEND
#else
void ppart (x, p, m, n, pk, pj)
    int            *p, m, n, *pk, *pj;
    float          *x;
#endif
{
int             i, j, f;
float           xx;

    f = myrandom (m, n);
    xx = x[p[f]];
    i = m;
    j = n;
up:for (i = i; i <= n; i++)
  if (xx < x[p[i]])
      goto down;
    i = n;
down:for (j = j; j >= m; j--)
  if (x[p[j]] < xx)
      goto change;
    j = m;
change:if (i < j)
    {
  exchge (&p[i], &p[j]);
  i++;
  j--;
  goto up;
    }
    else
    if (i, f)
    {
  exchge (&p[i], &p[f]);
  i++;
    }
    else
    if (f < j)
    {
  exchge (&p[f], &p[j]);
  j--;
    }
    *pk = i;
    *pj = j;
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
int myrandom (int m, int n)
_XFUNCPROTOEND
#else
int myrandom (m, n)
    int             m, n;
#endif
{
float           r;

    r = 0.0;
    r=franuni();
    r = m + (n - m) * r;
    return ((int) r);
}


/* Copyright (c) Colorado School of Mines, 1999.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
FRANUNI - Functions to generate a pseudo-random float uniformly distributed
  on [0,1); i.e., between 0.0 (inclusive) and 1.0 (exclusive)

franuni   return a random float

******************************************************************************
Function Prototypes:
float franuni (void);
void sranuni (int seed);

******************************************************************************
franuni:
Input:    (none)
Returned: pseudo-random float

sranuni:
seed    different seeds yield different sequences of random numbers.

******************************************************************************
Notes:
Adapted from subroutine uni in Kahaner, Moler, and Nash (1988).
This book references a set of unpublished notes by
Marsaglia.

According to the reference, this random
number generator "passes all known tests and has a period that is ...
approximately 10^19".

******************************************************************************
******************************************************************************
References:
"Numerical Methods and Software", D. Kahaner, C. Moler, S. Nash,
Prentice Hall, 1988.

Marsaglia G., "Comments on the perfect uniform random number generator",
Unpublished notes, Wash S. U.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 12/30/89
*****************************************************************************/


/* constants used to generate random numbers (16777216=2^24) */
#define CS 362436.0/16777216.0
#define CD 7654321.0/16777216.0
#define CM 16777213.0/16777216.0
#define NBITS 24

/* internal state variables */
static int i=16,j=4;
static float c=CS;
static float u[]={
  0.8668672834288,  0.3697986366357,  0.8008968294805,
  0.4173889774680,  0.8254561579836,  0.9640965269077,
  0.4508667414265,  0.6451309529668,  0.1645456024730,
  0.2787901807898,  0.06761531340295, 0.9663226330820,
  0.01963343943798, 0.02947398211399, 0.1636231515294,
  0.3976343250467,  0.2631008574685
};





#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
float franuni (void)
_XFUNCPROTOEND
#else
float franuni ()
#endif

{ float uni;

  /* basic generator is Fibonacci */
  uni = u[i]-u[j];
  if (uni<0.0) uni += 1.0;
  u[i] = uni;
  i--;
  if (i<0) i = 16;
  j--;
  if (j<0) j = 16;

  /* second generator is congruential */
  c -= CD;
  if (c<0.0) c += CM;

  /* combination generator */
  uni -= c;
  if (uni<0.0) uni += 1.0;
  return uni;
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void exchge (int *pa,int  *pb)
_XFUNCPROTOEND
#else
void exchge (pa, pb)
    int            *pa, *pb;
#endif
{
int             temp;

    temp = *pa;
    *pa = *pb;
    *pb = temp;
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void plotwiggle (void)
_XFUNCPROTOEND
#else
void plotwiggle ()
#endif
{
float           temp1, temp2;

    if (npreder && n3c > 1)
    {
  sreed ("in", pdata + data.n1[0] * data.n2, data.n1[0] * data.n2 * data.esize);
  n3c--;
  ipreder = 1;
    }
    else
  ipreder = 0;
    if (overplot == 0 || n3c == data.n3)
    {
  gl_erase ();

    }

    gl_vplotint (&position, &coordinate, &axis1, &axis2);
    gl_plotpram (&color, &coordinate);
    gl_color (plot.col[0]);
    gl_fat (plot.fat[0]);
    if (overplot)
    {
  plot.col[0]--;
/* DAMN IT IF THEY ASK FOR COLOR ZERO THEN THEY WANT COLOR ZERO, not YELLOW
 *  if (plot.col[0] == 0)
 *      plot.col[0] = 6;
 */
    }
    /* draw n2 traces */
    for (ix = ix0; ix <= ixmax; ix++)
    {
  q = pdata + data.n1[0] * ix;
  qe = q + data.n1[0] * data.n2;
  if (coordinate.transp)
  {
      for (it = it0; it <= itmax; it++)
      {
    position2[it] = data.o1[0] + it * data.d1[0];
    position1[it] = xpos[ix] + q[it] * scale;
      }
  }
  else
  {
      for (it = it0; it <= itmax; it++)
      {
    position1[it] = data.o1[0] + it * data.d1[0];
    position2[it] = xpos[ix] + q[it] * scale;
      }
  } 
      for (it = it0; it <= itmax; it++)
      {
          check1 (&position1[it], &position2[it], coordinate);
      }

  for (it = it0; it <= itmax; it++)
  {
      if (it == it0)
      {
    gl_umove (position1[it], position2[it]);
      }
      else
      {
    gl_udraw (position1[it], position2[it]);
      }
  }
  if (seemean)    /* plot mean lines of traces -jon */
  {
      it = it0;
      cposition1 = data.o1[0] + it * data.d1[0];
      cposition2 = xpos[ix];
      check (&cposition1, &cposition2, coordinate);
      gl_umove (cposition1, cposition2);
      it = itmax;
      cposition1 = data.o1[0] + it * data.d1[0];
      cposition2 = xpos[ix];
      check (&cposition1, &cposition2, coordinate);
      gl_udraw (cposition1, cposition2);
  }
  /* draw nx traces with polygon fillings */
  if (ipoly == 1)
  {
      iii = 0;
      if (q[it0] > 0.)
      {
    cposition1 = data.o1[0] + it0 * data.d1[0];
    cposition2 = xpos[ix];
    check (&cposition1, &cposition2, coordinate);
    px[iii] = cposition1;
    py[iii] = cposition2;
    iii += 1;
    px[iii] = position1[it0];
    py[iii] = position2[it0];
    iii += 1;
      }

      for (it = it0 + 1; it <= itmax - 1; it++)
      {

    /*
     * see if interpolation needed 
     */
    if (q[it] > 0 && q[it - 1] <= 0.)
    {
        cposition1 = data.o1[0] + it * data.d1[0];
        px[iii] = cposition1 - data.d1[0] * q[it] / (q[it] - q[it - 1]);
        py[iii] = xpos[ix];;
        check (&px[iii], &py[iii], coordinate);
        iii += 1;
        px[iii] = position1[it];
        py[iii] = position2[it];
        iii += 1;
    }
    else
    if (q[it] > 0. && iii != 0)
    {
        px[iii] = position1[it];
        py[iii] = position2[it];
        iii += 1;
    }
    /* polygon fillings of positive peak */
    else
    if (q[it] <= 0. && q[it - 1] > 0.)
    {
        cposition1 = data.o1[0] + it * data.d1[0];
        px[iii] = cposition1 - data.d1[0] * q[it] / (q[it] - q[it - 1]);
        py[iii] = xpos[ix];;
        check (&px[iii], &py[iii], coordinate);
        iii = iii + 1;
        gl_uarea (px, py, iii, fatp, ymask, xmask, 1);
        iii = 0;
    }
      }
      /* check last sample of trace */
      if (q[itmax] <= 0. && iii > 0.)
      {
    cposition1 = data.o1[0] + itmax * data.d1[0];
    px[iii] = cposition1 - data.d1[0] * q[itmax] / (q[itmax] - q[itmax - 1]);
    py[iii] = xpos[ix];
    check (&px[iii], &py[iii], coordinate);
    iii = iii + 1;
    gl_uarea (px, py, iii, fatp, ymask, xmask, 2);
    iii = 0;
      }
      else
      if (q[itmax] > 0. && iii != 0)
      {
    px[iii] = position1[itmax];
    py[iii] = position2[itmax];
    iii = iii + 1;
    cposition1 = data.o1[0] + itmax * data.d1[0];
    px[iii] = cposition1;
    py[iii] = xpos[ix];
    check (&px[iii], &py[iii], coordinate);
    iii = iii + 1;
    gl_uarea (px, py, iii, fatp, ymask, xmask, 3);
    iii = 0;
      }
      else
      if (q[itmax] > 0. && iii == 0)
      {
    cposition1 = data.o1[0] + itmax * data.d1[0];
    px[iii] = cposition1 - data.d1[0] * q[itmax] / (q[itmax] - q[itmax - 1]);
    py[iii] = xpos[ix];
    check (&px[iii], &py[iii], coordinate);
    iii = iii + 1;
    px[iii] = position1[itmax];
    py[iii] = position2[itmax];
    iii = iii + 1;
    cposition1 = data.o1[0] + itmax * data.d1[0];
    px[iii] = cposition1;
    py[iii] = xpos[ix];
    check (&px[iii], &py[iii], coordinate);
    iii = iii + 1;
    gl_uarea (px, py, iii, fatp, ymask, xmask, 4);
    iii = 0;
      }
  }
  if (ipreder)
  {
      for (it = it0; it <= itmax; it++)
      {
    value = q[it] * qe[it] - qe[it] * qe[it];
    cposition1 = data.o1[0] + it * data.d1[0];
    if (value > 0.0)
    {
        cposition2 = xpos[ix] + q[it] * scale;
        gl_umove (cposition1, cposition2);
        cposition2 = xpos[ix] + qe[it] * scale;
        gl_udraw (cposition1, cposition2);
    }
    else
    if (value < 0.0)
    {
        cposition1 = data.o1[0] + it * data.d1[0];
        cposition2 = xpos[ix] + qe[it] * scale;
        gl_umove (cposition1, cposition2);
        cposition1 = dth + data.o1[0] + it * data.d1[0];
        cposition2 = xpos[ix] + qe[it] * scale;
        gl_udraw (cposition1, cposition2);

    }
      }
  }
    }
    if (overplot == 0 || n3c == 1)
    {
  if (fastplt < 20)
  {
      /*plot timing lines */
      if ( timeline.col[counter] != 0 && fastplt < 10)
      {
              
    if (timeline.grid1)
                   if( timeline.g1num != 0.0)
        gl_plotgrid (&coordinate, &axis1, &timeline, counter);
    if (timeline.grid2)
                   if( timeline.g2num != 0.0)
        gl_plotgrid (&coordinate, &axis2, &timeline, counter);
      }
}
  gl_stdplot (&data, &coordinate, &axis1, &axis2, &grid, &title, counter, fastplt, wantframe,wantframenum);
  counter++;
    }

}




#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void rotate (float x,float *x2,float  min,float  max,struct datainfo data)
_XFUNCPROTOBEGIN
#else
void rotate (x, x2, min, max, data)
    float          x, *x2, min, max;
    struct datainfo data;
#endif
{
int             i, j;
float temp, temp1;
       temp = x;
       temp1 = (min + max) - x;
       *x2 = temp1;
}


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void transp (float *x, float *y, struct datainfo data)
_XFUNCPROTOEND
#else
void transp (x,y, data)
    float          *x, *y;
    struct datainfo data;
#endif
{
float          *xyexch;
int             i, tempalloc, tempval;

    tempalloc = data.n1[0];
    xyexch = (float *) calloc ((data.esize / 2) * (tempalloc), sizeof (float));

    for (i = 0; i < tempalloc; i++)
    {
  xyexch[i] = x[i];
  x[i] = y[i];
  y[i] = xyexch[i];
    }

}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void check1 (float *x,float *y,struct coordinfo coord)
_XFUNCPROTOEND
#else
void check1 (x, y, coord)
    float          *x, *y;
    struct coordinfo coord;
#endif
{
float           xyexch;

    if (coord.xreverse)
    {
  *x = coord.min1 + coord.max1 - *x;
    }
    if (coord.yreverse)
    {
  *y = coord.min2 + coord.max2 - *y;
    }
}



#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
void check (float *x, float *y, struct coordinfo coord)
_XFUNCPROTOEND
#else
void check (x, y, coord)
    float          *x, *y;
    struct coordinfo coord;
#endif
{
float           xyexch;

    if (coord.transp)
    {
  xyexch = *x;
  *x = *y;
  *y = xyexch;
    }
    if (coord.xreverse)
    {
  *x = coord.min1 + coord.max1 - *x;
    }
    if (coord.yreverse)
    {
  *y = coord.min2 + coord.max2 - *y;
    }
}
