/* Regular grid routines. */
/*************************************************************************

Copyright Rice University, 2008.
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

#include <trip/base.h>

#include "thegrid.h"

#ifndef _sf_thegrid_h

#define TOL 0.01
/*^*/

/** Axis: Basic Metadata element for regular grids. Defines a uniformly
    sampled axis, after the fashion of SEPlib77 or RSF. Struct part of
    RVLGrid Axis class. 
    *
    Innovative feature: each axis carries an integer token (axis.id)
    meant to signify its position in a global ordering of axes (hence
    its physical meaning)..
*/
typedef struct {
  /** number of gridpoints on axis */
  size_t n; 
  /** step between gridpoints */
  ireal d;
  /** coordinate of first gridpoint */
  ireal o;
  /** axis index */
  int id;
} axis;
/*^*/

/** Regular grid struct */
typedef struct {
  size_t dim;
  axis axes[RARR_MAX_NDIM];
} grid;
/*^*/

#endif

int init_default_axis(axis * a) 
/*< ** default constructor. Default values: n=1, d=1.0, o=0.0, id=1
 @param[out]  a (axis *) - axis to be initialized 
>*/
{
  a->n=0;
  a->d=1.0;
  a->o=0.0;
  a->id=0;
  return 0;
}

int init_axis(axis * a, size_t n, ireal d, ireal o) 
/*< ** main constructor - not terribly useful. id defaults to 1.
@param[out] a (axis *) - axis to be initialized 
@param[in] n (int) - number of gridpoints on axis
@param[in] d (ireal) - step
@param[in] o (ireal) - coord of first data point
>*/
{
  a->n=n;
  a->d=d;
  a->o=o;
  a->id=0;
  return 0;
}

int fprint_axis(FILE * fp, axis a) 
/*< ** print axis to stream 
@param[in] fp (FILE *) - output stream
f@param[in] a (axis) - axis to be printed
 >*/
{
    fprintf(fp,"axis: n=%d d=%e o=%e id=%d\n",(int) a.n,a.d,a.o,a.id);
  return 0;
}

int print_axis(axis a) 
/*< ** print axis to stdout
@param[in] a (axis) - axis to be printed
>*/
{
  return fprint_axis(stdout,a);
}

int compare_axis(axis a1, axis a2) 
/*< ** compare two axes for equality
@param[in] a1 (axis) - first axis 
@param[in] a2 (axis) - second axis
@return 0 if axes are same, else 1
>*/
{
  int err=0;
  err = err || (a1.n != a2.n);
  err = err || (fabs((double)(a1.d-a2.d)) > TOL*fmin((double)(a1.d),(double)(a2.d)));
  err = err || (fabs((double)(a1.o-a2.o)) > TOL*fmin((double)(a1.d),(double)(a2.d)));
  err = err || (a1.id != a2.id);
  return err;
}

int init_default_grid(grid * g) 
/*< ** default constructor. 
@param[out] g (grid *) - grid to be initialized
>*/
{
  int i;
  g->dim=0;
  for (i=0;i<RARR_MAX_NDIM;i++) {
    init_default_axis(&(g->axes[i]));
    g->axes[i].id=i;
  }
  return 0;
}

int init_grid(grid * g, size_t dim) 
/*< ** main constructor. initializes metadata (dim) but not data
    (axes). Since no dynamic allocation takes place, only real role is
    to sanity-check dim.
@param[out] g (grid *) - grid to be initialized
@param[in] dim (int) - dimension of grid, at most \ref RARR_MAX_NDIM
>*/
{
  int i;
  if (dim<1 || dim>RARR_MAX_NDIM) return E_BADINPUT;
  g->dim=dim;
  for (i=0;i<RARR_MAX_NDIM;i++) { 
    init_default_axis(&(g->axes[i]));
    g->axes[i].id=i;
  }
  return 0;
}

int fprint_grid(FILE * fp, grid a) 
/*< ** print grid to stream 
@param[in] a (grid) - grid to be printed
@param[in] fp (FILE *) - stream to which to print
>*/
{
  int i;
  fprintf(fp,"Grid data structure, consisting of %d axes:\n",(int) a.dim);
  for (i=0;i<a.dim;i++) fprint_axis(fp,a.axes[i]);
  return 0;
}

int print_grid(grid a) 
/*< * print grid to stdout 
@param[in] a (grid) - grid to be printed
>*/
{ return fprint_grid(stdout,a); }

int compare_grid(grid g1, grid g2) 
/*< * compare two grids by comparing their axes 
@param[in] g1 (grid) - first grid
@param[in] g2 (grid) - second grid
@return 0 if grids same, else 1
>*/
{
  int err=0;
  int i;
  if (g1.dim != g2.dim) return 1;
  for (i=0;i<g1.dim;i++) err = err || compare_axis(g1.axes[i],g2.axes[i]);
  return err;
}

int get_datasize_grid(grid g) 
/*< ** get number of gridpoints (product of n's)
@param[in] g (grid) - input grid
@return product of axis lengths
>*/
{
  _IPNT _n;
  int i;
  int sz=1;
  get_n(_n,g);
  for (i=0;i<g.dim;i++) sz*=_n[i];
  return sz;
}

int get_n(_IPNT n, grid g) 
/*< ** get axis length array 
@param[in] g (grid) - input grid
@param[out] n (IPNT) - axis lengths
>*/
{
  int i; 
  for (i=0;i<RARR_MAX_NDIM;i++) {
    n[i]=g.axes[i].n;
  }
  return 0;
}

int get_d(_RPNT d, grid g) 
/*< ** get grid origin coordinate array 
@param[in] g (grid) - input grid
@param[out] o (RPNT) - grid origin
>*/
{
  int i; 
  for (i=0;i<RARR_MAX_NDIM;i++) {
    d[i]=g.axes[i].d;
  }
  return 0;
}

int get_o(_RPNT o, grid g) 
/*< * get array of indices of grid origin in global grid 
@param[in] g (grid) - input grid
@param[out] gs (IPNT) - global indices of grid origin
>*/
{
  int i; 
  for (i=0;i<RARR_MAX_NDIM;i++) {
    o[i]=g.axes[i].o;
  }
  return 0;
}

int get_gs(_IPNT gs, grid g) 
/*< * get array of indices of grid origin in global grid 
@param[in] g (grid) - input grid
@param[out] gs (IPNT) - global indices of grid origin
>*/
{
  int i;
  for (i=0;i<g.dim;i++) {
    if (g.axes[i].o<0) gs[i]=(int)((g.axes[i].o-g.axes[i].d*TOL)/(g.axes[i].d));
    else gs[i]=(int)((g.axes[i].o+g.axes[i].d*TOL)/(g.axes[i].d));
  }
  return 0;
}

int get_ord(_IPNT od, grid g) 
/*< * returns axis order array, i.e. axis index as a function of id,
    rather than id as a function of axis index (which is stored). 
@param[in] g (grid) - input 9grid
@param[out] a (IPNT) - axis order array
>*/ 
{
  int i,j;
  for (i=0;i<RARR_MAX_NDIM;i++) {
    for (j=0;j<RARR_MAX_NDIM;j++) {
      if (g.axes[i].id == j) od[j]=i;
    }
  }
  return 0;
}
