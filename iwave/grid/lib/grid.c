#include "grid.h"

void init_default_axis(axis * a) {
  a->n=0;
  a->d=1.0;
  a->o=0.0;
  a->id=0;
}

void init_axis(axis * a, size_t n, ireal d, ireal o) {
  a->n=n;
  a->d=d;
  a->o=o;
  a->id=0;
}

void copy_axis(axis * tgt, const axis * src) {
  tgt->n = src->n;
  tgt->d = src->d;  
  tgt->o = src->o;
  tgt->id = src->id;
}

void fprint_axis(FILE * fp, axis a) {
  fprintf(fp,"axis: n=%ld d=%e o=%e id=%d\n",(long) a.n,a.d,a.o,a.id);
}

void print_axis(axis a) {
  fprint_axis(stdout,a);
}

#define MIN(a,b) ((a) < (b) ? (a) : (b))

int compare_axis(axis a1, axis a2) {
  int err=0;
  err = err || (a1.n != a2.n);
  err = err || (fabs((double)(a1.d-a2.d)) > TOL*MIN((double)(a1.d),(double)(a2.d)));
  err = err || (fabs((double)(a1.o-a2.o)) > TOL*MIN((double)(a1.d),(double)(a2.d)));
  err = err || (a1.id != a2.id);
  return err;
}

void init_default_grid(grid * g) {
  int i;
  g->dim=0;
  g->gdim=0;
  for (i=0;i<RARR_MAX_NDIM;i++) {
    init_default_axis(&(g->axes[i]));
    g->axes[i].id=i;
  }
}

int init_grid(grid * g, int dim, int gdim) {
  int i;
  if (dim<1 || dim>gdim || gdim>RARR_MAX_NDIM) return E_BADINPUT;
  g->dim=dim;
  g->gdim=gdim;
  for (i=0;i<RARR_MAX_NDIM;i++) { 
    init_default_axis(&(g->axes[i]));
    g->axes[i].id=i;
  }
  return 0;
}

void copy_grid(grid * tgt, const grid * src) {
  int i;
  tgt->dim  = src->dim;
  tgt->gdim = src->gdim;
  for (i=0;i<RARR_MAX_NDIM;i++) 
    copy_axis(&(tgt->axes[i]),&(src->axes[i]));
}

void fprint_grid(FILE * fp, grid a) {
  int i;
  fprintf(fp,"Grid data structure:\n");
  fprintf(fp,"%d axes\n",a.gdim);
  fprintf(fp,"%d physical axes\n",a.dim);
  for (i=0;i<a.gdim;i++) fprint_axis(fp,a.axes[i]);
}

void print_grid(grid a) { fprint_grid(stdout,a); }

int compare_grid(const grid g1, const grid g2) {
  int i;
  if (g1.gdim != g2.gdim) return 1;
  if (g1.dim  != g2.dim ) return 2;
  for (i=0;i<g1.gdim;i++) { if (compare_axis(g1.axes[i],g2.axes[i])) return 3+i; }
  return 0;
}

int compatible_grid(const grid g1, const grid g2) {
  int i;
  ireal rtest;
  int itest;
  /* compatibility: 
     dim must match
     gdim must match
     d's must match, to within specified tolerance
     difference of o's must be int multiple of d
  */
  if (g1.gdim != g2.gdim) return 1;
  if (g1.dim  != g2.dim ) return 2;
  for (i=0;i<g1.gdim;i++) { 
    if (iwave_abs(g1.axes[i].d-g2.axes[i].d) > TOL*g1.axes[i].d) return 3;
    rtest = g1.axes[i].o-g2.axes[i].o;
    itest = rtest/g1.axes[i].d;
    if ((iwave_abs((itest-1)*g1.axes[i].d-rtest) > TOL*g1.axes[i].d) &&
	(iwave_abs((itest+0)*g1.axes[i].d-rtest) > TOL*g1.axes[i].d) &&
	(iwave_abs((itest+1)*g1.axes[i].d-rtest) > TOL*g1.axes[i].d)) return 4;
  }
  return 0;
}

int get_datasize_grid(grid g) {
  _IPNT _n;
  int i;
  int sz=1;
  get_n(_n,g);
  for (i=0;i<g.dim;i++) sz*=_n[i];
  return sz;
}

int get_global_datasize_grid(grid g) {
  _IPNT _n;
  int i;
  int sz=1;
  get_n(_n,g);
  for (i=0;i<g.gdim;i++) sz*=_n[i];
  return sz;
}

ireal get_cellvol_grid(grid g) {
  ireal v = REAL_ONE;
  int i;
  RPNT _d;
  get_d(_d,g);
  for (i=0;i<g.dim;i++) v*=_d[i];
  return v;
}

ireal get_global_cellvol_grid(grid g) {
  ireal v = REAL_ONE;
  int i;
  RPNT _d;
  get_d(_d,g);
  for (i=0;i<g.gdim;i++) v*=_d[i];
  return v;
}

int get_panelnum_grid(grid g) {
  _IPNT _n;
  int i;
  int sz=1;
  get_n(_n,g);
  for (i=g.dim;i<g.gdim;i++) sz*=_n[i];
  return sz;
}

void get_n(_IPNT n, grid g) {
  int i; 
  for (i=0;i<RARR_MAX_NDIM;i++) {
    n[i]=g.axes[i].n;
  }
}
void get_d(_RPNT d, grid g) {
  int i; 
  for (i=0;i<RARR_MAX_NDIM;i++) {
    d[i]=g.axes[i].d;
  }
}
void get_o(_RPNT o, grid g) {
  int i; 
  for (i=0;i<RARR_MAX_NDIM;i++) {
    o[i]=g.axes[i].o;
  }
}
void get_gs(_IPNT gs, grid g) {
  int i;
  for (i=0;i<RARR_MAX_NDIM;i++) {
    if (g.axes[i].o<0) gs[i]=(int)((g.axes[i].o-g.axes[i].d*TOL)/(g.axes[i].d));
    else gs[i]=(int)((g.axes[i].o+g.axes[i].d*TOL)/(g.axes[i].d));
  }
}
void get_ord(_IPNT od, grid g) {
  int i,j;
  for (i=0;i<RARR_MAX_NDIM;i++) {
    for (j=0;j<RARR_MAX_NDIM;j++) {
      if (g.axes[i].id == j) od[j]=i;
    }
  }
}
