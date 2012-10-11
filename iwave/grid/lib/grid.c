#include "grid.h"

int init_default_axis(axis * a) {
  a->n=0;
  a->d=1.0;
  a->o=0.0;
  a->id=0;
  return 0;
}

int init_axis(axis * a, size_t n, ireal d, ireal o) {
  a->n=n;
  a->d=d;
  a->o=o;
  a->id=0;
  return 0;
}

int fprint_axis(FILE * fp, axis a) {
    fprintf(fp,"axis: n=%ld d=%e o=%e id=%d\n",(long) a.n,a.d,a.o,a.id);
  return 0;
}

int print_axis(axis a) {
  return fprint_axis(stdout,a);
}

int compare_axis(axis a1, axis a2) {
  int err=0;
  err = err || (a1.n != a2.n);
  err = err || (fabs((double)(a1.d-a2.d)) > TOL*fmin((double)(a1.d),(double)(a2.d)));
  err = err || (fabs((double)(a1.o-a2.o)) > TOL*fmin((double)(a1.d),(double)(a2.d)));
  err = err || (a1.id != a2.id);
  return err;
}

int init_default_grid(grid * g) {
  int i;
  g->dim=0;
  g->gdim=0;
  for (i=0;i<RARR_MAX_NDIM;i++) {
    init_default_axis(&(g->axes[i]));
    g->axes[i].id=i;
  }
  return 0;
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

int fprint_grid(FILE * fp, grid a) {
  int i;
  fprintf(fp,"Grid data structure:\n");
  fprintf(fp,"%d axes\n",a.gdim);
  fprintf(fp,"%d physical axes\n",a.dim);
  for (i=0;i<a.dim;i++) fprint_axis(fp,a.axes[i]);
  return 0;
}

int print_grid(grid a) { return fprint_grid(stdout,a); }

int compare_grid(grid g1, grid g2) {
  int err=0;
  int i;
  if (g1.gdim != g2.gdim) return 1;
  for (i=0;i<g1.gdim;i++) err = err || compare_axis(g1.axes[i],g2.axes[i]);
  return err;
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

int get_panelnum_grid(grid g) {
  _IPNT _n;
  int i;
  int sz=1;
  get_n(_n,g);
  for (i=g.dim;i<g.gdim;i++) sz*=_n[i];
  return sz;
}

int get_n(_IPNT n, grid g) {
  int i; 
  for (i=0;i<RARR_MAX_NDIM;i++) {
    n[i]=g.axes[i].n;
  }
  return 0;
}
int get_d(_RPNT d, grid g) {
  int i; 
  for (i=0;i<RARR_MAX_NDIM;i++) {
    d[i]=g.axes[i].d;
  }
  return 0;
}
int get_o(_RPNT o, grid g) {
  int i; 
  for (i=0;i<RARR_MAX_NDIM;i++) {
    o[i]=g.axes[i].o;
  }
  return 0;
}
int get_gs(_IPNT gs, grid g) {
  int i;
  for (i=0;i<RARR_MAX_NDIM;i++) {
    if (g.axes[i].o<0) gs[i]=(int)((g.axes[i].o-g.axes[i].d*TOL)/(g.axes[i].d));
    else gs[i]=(int)((g.axes[i].o+g.axes[i].d*TOL)/(g.axes[i].d));
  }
  return 0;
}
int get_ord(_IPNT od, grid g) {
  int i,j;
  for (i=0;i<RARR_MAX_NDIM;i++) {
    for (j=0;j<RARR_MAX_NDIM;j++) {
      if (g.axes[i].id == j) od[j]=i;
    }
  }
  return 0;
}
