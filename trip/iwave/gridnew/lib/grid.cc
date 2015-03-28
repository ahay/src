#include "grid.h"
#include "except.hh"

void init_default_axis(axis * a) {
  a->n=0;
  a->d=1.0;
  a->o=0.0;
  a->id=-1;
}

void init_axis(axis * a, size_t n, ireal d, ireal o) {
  a->n=n;
  a->d=d;
  a->o=o;
  a->id=-1;
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

void fprint_num_axis(FILE * fp, int i, axis a) {
  fprintf(fp,"axis: n%d=%ld d%d=%e o%d=%e id%d=%d\n",i+1,(long)a.n,i+1,a.d,i+1,a.o,i+1,a.id);
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
  for (i=0;i<IARR_MAX_NDIM;i++) {
    init_default_axis(&(g->axes[i]));
    //    g->axes[i].id=i;
  }
}

int init_grid(grid * g, int dim, int gdim) {
  int i;
  if (dim<1 || dim>gdim || gdim>IARR_MAX_NDIM) return E_BADINPUT;
  g->dim=dim;
  g->gdim=gdim;
  for (i=0;i<IARR_MAX_NDIM;i++) { 
    init_default_axis(&(g->axes[i]));
    //    g->axes[i].id=i;
  }
  return 0;
}

void copy_grid(grid * tgt, const grid * src) {
  int i;
  tgt->dim  = src->dim;
  tgt->gdim = src->gdim;
  for (i=0;i<IARR_MAX_NDIM;i++) 
    copy_axis(&(tgt->axes[i]),&(src->axes[i]));
}

void fprint_grid(FILE * fp, grid a) {
  int i;
  fprintf(fp,"Grid data structure:\n");
  fprintf(fp,"gdim = %d axes\n",a.gdim);
  fprintf(fp,"dim  = %d physical axes\n",a.dim);
  for (i=0;i<a.gdim;i++) fprint_num_axis(fp,i,a.axes[i]);
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
     NOTE: n's do not need to match, nor do axis id's.
  */
  if (g1.gdim != g2.gdim) return 1;
  if (g1.dim  != g2.dim ) return 2;
  for (i=0;i<g1.gdim;i++) { 
    if (iwave_abs(g1.axes[i].d-g2.axes[i].d) > TOL*g1.axes[i].d) return 3;
    rtest = g1.axes[i].o-g2.axes[i].o;
    itest = (int) (rtest/g1.axes[i].d);
    if ((iwave_abs((itest-1)*g1.axes[i].d-rtest) > TOL*g1.axes[i].d) &&
	(iwave_abs((itest+0)*g1.axes[i].d-rtest) > TOL*g1.axes[i].d) &&
	(iwave_abs((itest+1)*g1.axes[i].d-rtest) > TOL*g1.axes[i].d)) return 4;
  }
  return 0;
}

// dimension includes internal extended axes if any - i.e. not same
// as spatial dimension, necessarily
int get_dimension_grid(grid g) {
  _IPNT _id;
  get_id(_id,g);
  int _dim=0;
  for (int i=0;i<IARR_MAX_NDIM;i++) {
    // case 1: physical axes
    if ((-1 < _id[i]) && (_id[i]<g.dim)) _dim++;
    // case 2: internal extended axes
    if ((EXTINT-1 < _id[i])) _dim++;
  }
  // sanity check
  if (_dim < g.dim || _dim > g.gdim) {
    RVL::RVLException e;
    e<<"ERROR: grid::get_dimension_grid\n";
    e<<"  nonsense dim calc: obtained = "<<_dim<<"\n";
    e<<"  grid params:\n";
    e<<"  dim="<<g.dim<<" gdim="<<g.gdim<<"\n";
    for (int i=0;i<IARR_MAX_NDIM;i++) 
      e<<"axis "<<i<<": n="<<g.axes[i].n<<" d="<<g.axes[i].d
       <<" o="<<g.axes[i].o<<" id="<<g.axes[i].id<<"\n";
    throw e;
  }
  return _dim;
}
    
// modified 19.03.15 to account for internal extended axes
// and to properly account for axis id in general
int get_datasize_grid(grid g) {
  _IPNT _n;
  _IPNT _id;
  int i;
  int sz=1;
  int id_check=0;
  get_n(_n,g);
  get_id(_id,g);
  for (i=0;i<IARR_MAX_NDIM;i++) {
    // case 1: physical axes
    if ((-1 < _id[i]) && (_id[i]<g.dim)) { 
      sz*=_n[i];
      id_check++;
    }
  }
  // if no ids or not the right number,
  // recompute using first dim axes
  if (id_check != g.dim) {
    sz=1;
    for (i=0;i<g.dim;i++) sz*=_n[i];
  }
  return sz;
}

// modified 19.03.15 to account for internal extended axes
// and to properly account for axis id in general
int get_extended_datasize_grid(grid g) {
  _IPNT _n;
  _IPNT _id;
  int i;
  int sz=1;
  get_n(_n,g);
  get_id(_id,g);
  for (i=0;i<IARR_MAX_NDIM;i++) {
    // case 1: physical axes
    if ((-1 < _id[i]) && (_id[i]<g.dim)) sz*=_n[i];
    // case 2: internal extended axes
    if ((EXTINT-1 < _id[i])) sz*=_n[i];
  }
  // sanity check: should be at least as many ext'd axes
  // as spatial only - if not (i.e. if size is less) then  
  // proper id info not supplied, and this computation makes
  // no sense
  if (sz < get_datasize_grid(g)) {
    RVL::RVLException e;
    e<<"ERROR: grid::get_extended_datasize_grid\n";
    e<<"  axis id info does not support computation of ext'd size\n";
    e<<"  grid params:\n";
    e<<"  dim="<<g.dim<<" gdim="<<g.gdim<<"\n";
    for (int i=0;i<IARR_MAX_NDIM;i++) 
      e<<"axis "<<i<<": n="<<g.axes[i].n<<" d="<<g.axes[i].d
       <<" o="<<g.axes[i].o<<" id="<<g.axes[i].id<<"\n";
    throw e;  
  }
  return sz;
}

size_t get_global_datasize_grid(grid g) {
  _IPNT _n;
  int i;
  size_t sz=1;
  get_n(_n,g);
  for (i=0;i<g.gdim;i++) sz*=_n[i];
  return sz;
}

// modified 19.03.15 to account for internal extended axes
// and to properly account for axis id in general
ireal get_cellvol_grid(grid g) {
  ireal v = REAL_ONE;
  int i;
  RPNT _d;
  IPNT _id;
  get_d(_d,g);
  get_id(_id,g);
  int id_check = 0;
  for (i=0;i<IARR_MAX_NDIM;i++) {
    // case 1: physical axes
    if ((-1 < _id[i]) && (_id[i]<g.dim)) { 
      v*=_d[i];
      id_check++;
    }
  }
  // recompute if not enough axes id'd as spatial
  // in that case assume that the first dim are the
  // spatial axes
  if (id_check < g.dim) {
    v = REAL_ONE;
    for (i=0;i<g.dim;i++) v*=_d[i];
  }
  return v;
}

// modified 19.03.15 to account for internal extended axes
// and to properly account for axis id in general
ireal get_extended_cellvol_grid(grid g) {
  ireal v = REAL_ONE;
  int i;
  RPNT _d;
  IPNT _id;
  get_d(_d,g);
  get_id(_id,g);
  int id_check = 0;
  for (i=0;i<IARR_MAX_NDIM;i++) {
    // case 1: physical axes
    if ((-1 < _id[i]) && (_id[i]<g.dim)) {
      v*=_d[i];
      id_check++;
    }
    // case 2: internal extended axes
    if ((EXTINT-1 < _id[i])) {
      v*=_d[i];
      id_check++;
    }
  }
  // there must be at least as many ext'd axes as spatial axes
  // if you didn't find that many, then this grid is not
  // adequately decorated
  if (id_check < g.dim) {
    RVL::RVLException e;
    e<<"ERROR: grid::get_extended_datasize_grid\n";
    e<<"  axis id info does not support computation of ext'd size\n";
    e<<"  grid params:\n";
    e<<"  dim="<<g.dim<<" gdim="<<g.gdim<<"\n";
    for (int i=0;i<IARR_MAX_NDIM;i++) 
      e<<"axis "<<i<<": n="<<g.axes[i].n<<" d="<<g.axes[i].d
       <<" o="<<g.axes[i].o<<" id="<<g.axes[i].id<<"\n";
    throw e;  
  }
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

// modified 19.03.15 to account for internal extended axes
// and to properly account for axis id in general
// NOTE: if id's not assigned, will return 0 - impossible
// to count without id's!!!
int get_panelnum_grid(grid g) {
  _IPNT _n;
  _IPNT _id;
  int i;
  int sz=1;
  int id_check=0;
  get_n(_n,g);
  get_id(_id,g);
  for (i=0;i<IARR_MAX_NDIM;i++) {
    // id=dim - time axis id<EXTINT external extended axis
    if ((_id[i] > g.dim) && (_id[i] < EXTINT)) {
      sz*=_n[i];
      id_check++;
    }
  }
  return sz;
}

void get_n(_IPNT n, grid g) {
  int i; 
  for (i=0;i<IARR_MAX_NDIM;i++) {
    n[i]=g.axes[i].n;
  }
}

void get_d(_RPNT d, grid g) {
  int i; 
  for (i=0;i<IARR_MAX_NDIM;i++) {
    d[i]=g.axes[i].d;
  }
}

void get_o(_RPNT o, grid g) {
  int i; 
  for (i=0;i<IARR_MAX_NDIM;i++) {
    o[i]=g.axes[i].o;
  }
}

void get_gs(_IPNT gs, grid g) {
  int i;
  for (i=0;i<IARR_MAX_NDIM;i++) {
    if (g.axes[i].o<0) gs[i]=(int)((g.axes[i].o-g.axes[i].d*TOL)/(g.axes[i].d));
    else gs[i]=(int)((g.axes[i].o+g.axes[i].d*TOL)/(g.axes[i].d));
  }
}

void get_ge(_IPNT ge, grid g) {
  int i;
  IPNT gs;
  IPNT nn;
  get_gs(gs, g);
  get_n(nn, g);
  for (i=0;i<IARR_MAX_NDIM;i++) {
    /*
    if (g.axes[i].o<0) gs[i]=(int)((g.axes[i].o-g.axes[i].d*TOL)/(g.axes[i].d));
    else gs[i]=(int)((g.axes[i].o+g.axes[i].d*TOL)/(g.axes[i].d));
    */
    ge[i]=gs[i]+nn[i]-1;
  }
}

void get_id(_IPNT id, grid g) {
  int i;
  for (i=0;i<IARR_MAX_NDIM;i++) {
    id[i]=g.axes[i].id;
  }
}

void get_ord(_IPNT od, grid g) {
  int i,j;
  for (i=0;i<IARR_MAX_NDIM;i++) {
    for (j=0;j<IARR_MAX_NDIM;j++) {
      if (g.axes[i].id == j) od[j]=i;
    }
  }
}

bool grid_union(grid * g, axis const * ax) {
  for (int i=0; i<g->gdim; i++) {
    if (g->axes[i].id == ax->id) {
      // check for commensurable
      if (fabs(g->axes[i].d-ax->d) > TOL*g->axes[i].d) {
	fprintf(stderr,"d incompat: g = \n");
	fprint_grid(stderr,*g);
	fprintf(stderr,"ax = \n");
	fprint_axis(stderr,*ax);
	return false;
      }
      float dxo = g->axes[i].o - ax->o;
      float dx  = g->axes[i].d;
      if (dxo - dx*((int)((dxo/dx)+2.0f*TOL)) > TOL) {
	fprintf(stderr,"Error: grid_union\n");
	fprintf(stderr,"o incompat: g = \n");
	fprint_grid(stderr,*g);
	fprintf(stderr,"ax = \n");
	fprint_axis(stderr,*ax);
	return false;
      }
      // else, reset min, max, n
      float xmin = iwave_min(g->axes[i].o, ax->o);
      float xmax = iwave_max(g->axes[i].o+(g->axes[i].n-1)*dx, (ax->o+(ax->n-1)*dx));
      g->axes[i].o = xmin;
      g->axes[i].n = 1+(int)((xmax-xmin+TOL)/dx);
      return true;
    }
  }
  // insert new axis in id order
  g->gdim++;
  // fprintf(stderr,"grid_union -> inserting new axis %d\n",g->gdim);
  if (g->gdim > RARR_MAX_NDIM) {
    fprintf(stderr,"Error: grid_union\n");
    fprintf(stderr,"  attempt to exceed max # of axes = %d\n",RARR_MAX_NDIM);
    return false;
  }
  int idx = 0;
  while (g->axes[idx].id < ax->id && g->axes[idx].id > -1) idx++;
  // if at uninit id, just copy
  if (g->axes[idx].id < 0) {
    // fprintf(stderr,"grid_union: copy new axis onto axis %d\n",idx);
    copy_axis(&(g->axes[idx]),ax);
    // fprintf(stderr,"grid_union: resulting grid\n");
    // fprint_grid(stderr,*g);
  }
  else {
    // else shift all subsequent axes out of way
    // fprintf(stderr,"grid_union: gdim=%d new ax id =%d idx=%d\n",g->gdim,ax->id,idx);
    for (int i=g->gdim-1; i>idx; i--) {
      // fprintf(stderr,"grid_union: copy axis %d onto axis %d\n",i-1,i);
      copy_axis(&(g->axes[i]),&(g->axes[i-1]));
    }
    // fprintf(stderr,"grid_union: copy new axis onto axis %d\n",idx);
    copy_axis(&(g->axes[idx]),ax);
    // fprintf(stderr,"grid_union: resulting grid\n");
    // fprint_grid(stderr,*g);
  }
  return true;
}

bool init_step(grid g, IPNT step, bool fwd) {
  // sanity check dim, gdim
  if (!(g.dim<g.gdim) || g.gdim > IARR_MAX_NDIM) {
    return false;
  }
  // the global origin is ex def a member of every grid
  /*
  for (int i=g.dim; i< g.gdim; i++) {
    step[i] = (int)((g.axes[i].o+TOL)/g.axes[i].d);
  }
  */
  get_gs(step, g);
  // time axis - index=dim: last data point if fwd = false (reverse)
  if (!fwd) {
    step[g.dim]+= g.axes[g.dim].n-1;
  }
  return true;
}

bool next_step(grid g, IPNT step) {

  bool more = true;

  IPNT istart;
  IPNT istop;

  get_gs(istart,g);
  get_ge(istop,g);

  for (int i=g.dim+1; i<g.gdim; i++) {

    if (step[i]<istop[i]) {
      step[i]++;
      for (int j=g.dim+1;j<i;j++) {
	step[j]=istart[j];
      }
      return more;
    }
  }
  more = false;
  return more;
}
