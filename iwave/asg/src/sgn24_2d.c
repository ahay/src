/* 
time-step functions
WWS 10.04.11 
After D.S. 08.11.09
After IST 2008-9
********************************************************************************
Implementations of 2-k staggered grid schemes in 1D, 2D, and 3D.
*/

#include "sgn.h"

/* if IPNT ib = IPNT ia + off, then ib indexes the same point in 
   space in RARR b as ia does in RARR a

   thus a
*/
void rel_offset(RARR * a, RARR * b, IPNT off) {
  int idim;
  for (idim = 0; idim<RARR_MAX_NDIM; idim++) 
    off[idim] = a->_dims[idim].gs0 - b->_dims[idim].gs0;
}

/*============================================================================*/

#define NSTORE 20000
static ireal _delta[NSTORE];

//#define ORIG
#undef ORIG

#define VEC
//#undef VEC

#ifdef _OPENMP
#include <omp.h>
#endif /*_OPENMP*/
/*----------------------------------------------------------------------------*/
/*
1/24 constant.
*/
#define C24 4.166666666666666666666667e-2
/*----------------------------------------------------------------------------*/

int sgn_gts2d_24p01sv(RDOM *, RDOM *, RDOM *, void *, int);
int sgn_gts2d_24p01(RDOM *, RDOM *, RDOM *, void *, int);
int sgn_gts2d_24v0(RDOM *, RDOM *, RDOM *, void *, int);
int sgn_gts2d_24v1(RDOM *, RDOM *, RDOM *, void *, int);

int sgn_ts2d_24p01(RDOM *, void *);
int sgn_ts2d_24v0(RDOM *, void *);
int sgn_ts2d_24v1(RDOM *, void *);

// migration to new time step model, 02.12
//int sgn_gts2d_24(RDOM *dom, RDOM *rdom, RDOM *cdom, int iarr, void *pars,int _fwd) {
int sgn_gts2d_24(RDOM *dom, int iarr, void *pars) {

  int _fwd=1;
  SGN_TS_PARS * sgnpars = (SGN_TS_PARS *) pars;  
  /* sanity test - storage mode */
#ifdef STORE
  if ((dom->_s)[D_P0]._dims[0].n > NSTORE)  {
    fprintf(stderr,"NOTIMESTEP 1\n");
    return E_NOTIMESTEP;
  }
#endif
  // Bill's vector farting around
  // 02.12: changed all interfaces back to single domain
#ifndef ORIG
#ifdef VEC
  if ( iarr == D_P0 ) return sgn_gts2d_24p01sv(dom,dom,dom, sgnpars,_fwd);
#else
  if ( iarr == D_P0 ) return sgn_gts2d_24p01(dom,dom,dom, sgnpars,_fwd);
#endif
  if ( iarr == D_V0 ) return sgn_gts2d_24v0(dom, dom, dom, sgnpars, _fwd);
  if ( iarr == D_V1 ) return sgn_gts2d_24v1(dom, dom, dom, sgnpars, _fwd);
#else
  // original Igor code
  if ( iarr == D_P0 ) return sgn_ts2d_24p01(dom, sgnpars);    
  if ( iarr == D_V0 ) return sgn_ts2d_24v0(dom, sgnpars);
  if ( iarr == D_V1 ) return sgn_ts2d_24v1(dom, sgnpars);
#endif
  // in all other cases, no-op
  return 0;
}

/*----------------------------------------------------------------------------*/

#ifdef VEC

int sgn_gts2d_24p01sv(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd) {
  int nx, ny, gxs, gys, iy;
  register ireal * restrict _px, * restrict _py,
    * restrict _mpx, * restrict _epx, * restrict _epy, * restrict _vx3,
    * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0;
  register ireal lax, lay, dt2, etaydt;
  //  register int i,j;
  register int i;
  register ireal *tmp;  /*pointer for swapping p with mp during imaging accumulation*/ 
  register ireal *_rmpx;  /*pointer for stroe scaling multipliers*/ 
  //register int * ia, *ja;
  //ireal _ddel;
  //int _mi, _mrp, _mp0;

  RARR *s, *rs, *cs;
  //  ireal delta __attribute__ ((__aligned__(16)));
  //int err = 0;
  //fprintf(stderr,"sgn_gts2d_24p01,VEC, _fwd %d, err %d\n ",_fwd, err);
  //scanf("%d", &err); 
  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
  
  nx = s[D_P0]._dims[0].n;
  ny = s[D_P0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
  
  gxs = s[D_P0]._dims[0].gs;
  gys = s[D_P0]._dims[1].gs;
  
  lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;
  /* commont out by D.S. 01.05.11 */
  /*
    if (!_fwd) {
    lax=-lax;
    lay=-lay;
  }
  */
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
  
  _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs) * s[D_P0 ]._dims[0].n0;
  _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs) * s[D_P1 ]._dims[0].n0;
  _mpx = cs[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs) * s[D_MP0]._dims[0].n0;
  _vx3 = rs[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs) * s[D_V0 ]._dims[0].n0 + 1;
  _vy2 = rs[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs) * s[D_V1 ]._dims[0].n0;
  _vy3 = _vy2 + s[D_V1]._dims[0].n0;
  _vy1 = _vy2 - s[D_V1]._dims[0].n0;
  _vy0 = _vy1 - s[D_V1]._dims[0].n0;
  _epx = s[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);              /* 1D */
  _epy = s[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs);        /* 1D */
  
  // this will work fine for pert about scalars - what does it mean in general?
  _rmpx = rs[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs) * s[D_MP0]._dims[0].n0;

  if (!_fwd) {  printf("mpx:\n");
    for (i=0;i<10;i++) 
      printf("mpx[%d]=%20.14e\n",i,*(_rmpx+241*s[D_MP0]._dims[0].n0 + i));
  }
  
  for ( iy = 0; iy < ny; iy ++) {
    etaydt = (*_epy) * dt2;
    
    /* swap pointers when _fwd = 0 (imaging accumulation)*/
    if(!_fwd) { 
      tmp = _px;
      _px=_mpx; 
      _mpx=tmp;
    }
    
    for ( i = 0; i < nx; ++i )  
      _delta[i] = ((_vx3[i] - _vx3[i-3] + (_vx3[i-2] - _vx3[i-1]) * (ireal)27.0) * lax + 
		   (_vy3[i] - _vy0[i  ] + (_vy1[i  ] - _vy2[i  ]) * (ireal)27.0) * lay);
    
    /* analysis:
       p = (a + b)/c
       only a depends on p
       so can write
       p *= (1+epx*dt2);
       p += delta*mpx;
       p /= 1+epx*dt2;
       for operator version of mpx: (note that ia deps on both indices!)
       for (j=ia[i];j<ia[i+1];j++) p+=delta[i-j-ia[i]+ja[i]]*mpx[j];
    */
    
    if(_fwd)           
      /*
      if (ia && ja) {
	for ( i = 0; i < nx; ++i ) { 
	  _px[i] *= (1.0f - _epx[i] * dt2);
	  _py[i] *= (1.0f - etaydt);
	  _ddel=0.0f;
	  for (j=ia[_mi+i];j<ia[_mi+i+1];j++) 
	    _ddel+=_delta[i-j+ja[_mi+i]]*_mpx[j];
	  _px[i] = (_px[i]+_ddel)/(1.0f+_epx[i]*dt2);
	  _py[i] = (_py[i]+_ddel)/(1.0f + etaydt);	 
	}
      }

      else {
      */
	for ( i = 0; i < nx; ++i ) { 
	  //_px[i] = _px[i]  + _delta[i] * _mpx[i] / ((ireal)1.0 + etaxdt);
	  _px[i] = (_px[i] * (1.0f - _epx[i] * dt2) + _delta[i]*_mpx[i]) / (1.0f + _epx[i] * dt2);
	  //_py[i] = _py[i]  + _delta[i]*_mpx[i] / ((ireal)1.0 + etaydt);
	  _py[i] = (_py[i] * (1.0f - etaydt) + _delta[i]*_mpx[i]) / (1.0f + etaydt);	 
	}
    //      }	
    else
      for ( i = 0; i < nx; ++i ) 
	_px[i] = _px[i]  + _delta[i]*_mpx[i]/_rmpx[i];
    // _px[i] = (_px[i] * ((ireal)1.0 - etaxdt) + _delta[i]/_rmpx[i]) / ((ireal)1.0 + etaxdt);
    
    /* swap pointers back */
    if(!_fwd){ 
      tmp = _mpx;
      _mpx=_px; 
      _px=tmp;
    }
    
    _px  += s[D_P0 ]._dims[0].n0;
    _py  += s[D_P1 ]._dims[0].n0; 
    /*
    if (ia && ja) {
      _mrp += s[D_MP0]._dims[0].n0;
      _mi= ia[_mp0 + _mrp]- ia[_mp0];
      _mpx = cs[D_MP0]._s + _mi;
    }
    else
    */ 
      _mpx += s[D_MP0]._dims[0].n0;
    _vx3 += s[D_V0 ]._dims[0].n0; 
    _vy0 += s[D_V1 ]._dims[0].n0;
    _vy1 += s[D_V1 ]._dims[0].n0; 
    _vy2 += s[D_V1 ]._dims[0].n0; 
    _vy3 += s[D_V1 ]._dims[0].n0; 
    _epy ++;
  }

  return 0;
}

int sgn_gts2d_24p01(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd) {
  int nx, ny, gxs, gys, iy, tsz, tid;
  register ireal * restrict _px, * restrict _py,
    * restrict _mpx, * restrict _epx, * restrict _epy, * restrict _vx3,
    * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0;
  register ireal lax, lay, dt2, etaxdt, etaydt;
  register int i;
  register ireal *tmp;  /*pointer for swapping p with mp during imaging accumulation*/ 
  register ireal *_rmpx;  /*pointer for stroe scaling multipliers*/ 
  RARR *s, *rs, *cs;
#ifndef STORE
  ireal delta __attribute__ ((__aligned__(16)));
#endif
  //int err = 0;
  //fprintf(stderr,"sgn_gts2d_24p01,VEC, _fwd %d, err %d\n ",_fwd, err);
  //scanf("%d", &err); 
  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
  
  nx = s[D_P0]._dims[0].n;
  ny = s[D_P0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
  
  gxs = s[D_P0]._dims[0].gs;
  gys = s[D_P0]._dims[1].gs;
  
  lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
  lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
  
#pragma omp parallel private(tsz,tid,iy,_px,_py,_mpx,_epx,_epy,_vx3,_vy3,_vy2,_vy1,_vy0,delta,etaxdt,etaydt,i)
  {
#ifdef _OPENMP
    tsz = omp_get_num_threads();
    tid = omp_get_thread_num();
#else
    tsz = 1;
    tid = 0;
#endif
    
    _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0;
    _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
    _mpx = cs[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs + tid) * s[D_MP0]._dims[0].n0;
    _vx3 = rs[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0 + 1;
    _vy2 = rs[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
    _vy3 = _vy2 + s[D_V1]._dims[0].n0;
    _vy1 = _vy2 - s[D_V1]._dims[0].n0;
    _vy0 = _vy1 - s[D_V1]._dims[0].n0;
    _epx = s[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);              /* 1D */
    _epy = s[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs + tid);        /* 1D */
    
    _rmpx = rs[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs + tid) * s[D_MP0]._dims[0].n0;

    for ( iy = tid; iy < ny; iy += tsz )
      {
	etaydt = (*_epy) * dt2;

	/* swap pointers when _fwd = 0 (imaging accumulation)*/
	if(!_fwd) { 
	  tmp = _px;
	  _px=_mpx; 
	  _mpx=tmp;
	}
	
#ifdef STORE

	for ( i = 0; i < nx; ++i )
	  {
	    _delta[i] = ((_vx3[i] - _vx3[i-3] + (_vx3[i-2] - _vx3[i-1]) * (ireal)27.0) * lax + 
			(_vy3[i] - _vy0[i  ] + (_vy1[i  ] - _vy2[i  ]) * (ireal)27.0) * lay) * _mpx[i];
	  }
	for ( i = 0; i < nx; ++i )
	  {
	    etaxdt = _epx[i] * dt2; 
	    if(_fwd)           
	      _px[i] = _px[i]  + _delta[i] / ((ireal)1.0 + etaxdt);
	      //_px[i] = (_px[i] * ((ireal)1.0 - etaxdt) + _delta[i]) / ((ireal)1.0 + etaxdt);
	    else {
	      _px[i] = _px[i]  + _delta[i]/_rmpx[i];
	      // _px[i] = (_px[i] * ((ireal)1.0 - etaxdt) + _delta[i]/_rmpx[i]) / ((ireal)1.0 + etaxdt);
	    }
	  }                
	if(_fwd){
	  for ( i = 0; i < nx; ++i )
	    {
	      _py[i] = _py[i]  + _delta[i] / ((ireal)1.0 + etaydt);
	      //_py[i] = (_py[i] * ((ireal)1.0 - etaydt) + _delta[i]) / ((ireal)1.0 + etaydt);
	    }
	}
#else            
	for ( i = 0; i < nx; ++i )
	  {
	    etaxdt = _epx[i] * dt2;
	    delta = ((_vx3[i] - _vx3[i-3] + (_vx3[i-2] - _vx3[i-1]) * (ireal)27.0) * lax + 
		     (_vy3[i] - _vy0[i  ] + (_vy1[i  ] - _vy2[i  ]) * (ireal)27.0) * lay) * _mpx[i];
	    if(_fwd)
	      _px[i] = _px[i]  + delta/ ((ireal)1.0 + etaxdt);
	    //_px[i] = (_px[i] * ((ireal)1.0 - etaxdt) + delta) / ((ireal)1.0 + etaxdt);
	    else {
	      _px[i] = _px[i]  + delta/_rmpx[i];
	      // _px[i] = (_px[i] * ((ireal)1.0 - etaxdt) + delta/_rmpx[i]) / ((ireal)1.0 + etaxdt);
	    }
	  }                
	if(_fwd){
	  for ( i = 0; i < nx; ++i )
	    {
	      delta = ((_vx3[i] - _vx3[i-3] + (_vx3[i-2] - _vx3[i-1]) * (ireal)27.0) * lax + 
		       (_vy3[i] - _vy0[i  ] + (_vy1[i  ] - _vy2[i  ]) * (ireal)27.0) * lay) * _mpx[i];
	      
	      _py[i] = _py[i] + delta/ ((ireal)1.0 + etaydt);
	      // _py[i] = (_py[i] * ((ireal)1.0 - etaydt) + delta) / ((ireal)1.0 + etaydt);
	    }
	}
#endif
	/* swap pointers back */
	if(!_fwd){ 
	  tmp = _mpx;
	  _mpx=_px; 
	  _px=tmp;
	}
	
	_px  += tsz * s[D_P0 ]._dims[0].n0;
	_py  += tsz * s[D_P1 ]._dims[0].n0; 
	_mpx += tsz * s[D_MP0]._dims[0].n0;
	_vx3 += tsz * s[D_V0 ]._dims[0].n0; 
	_vy0 += tsz * s[D_V1 ]._dims[0].n0;
	_vy1 += tsz * s[D_V1 ]._dims[0].n0; 
	_vy2 += tsz * s[D_V1 ]._dims[0].n0; 
	_vy3 += tsz * s[D_V1 ]._dims[0].n0; 
	_epy += tsz;
      }
  } /* omp parallel */

  return 0;
}

/*----------------------------------------------------------------------------*/

int sgn_gts2d_24v0(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
  int nx, ny, gxs, gys, iy;
  register ireal * restrict _vx, * restrict _mvx,
    * restrict _evx, * restrict _px3;
  register ireal lax, dt2;
  register ireal *tmp;  /*pointer for swapping p with mp during imaging accumulation*/ 
  register ireal *_rmvx;  /*pointer for stroe scaling multipliers*/ 
  register int i; 
  RARR *s, *rs, *cs;
  //int err = 0;

  //scanf("%d", &err); 
     
  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;

  lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
  /* add by D.S. 03.05.11 */
  if (!_fwd) {
    lax=-lax;
  }
  /*   
  fprintf(stderr,"\n*************** enter sgn_gts2d_24v0,VEC, _fwd=%d lax=%g dt2=%g\n",_fwd,lax,dt2);
  fprintf(stderr,"OUTPUT D_V0:\n");
  ra_print(&(s[D_V0]),stderr);
  fprintf(stderr,"INPUT D_P0:\n");
  ra_print(&(rs[D_P0]),stderr);
  fprintf(stderr,"MULTIPLIER D_MV0:\n");
  ra_print(&(cs[D_MV0]),stderr);
  fprintf(stderr,"\n*************** sgn_gts2d_24v0 loops \n",_fwd);
  */
  nx = s[D_V0]._dims[0].n;
  ny = s[D_V0]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
    
  gxs = s[D_V0]._dims[0].gs;
  gys = s[D_V0]._dims[1].gs;
    
   _vx  = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs) * s[D_V0 ]._dims[0].n0;
  _mvx = cs[D_MV0]._s + (gxs - s[D_MV0]._dims[0].gs) + (gys - s[D_MV0]._dims[1].gs) * s[D_MV0]._dims[0].n0;
  _px3 = rs[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs) * s[D_P0 ]._dims[0].n0 + 2;
  _evx = rs[D_EV0]._s + (gxs - s[D_EV0]._dims[0].gs);              /* 1D */
  _rmvx = rs[D_MV0]._s + (gxs - s[D_MV0]._dims[0].gs) + (gys - s[D_MV0]._dims[1].gs) * s[D_MV0]._dims[0].n0;      

  for ( iy = 0; iy < ny; iy ++ ) {            
    /* swap pointers when _fwd = 0 (imaging accumulation)*/
    if(!_fwd) { 
      tmp = _vx;
      _vx = _mvx; 
      _mvx = tmp;
    }
    
    if(_fwd)
      for ( i = 0;  i < nx; ++i ) 
	//	      _vx[i] = _vx[i] + (_px3[i] - _px3[i-3] + (_px3[i-2] - _px3[i-1]) * (ireal)27.0) * lax * _mvx[i]/((ireal)1.0 + etaxdt);    
	_vx[i] = (_vx[i] * (1.0f - _evx[i] * dt2) + (_px3[i] - _px3[i-3] + (_px3[i-2] - _px3[i-1]) * 27.0f) * lax * _mvx[i]) / (1.0f + _evx[i] * dt2);
    else 
      for ( i = 0;  i < nx; ++i ) 
	//	      _vx[i] = _vx[i]  + (_px3[i] - _px3[i-3] + (_px3[i-2] - _px3[i-1]) * (ireal)27.0) * lax * _mvx[i] /_rmvx[i];
	_vx[i] = (_vx[i] * (1.0f - _evx[i] * dt2) + (_px3[i] - _px3[i-3] + (_px3[i-2] - _px3[i-1]) * 27.0f) * lax * _mvx[i]/_rmvx[i]) / (1.0f + _evx[i] * dt2);
    
    /* swap pointers back */
    if(!_fwd) { 
      tmp = _mvx;
      _mvx=_vx; 
      _vx=tmp;
    }
    _vx += s[D_V0 ]._dims[0].n0;
    _mvx += s[D_MV0]._dims[0].n0;
    _px3 += s[D_P0 ]._dims[0].n0;

  } 

  /*
  fprintf(stderr,"\n*************** exit sgn_gts2d_24v0,VEC, _fwd %d \n",_fwd);
  fprintf(stderr,"OUTPUT D_V0:\n");
  ra_print(&(s[D_V0]),stderr);
  fprintf(stderr,"INPUT D_P0:\n");
  ra_print(&(rs[D_P0]),stderr);
  fprintf(stderr,"MULTIPLIER D_MV0:\n");
  ra_print(&(cs[D_MV0]),stderr);
  fprintf(stderr,"\n*************** sgn_gts2d_24v0 exit \n",_fwd);
  */

  return 0;
}

/*----------------------------------------------------------------------------*/

int sgn_gts2d_24v1(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
  int nx, ny, gxs, gys, iy;
  register ireal * restrict _vy, * restrict _mvy, * restrict _evy,
    * restrict _py3, * restrict _py2, * restrict _py1, * restrict _py0;
  register ireal lay, dt2, etaydt;
  register int i;
  register ireal *tmp;  /*pointer for swapping p with mp during imaging accumulation*/ 
  register ireal *_rmvy;  /*pointer for stroe scaling multipliers*/ 
  RARR *s, *rs, *cs;
  //int err = 0;
  //fprintf(stderr,"sgn_gts2d_24v1,VEC _fwd %d, err %d\n ",_fwd, err);

  s = dom->_s;
  rs = rdom->_s;
  cs = cdom->_s;
    
  nx = s[D_V1]._dims[0].n;
  ny = s[D_V1]._dims[1].n;
  if ( nx * ny == 0 ) return 0;
    
  gxs = s[D_V1]._dims[0].gs;
  gys = s[D_V1]._dims[1].gs;
    
  lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;
  dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
  /* add by D.S. 03.05.11 */
  if (!_fwd) {
    lay=-lay;
  } 
     
  _vy  = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs) * s[D_V1 ]._dims[0].n0;
  _mvy = cs[D_MV1]._s + (gxs - s[D_MV1]._dims[0].gs) + (gys - s[D_MV1]._dims[1].gs) * s[D_MV1]._dims[0].n0;
  _py1 = rs[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs) * s[D_P1 ]._dims[0].n0;
  _py2 = _py1 + s[D_P1]._dims[0].n0;
  _py3 = _py2 + s[D_P1]._dims[0].n0;
  _py0 = _py1 - s[D_P1]._dims[0].n0;
  _evy = rs[D_EV1]._s + (gys - s[D_EV1]._dims[0].gs);        /* 1D */
  
  _rmvy = rs[D_MV1]._s + (gxs - s[D_MV1]._dims[0].gs) + (gys - s[D_MV1]._dims[1].gs) * s[D_MV1]._dims[0].n0;
  
  for ( iy = 0; iy < ny; iy ++) {            
      
    etaydt = (*_evy) * dt2;
    /* swap pointers when _fwd = 0 (imaging accumulation)*/
    if(!_fwd) { 
      tmp = _vy;
      _vy = _mvy; 
      _mvy = tmp;
    }
    
    if(_fwd)
      for ( i = 0; i < nx; ++i ) 
	//  _vy[i] = _vy[i] + (_py3[i] - _py0[i] + (_py1[i] - _py2[i]) * (ireal)27.0) * lay * _mvy[i]/ ((ireal)1.0 + etaydt);
	_vy[i] = (_vy[i] * (1.0f - etaydt) + (_py3[i] - _py0[i] + (_py1[i] - _py2[i]) * 27.0f) * lay * _mvy[i]) / (1.0f + etaydt);
    else 
      for ( i = 0; i < nx; ++i ) 
	//  _vy[i] = _vy[i]  + (_py3[i] - _py0[i] + (_py1[i] - _py2[i]) * (ireal)27.0) * lay * _mvy[i]/_rmvy[i];
	_vy[i] = (_vy[i] * (1.0 - etaydt) + (_py3[i] - _py0[i] + (_py1[i] - _py2[i]) * 27.0f) * lay * _mvy[i]/_rmvy[i]) / (1.0f + etaydt);
    
    /* swap pointers back */
    if(!_fwd) { 
      tmp = _mvy;
      _mvy = _vy; 
      _vy = tmp;
    }
    
    _vy += s[D_V1 ]._dims[0].n0;
    _mvy += s[D_MV1]._dims[0].n0;
    _py0 += s[D_P1 ]._dims[0].n0;
    _py1 += s[D_P1 ]._dims[0].n0;
    _py2 += s[D_P1 ]._dims[0].n0;
    _py3 += s[D_P1 ]._dims[0].n0;
    _evy ++;
    
  }

  return 0;
}
/*---- END VECTOR BRANCH -----------------------------------------------------*/
/*---- BEGIN POINTER BRANCH --------------------------------------------------*/

#else 

/*----------------------------------------------------------------------------*/

int sgn_gts2d_24p01(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
    int nx, ny, gxs, gys, px_a, py_a, mpx_a, vx_a, vy_a, iy, tsz, tid;
    register ireal * restrict _px, * restrict _py, * restrict _pxend,
                  * restrict _mpx, * restrict _epx, * restrict _epy, * restrict _vx3,
                  * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0;
    register ireal lax, lay, dt2, vx3, vx2, vx1, vx0, delta, etaxdt, etaydt;
    register ireal *tmp;  /*pointer for swapping p with mp during imaging accumulation*/ 
    register ireal *_rmpx;  /*pointer for stroe scaling multipliers*/ 
    RARR *s, *rs, *cs;
    //int err = 0;
    //fprintf(stderr,"sgn_gts2d_24p01,POINTER, _fwd %d, err %d\n ",_fwd, err);
    //scanf("%d", &err); 
    
    s = dom->_s;
    rs = rdom->_s;
    cs = cdom->_s;

    nx = s[D_P0]._dims[0].n;
    ny = s[D_P0]._dims[1].n;
    if ( nx * ny == 0 ) return 0;
            
    gxs = s[D_P0]._dims[0].gs;
    gys = s[D_P0]._dims[1].gs;
    
    lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
    lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;
    dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;

    #pragma omp parallel private(tsz,tid,iy,_px,_py,_pxend,_mpx,_epx,_epy,_vx3,_vy3,_vy2,_vy1,_vy0,vx3,vx2,vx1,vx0,delta,etaxdt,etaydt)
    {
        #ifdef _OPENMP
        tsz = omp_get_num_threads();
        tid = omp_get_thread_num();
        #else
        tsz = 1;
        tid = 0;
        #endif
        
        #pragma omp single
        {            
            px_a  = tsz * s[D_P0 ]._dims[0].n0 - nx;
            py_a  = tsz * s[D_P1 ]._dims[0].n0 - nx;
            mpx_a = tsz * s[D_MP0]._dims[0].n0 - nx;
            vx_a  = tsz * s[D_V0 ]._dims[0].n0 - nx;
            vy_a  = tsz * s[D_V1 ]._dims[0].n0 - nx;
        }
        #pragma omp barrier

        _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0;
        _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
        _mpx = cs[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs + tid) * s[D_MP0]._dims[0].n0;
        _vx3 = rs[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0 + 1;
        _vy2 = rs[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
        _vy3 = _vy2 + s[D_V1]._dims[0].n0;
        _vy1 = _vy2 - s[D_V1]._dims[0].n0;
        _vy0 = _vy1 - s[D_V1]._dims[0].n0;
        _epx = rs[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);              /* 1D */
        _epy = rs[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs + tid);        /* 1D */

	_rmpx = rs[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs + tid) * s[D_MP0]._dims[0].n0;

        for ( iy = tid; iy < ny; iy += tsz )
        {
            vx2 = _vx3[-1]; vx1 = _vx3[-2]; vx0 = _vx3[-3];
            etaydt = (*_epy) * dt2;
            
	    /* swap pointers when _fwd = 0 (imaging accumulation)*/
	    if(!_fwd) { 
	      tmp = _px;
	      _px=_mpx; 
	      _mpx=tmp;
	    }

            for ( _pxend = _px + nx; _px < _pxend; )
            {
                vx3 = *_vx3++;
                etaxdt = (*_epx++) * dt2;
                
                delta = ((vx3 - vx0 + (vx1 - vx2) * 27.0) * lax + 
                         ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay) * (*_mpx++);
		if(_fwd)
		  (*_px) = (*_px) + delta/(1.0 + etaxdt);
		  // (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
		else {
		  (*_px) = (*_px) + delta/(*_rmpx++);
		  // (*_px) = ((*_px) * (1.0 - etaxdt) + delta/(*_rmpx++)) / (1.0 + etaxdt);
		}
                _px++;
		
                if(_fwd){
		  (*_py) = (*_py) + delta/(1.0 + etaydt);
		  // (*_py) = ((*_py) * (1.0 - etaydt) + delta) / (1.0 + etaydt);
		  _py++;
                }

                vx0 = vx1; vx1 = vx2; vx2 = vx3;
            }
	    /* swap pointers back */
	    if(!_fwd){ 
	      tmp = _mpx;
	      _mpx=_px; 
	      _px=tmp;
	    }
            
            _px += px_a; _py += py_a; _mpx += mpx_a; _vx3 += vx_a; _vy0 += vy_a;
            _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; _epx -= nx; _epy += tsz;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_gts2d_24v0(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
    int nx, ny, gxs, gys, vx_a, mvx_a, px_a, iy, tsz, tid;
    register ireal * restrict _vx, * restrict _vxend, * restrict _mvx,
                  * restrict _evx, * restrict _px3;
    register ireal lax, dt2, px3, px2, px1, px0, etaxdt;
    register ireal *tmp;  /*pointer for swapping p with mp during imaging accumulation*/ 
    register ireal *_rmvx;  /*pointer for stroe scaling multipliers*/ 
    RARR *s, *rs, *cs;
    //int err = 0;
    //fprintf(stderr,"sgn_gts2d_24v0,POINTER, _fwd %d, err %d\n ",_fwd, err);
    //scanf("%d", &err); 
    

    s = dom->_s;
    rs = rdom->_s;
    cs = cdom->_s;
   
    nx = s[D_V0]._dims[0].n;
    ny = s[D_V0]._dims[1].n;
    if ( nx * ny == 0 ) return 0;
    
    gxs = s[D_V0]._dims[0].gs;
    gys = s[D_V0]._dims[1].gs;
    
    lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
    dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
    #pragma omp parallel private(tsz,tid,iy,_vx,_vxend,_mvx,_evx,_px3,px3,px2,px1,px0,etaxdt)
    {
        #ifdef _OPENMP
        tsz = omp_get_num_threads();
        tid = omp_get_thread_num();
        #else
        tsz = 1;
        tid = 0;
        #endif
        
        #pragma omp single
        {
            vx_a  = tsz * s[D_V0 ]._dims[0].n0 - nx;
            mvx_a = tsz * s[D_MV0]._dims[0].n0 - nx;
            px_a  = tsz * s[D_P0 ]._dims[0].n0 - nx;
        }
        #pragma omp barrier
        
        _vx  = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0;
        _mvx = cs[D_MV0]._s + (gxs - s[D_MV0]._dims[0].gs) + (gys - s[D_MV0]._dims[1].gs + tid) * s[D_MV0]._dims[0].n0;
        _px3 = rs[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0 + 2;
        _evx = rs[D_EV0]._s + (gxs - s[D_EV0]._dims[0].gs);              /* 1D */
       
	_rmvx = rs[D_MV0]._s + (gxs - s[D_MV0]._dims[0].gs) + (gys - s[D_MV0]._dims[1].gs + tid) * s[D_MV0]._dims[0].n0;
 
        for ( iy = tid; iy < ny; iy += tsz )
        {
            px2 = _px3[-1]; px1 = _px3[-2]; px0 = _px3[-3];
            
	    /* swap pointers when _fwd = 0 (imaging accumulation)*/
	    if(!_fwd) { 
	      tmp = _vx;
	      _vx = _mvx; 
	      _mvx = tmp;
	    }

            for ( _vxend = _vx + nx; _vx < _vxend; )
            {
                px3 = *_px3++;
                etaxdt = (*_evx++) * dt2;
                if(_fwd)
		  (*_vx) = (*_vx) + (px3 - px0 + (px1 - px2) * 27.0) * lax * (*_mvx++)/(1.0 + etaxdt);
		  // (*_vx) = ((*_vx) * (1.0 - etaxdt) + (px3 - px0 + (px1 - px2) * 27.0) * lax * (*_mvx++)) / (1.0 + etaxdt);
		else {
		  (*_vx) = (*_vx) + (px3 - px0 + (px1 - px2) * 27.0) * lax * (*_mvx++) / (*_rmvx++);
		  //(*_vx) = ((*_vx) * (1.0 - etaxdt) + (px3 - px0 + (px1 - px2) * 27.0) * lax * (*_mvx++)/(*_rmvx++)) / (1.0 + etaxdt);
		}
                _vx++;
                
                px0 = px1; px1 = px2; px2 = px3;
            }
	    /* swap pointers back */
	    if(!_fwd) { 
	      tmp = _mvx;
	      _mvx=_vx; 
	      _vx=tmp;
	    }
            _vx += vx_a; _mvx += mvx_a; _px3 += px_a; _evx -= nx;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_gts2d_24v1(RDOM *dom, RDOM *rdom, RDOM *cdom, void *pars, int _fwd)
{
    int nx, ny, gxs, gys, vy_a, mvy_a, py_a, iy, tsz, tid;
    register ireal * restrict _vy, * restrict _vyend, * restrict _mvy, * restrict _evy,
                  * restrict _py3, * restrict _py2, * restrict _py1, * restrict _py0;
    register ireal lay, dt2, etaydt;
    register ireal *tmp;  /*pointer for swapping p with mp during imaging accumulation*/ 
    register ireal *_rmvy;  /*pointer for stroe scaling multipliers*/ 
    RARR *s, *rs, *cs;
    //int err = 0;
    //fprintf(stderr,"sgn_gts2d_24v1,POINTER, _fwd %d, err %d\n ",_fwd, err);
    //scanf("%d", &err); 
    

    s = dom->_s;
    rs = rdom->_s;
    cs = cdom->_s;
       
    nx = s[D_V1]._dims[0].n;
    ny = s[D_V1]._dims[1].n;
    if ( nx * ny == 0 ) return 0;
    
    gxs = s[D_V1]._dims[0].gs;
    gys = s[D_V1]._dims[1].gs;
    
    lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;
    dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
    #pragma omp parallel private(tsz,tid,iy,_vy,_vyend,_mvy,_evy,_py3,_py2,_py1,_py0,etaydt)
    {
        #ifdef _OPENMP
        tsz = omp_get_num_threads();
        tid = omp_get_thread_num();
        #else
        tsz = 1;
        tid = 0;
        #endif
        
        #pragma omp single
        {
            vy_a  = tsz * s[D_V1 ]._dims[0].n0 - nx;
            mvy_a = tsz * s[D_MV1]._dims[0].n0 - nx;
            py_a  = tsz * s[D_P1 ]._dims[0].n0 - nx;
        }
        #pragma omp barrier
        
        _vy  = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
        _mvy = cs[D_MV1]._s + (gxs - s[D_MV1]._dims[0].gs) + (gys - s[D_MV1]._dims[1].gs + tid) * s[D_MV1]._dims[0].n0;
        _py1 = rs[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
        _py2 = _py1 + s[D_P1]._dims[0].n0;
        _py3 = _py2 + s[D_P1]._dims[0].n0;
        _py0 = _py1 - s[D_P1]._dims[0].n0;
        _evy = rs[D_EV1]._s + (gys - s[D_EV1]._dims[0].gs + tid);        /* 1D */
     
	_rmvy = rs[D_MV1]._s + (gxs - s[D_MV1]._dims[0].gs) + (gys - s[D_MV1]._dims[1].gs + tid) * s[D_MV1]._dims[0].n0;
   
        for ( iy = tid; iy < ny; iy += tsz )
        {            
            etaydt = (*_evy) * dt2;
            /* swap pointers when _fwd = 0 (imaging accumulation)*/
	    if(!_fwd) { 
	      tmp = _vy;
	      _vy = _mvy; 
	      _mvy = tmp;
	    }
	  
            for ( _vyend = _vy + nx; _vy < _vyend; )
            {
	      if(_fwd)
		(*_vy) = (*_vy) + ((*_py3++) - (*_py0++) +  ((*_py1++) - (*_py2++)) * 27.0) * lay * (*_mvy++)/(1.0 + etaydt); 
	      // (*_vy) = ((*_vy) * (1.0 - etaydt) + ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay * (*_mvy++)) / (1.0 + etaydt);
	      else {
		(*_vy) = (*_vy) + ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay * (*_mvy++)/(*_rmvy++);
		// (*_vy) = ((*_vy) * (1.0 - etaydt) + ((*_py3++) - (*_py0++) + ((*_py1++) - (*_py2++)) * 27.0) * lay * (*_mvy++)/(*_rmvy++)) / (1.0 + etaydt);
	      }
                _vy++;
            }

	    /* swap pointers back */
	    if(!_fwd) { 
	      tmp = _mvy;
	      _mvy = _vy; 
	      _vy = tmp;
	    }
	    
            _vy += vy_a; _mvy += mvy_a; _py0 += py_a; _py1 += py_a; _py2 += py_a; _py3 += py_a; _evy += tsz;
        }
    } /* omp parallel */

    return 0;
}
/*---- END POINTER BRANCH ----------------------------------------------------*/
#endif

int asg_ts2k(RDOM * unext, RDOM * ucurr, RDOM * coef, int iarr, void *pars) {
  return E_NOTIMESTEP;
}

/*********************************** original sgn functions ****************************************/

/*----------------------------------------------------------------------------*/

int sgn_ts2d_24p01(RDOM *dom, void *pars)
{
    int nx, ny, gxs, gys, px_a, py_a, mpx_a, vx_a, vy_a, iy, tsz, tid;
    register ireal * restrict _px, * restrict _py, * restrict _pxend,
                  * restrict _mpx, * restrict _epx, * restrict _epy, * restrict _vx3,
                  * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0;
    register ireal lax, lay, dt2, vx3, vx2, vx1, vx0, delta, etaxdt, etaydt;
    RARR *s;
    
    s = dom->_s;

    nx = s[D_P0]._dims[0].n;
    ny = s[D_P0]._dims[1].n;
    if ( nx * ny == 0 ) return 0;
            
    gxs = s[D_P0]._dims[0].gs;
    gys = s[D_P0]._dims[1].gs;
    
    lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
    lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;
    dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;

    #pragma omp parallel private(tsz,tid,iy,_px,_py,_pxend,_mpx,_epx,_epy,_vx3,_vy3,_vy2,_vy1,_vy0,vx3,vx2,vx1,vx0,delta,etaxdt,etaydt)
    {
        #ifdef _OPENMP
        tsz = omp_get_num_threads();
        tid = omp_get_thread_num();
        #else
        tsz = 1;
        tid = 0;
        #endif
        
        #pragma omp single
        {            
            px_a  = tsz * s[D_P0 ]._dims[0].n0 - nx;
            py_a  = tsz * s[D_P1 ]._dims[0].n0 - nx;
            mpx_a = tsz * s[D_MP0]._dims[0].n0 - nx;
            vx_a  = tsz * s[D_V0 ]._dims[0].n0 - nx;
            vy_a  = tsz * s[D_V1 ]._dims[0].n0 - nx;
        }
        #pragma omp barrier

        _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0;
        _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
        _mpx = s[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs + tid) * s[D_MP0]._dims[0].n0;
        _vx3 = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0 + 1;
        _vy2 = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
        _vy3 = _vy2 + s[D_V1]._dims[0].n0;
        _vy1 = _vy2 - s[D_V1]._dims[0].n0;
        _vy0 = _vy1 - s[D_V1]._dims[0].n0;
        _epx = s[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);              /* 1D */
        _epy = s[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs + tid);        /* 1D */

        for ( iy = tid; iy < ny; iy += tsz )
        {
            vx2 = _vx3[-1]; vx1 = _vx3[-2]; vx0 = _vx3[-3];
            etaydt = (*_epy) * dt2;
            
            for ( _pxend = _px + nx; _px < _pxend; )
            {
                vx3 = *_vx3++;
                etaxdt = (*_epx++) * dt2;
                
                delta = ((vx3 - vx0 + (vx1 - vx2) * 27.0) * lax + 
                         ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay) * (*_mpx++);
            
                (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
                _px++;
                
                (*_py) = ((*_py) * (1.0 - etaydt) + delta) / (1.0 + etaydt);
                _py++;
                
                vx0 = vx1; vx1 = vx2; vx2 = vx3;
            }
            
            _px += px_a; _py += py_a; _mpx += mpx_a; _vx3 += vx_a; _vy0 += vy_a;
            _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; _epx -= nx; _epy += tsz;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts2d_24p(RDOM *dom, void *pars)
{
    int nx, ny, gxs, gys, px_a, py_a, mpx_a, vx_a, vy_a, iy, tsz, tid;
    register ireal * restrict _px, * restrict _py, * restrict _pxend,
                  * restrict _mpx, * restrict _vx3,
                  * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0;
    register ireal lax, lay, vx3, vx2, vx1, vx0, delta;
    RARR *s;

    s = dom->_s;
    
    nx = s[D_P0]._dims[0].n;
    ny = s[D_P0]._dims[1].n;
    if ( nx * ny == 0 ) return 0;
    
    gxs = s[D_P0]._dims[0].gs;
    gys = s[D_P0]._dims[1].gs;
    
    lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
    lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;
    
    #pragma omp parallel private(tsz,tid,iy,_px,_py,_pxend,_mpx,_vx3,_vy3,_vy2,_vy1,_vy0,vx3,vx2,vx1,vx0,delta)
    {
        #ifdef _OPENMP
        tsz = omp_get_num_threads();
        tid = omp_get_thread_num();
        #else
        tsz = 1;
        tid = 0;
        #endif
        
        #pragma omp single
        {            
            px_a  = tsz * s[D_P0 ]._dims[0].n0 - nx;
            py_a  = tsz * s[D_P1 ]._dims[0].n0 - nx;
            mpx_a = tsz * s[D_MP0]._dims[0].n0 - nx;
            vx_a  = tsz * s[D_V0 ]._dims[0].n0 - nx;
            vy_a  = tsz * s[D_V1 ]._dims[0].n0 - nx;
        }
        #pragma omp barrier
        
        _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0;
        _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
        _mpx = s[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs + tid) * s[D_MP0]._dims[0].n0;
        _vx3 = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0 + 1;
        _vy2 = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
        _vy3 = _vy2 + s[D_V1]._dims[0].n0;
        _vy1 = _vy2 - s[D_V1]._dims[0].n0;
        _vy0 = _vy1 - s[D_V1]._dims[0].n0;
        
        for ( iy = tid; iy < ny; iy += tsz )
        {
            vx2 = _vx3[-1]; vx1 = _vx3[-2]; vx0 = _vx3[-3];
            
            for ( _pxend = _px + nx; _px < _pxend; )
            {
                vx3 = *_vx3++;
                
                delta = ((vx3 - vx0 + (vx1 - vx2) * 27.0) * lax + 
                         ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay) * (*_mpx++);
                (*_px++) += delta;
                (*_py++) += delta;
                
                vx0 = vx1; vx1 = vx2; vx2 = vx3;
            }
            
            _px += px_a; _py += py_a; _mpx += mpx_a; _vx3 += vx_a; 
            _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts2d_24p0(RDOM *dom, void *pars)
{
    int nx, ny, gxs, gys, px_a, py_a, mpx_a, vx_a, vy_a, iy, tsz, tid;
    register ireal * restrict _px, * restrict _py, * restrict _pxend, * restrict _mpx,
                  * restrict _epx, * restrict _vx3,
                  * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0;
    register ireal lax, lay, dt2, vx3, vx2, vx1, vx0, delta, etaxdt;
    RARR *s;
    
    s = dom->_s;
    
    nx = s[D_P0]._dims[0].n;
    ny = s[D_P0]._dims[1].n;
    if ( nx * ny == 0 ) return 0;
    
    gxs = s[D_P0]._dims[0].gs;
    gys = s[D_P0]._dims[1].gs;
    
    lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
    lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;
    dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
    #pragma omp parallel private(tsz,tid,iy,_px,_py,_pxend,_mpx,_epx,_vx3,_vy3,_vy2,_vy1,_vy0,vx3,vx2,vx1,vx0,delta,etaxdt)
    {
        #ifdef _OPENMP
        tsz = omp_get_num_threads();
        tid = omp_get_thread_num();
        #else
        tsz = 1;
        tid = 0;
        #endif
        
        #pragma omp single
        {            
            px_a  = tsz * s[D_P0 ]._dims[0].n0 - nx;
            py_a  = tsz * s[D_P1 ]._dims[0].n0 - nx;
            mpx_a = tsz * s[D_MP0]._dims[0].n0 - nx;
            vx_a  = tsz * s[D_V0 ]._dims[0].n0 - nx;
            vy_a  = tsz * s[D_V1 ]._dims[0].n0 - nx;
        }
        #pragma omp barrier
        
        _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0;
        _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
        _mpx = s[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs + tid) * s[D_MP0]._dims[0].n0;
        _vx3 = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0 + 1;
        _vy2 = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
        _vy3 = _vy2 + s[D_V1]._dims[0].n0;
        _vy1 = _vy2 - s[D_V1]._dims[0].n0;
        _vy0 = _vy1 - s[D_V1]._dims[0].n0;
        _epx = s[D_EP0]._s + (gxs - s[D_EP0]._dims[0].gs);              /* 1D */
        
        for ( iy = tid; iy < ny; iy += tsz )
        {
            vx2 = _vx3[-1]; vx1 = _vx3[-2]; vx0 = _vx3[-3];
            
            for ( _pxend = _px + nx; _px < _pxend; )
            {
                vx3 = *_vx3++;
                etaxdt = (*_epx++) * dt2;
                
                delta = ((vx3 - vx0 + (vx1 - vx2) * 27.0) * lax + 
                         ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay) * (*_mpx++);
                (*_py++) += delta;            
                
                (*_px) = ((*_px) * (1.0 - etaxdt) + delta) / (1.0 + etaxdt);
                _px++;
                
                vx0 = vx1; vx1 = vx2; vx2 = vx3;
            }
            
            _px += px_a; _py += py_a; _mpx += mpx_a; _vx3 += vx_a;
            _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; _epx -= nx;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts2d_24p1(RDOM *dom, void *pars)
{
    int nx, ny, gxs, gys, px_a, py_a, mpx_a, vx_a, vy_a, iy, tsz, tid;
    register ireal * restrict _px, * restrict _py, * restrict _pxend, * restrict _mpx,
                  * restrict _epy, * restrict _vx3, 
                  * restrict _vy3, * restrict _vy2, * restrict _vy1, * restrict _vy0;
    register ireal lax, lay, dt2, vx3, vx2, vx1, vx0, delta, etaydt;
    RARR *s;
    
    s = dom->_s;
    
    nx = s[D_P0]._dims[0].n;
    ny = s[D_P0]._dims[1].n;
    if ( nx * ny == 0 ) return 0;
    
    gxs = s[D_P0]._dims[0].gs;
    gys = s[D_P0]._dims[1].gs;
    
    lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
    lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;
    dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
    #pragma omp parallel private(tsz,tid,iy,_px,_py,_pxend,_mpx,_epy,_vx3,_vy3,_vy2,_vy1,_vy0,vx3,vx2,vx1,vx0,delta,etaydt)
    {
        #ifdef _OPENMP
        tsz = omp_get_num_threads();
        tid = omp_get_thread_num();
        #else
        tsz = 1;
        tid = 0;
        #endif
        
        #pragma omp single
        {            
            px_a  = tsz * s[D_P0 ]._dims[0].n0 - nx;
            py_a  = tsz * s[D_P1 ]._dims[0].n0 - nx;
            mpx_a = tsz * s[D_MP0]._dims[0].n0 - nx;
            vx_a  = tsz * s[D_V0 ]._dims[0].n0 - nx;
            vy_a  = tsz * s[D_V1 ]._dims[0].n0 - nx;
        }
        #pragma omp barrier
        
        _px  = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0;
        _py  = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
        _mpx = s[D_MP0]._s + (gxs - s[D_MP0]._dims[0].gs) + (gys - s[D_MP0]._dims[1].gs + tid) * s[D_MP0]._dims[0].n0;
        _vx3 = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0 + 1;
        _vy2 = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
        _vy3 = _vy2 + s[D_V1]._dims[0].n0;
        _vy1 = _vy2 - s[D_V1]._dims[0].n0;
        _vy0 = _vy1 - s[D_V1]._dims[0].n0;
        _epy = s[D_EP1]._s + (gys - s[D_EP1]._dims[0].gs + tid);        /* 1D */
            
        for ( iy = tid; iy < ny; iy += tsz )
        {
            vx2 = _vx3[-1]; vx1 = _vx3[-2]; vx0 = _vx3[-3];
            etaydt = (*_epy) * dt2;
            
            for ( _pxend = _px + nx; _px < _pxend; )
            {
                vx3 = *_vx3++;
                
                delta = ((vx3 - vx0 + (vx1 - vx2) * 27.0) * lax + 
                         ((*_vy3++) - (*_vy0++) + ((*_vy1++) - (*_vy2++)) * 27.0) * lay) * (*_mpx++);
                (*_px++) += delta;            
                
                (*_py) = ((*_py) * (1.0 - etaydt) + delta) / (1.0 + etaydt);
                _py++;
                
                vx0 = vx1; vx1 = vx2; vx2 = vx3;
            }
            
            _px += px_a; _py += py_a; _mpx += mpx_a; _vx3 += vx_a; 
            _vy0 += vy_a; _vy1 += vy_a; _vy2 += vy_a; _vy3 += vy_a; _epy += tsz;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts2d_24v0(RDOM *dom, void *pars)
{
    int nx, ny, gxs, gys, vx_a, mvx_a, px_a, iy, tsz, tid;
    register ireal * restrict _vx, * restrict _vxend, * restrict _mvx,
                  * restrict _evx, * restrict _px3;
    register ireal lax, dt2, px3, px2, px1, px0, etaxdt;
    RARR *s;
    
    s = dom->_s;
    
    nx = s[D_V0]._dims[0].n;
    ny = s[D_V0]._dims[1].n;
    if ( nx * ny == 0 ) return 0;
    
    gxs = s[D_V0]._dims[0].gs;
    gys = s[D_V0]._dims[1].gs;
    
    lax = ((SGN_TS_PARS*)pars)->lam[0] * C24;
    dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
    #pragma omp parallel private(tsz,tid,iy,_vx,_vxend,_mvx,_evx,_px3,px3,px2,px1,px0,etaxdt)
    {
        #ifdef _OPENMP
        tsz = omp_get_num_threads();
        tid = omp_get_thread_num();
        #else
        tsz = 1;
        tid = 0;
        #endif
        
        #pragma omp single
        {
            vx_a  = tsz * s[D_V0 ]._dims[0].n0 - nx;
            mvx_a = tsz * s[D_MV0]._dims[0].n0 - nx;
            px_a  = tsz * s[D_P0 ]._dims[0].n0 - nx;
        }
        #pragma omp barrier
        
        _vx  = s[D_V0 ]._s + (gxs - s[D_V0 ]._dims[0].gs) + (gys - s[D_V0 ]._dims[1].gs + tid) * s[D_V0 ]._dims[0].n0;
        _mvx = s[D_MV0]._s + (gxs - s[D_MV0]._dims[0].gs) + (gys - s[D_MV0]._dims[1].gs + tid) * s[D_MV0]._dims[0].n0;
        _px3 = s[D_P0 ]._s + (gxs - s[D_P0 ]._dims[0].gs) + (gys - s[D_P0 ]._dims[1].gs + tid) * s[D_P0 ]._dims[0].n0 + 2;
        _evx = s[D_EV0]._s + (gxs - s[D_EV0]._dims[0].gs);              /* 1D */
        
        for ( iy = tid; iy < ny; iy += tsz )
        {
            px2 = _px3[-1]; px1 = _px3[-2]; px0 = _px3[-3];
            
            for ( _vxend = _vx + nx; _vx < _vxend; )
            {
                px3 = *_px3++;
                etaxdt = (*_evx++) * dt2;
                
                (*_vx) = ((*_vx) * (1.0 - etaxdt) + (px3 - px0 + (px1 - px2) * 27.0) * lax * (*_mvx++)) / (1.0 + etaxdt);
                _vx++;
                
                px0 = px1; px1 = px2; px2 = px3;
            }
            
            _vx += vx_a; _mvx += mvx_a; _px3 += px_a; _evx -= nx;
        }
    } /* omp parallel */

    return 0;
}
/*----------------------------------------------------------------------------*/

int sgn_ts2d_24v1(RDOM *dom, void *pars)
{
    int nx, ny, gxs, gys, vy_a, mvy_a, py_a, iy, tsz, tid;
    register ireal * restrict _vy, * restrict _vyend, * restrict _mvy, * restrict _evy,
                  * restrict _py3, * restrict _py2, * restrict _py1, * restrict _py0;
    register ireal lay, dt2, etaydt;
    RARR *s;
    
    s = dom->_s;
    
    nx = s[D_V1]._dims[0].n;
    ny = s[D_V1]._dims[1].n;
    if ( nx * ny == 0 ) return 0;
    
    gxs = s[D_V1]._dims[0].gs;
    gys = s[D_V1]._dims[1].gs;
    
    lay = ((SGN_TS_PARS*)pars)->lam[1] * C24;
    dt2 = ((SGN_TS_PARS*)pars)->dt / 2.0;
    
    #pragma omp parallel private(tsz,tid,iy,_vy,_vyend,_mvy,_evy,_py3,_py2,_py1,_py0,etaydt)
    {
        #ifdef _OPENMP
        tsz = omp_get_num_threads();
        tid = omp_get_thread_num();
        #else
        tsz = 1;
        tid = 0;
        #endif
        
        #pragma omp single
        {
            vy_a  = tsz * s[D_V1 ]._dims[0].n0 - nx;
            mvy_a = tsz * s[D_MV1]._dims[0].n0 - nx;
            py_a  = tsz * s[D_P1 ]._dims[0].n0 - nx;
        }
        #pragma omp barrier
        
        _vy  = s[D_V1 ]._s + (gxs - s[D_V1 ]._dims[0].gs) + (gys - s[D_V1 ]._dims[1].gs + tid) * s[D_V1 ]._dims[0].n0;
        _mvy = s[D_MV1]._s + (gxs - s[D_MV1]._dims[0].gs) + (gys - s[D_MV1]._dims[1].gs + tid) * s[D_MV1]._dims[0].n0;
        _py1 = s[D_P1 ]._s + (gxs - s[D_P1 ]._dims[0].gs) + (gys - s[D_P1 ]._dims[1].gs + tid) * s[D_P1 ]._dims[0].n0;
        _py2 = _py1 + s[D_P1]._dims[0].n0;
        _py3 = _py2 + s[D_P1]._dims[0].n0;
        _py0 = _py1 - s[D_P1]._dims[0].n0;
        _evy = s[D_EV1]._s + (gys - s[D_EV1]._dims[0].gs + tid);        /* 1D */
        
        for ( iy = tid; iy < ny; iy += tsz )
        {            
            etaydt = (*_evy) * dt2;
            
            for ( _vyend = _vy + nx; _vy < _vyend; )
            {
                (*_vy) = ((*_vy) * (1.0 - etaydt) + ((*_py3++) - (*_py0++) + 
                          ((*_py1++) - (*_py2++)) * 27.0) * lay * (*_mvy++)) / (1.0 + etaydt);
                _vy++;
            }
            
            _vy += vy_a; _mvy += mvy_a; _py0 += py_a; _py1 += py_a; _py2 += py_a; _py3 += py_a; _evy += tsz;
        }
    } /* omp parallel */

    return 0;
}
