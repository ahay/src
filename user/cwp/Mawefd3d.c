/* 3D acoustic time-domain FD modeling 
   2Nth order in space, 2nd order in time
   with optimized FD scheme option and hybrid one-way ABC option 
   adj flag indicates backwards source injection, not exact adjoint propagator
*/
/*
  Copyright (C) 2014 Colorado School of Mines
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
#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
#endif
#include "fdutil.h"

#ifndef SF_HAS_SSE
#undef __SSE__
#endif
#ifndef SF_HAS_AVX
#undef __AVX__
#endif

#ifdef sun
#define restrict
#endif

#if defined __SSE__ || defined __AVX__
#include <immintrin.h>
#endif

static float* compute_fdcoef(int nop, float dz, float dx, float dy,
                             bool is_optimized, const int d_order);
static int factorial(int n);
static float* normal_fdcoef(int nop, const int d_order);
static float* optimal_fdcoef(int nop, const int d_order);

static void step_forward(float*** restrict u0, float*** restrict u1,
                         float*** restrict vel, float*** restrict rho,
                         float* restrict fdcoef_d2, float* restrict fdcoef_d1,
                         int nop, int nzpad, int nxpad, int nypad);
static float* damp_make(int ntransit);
static void apply_abc(float*** uu2, float*** uu1, int nz, int nx, int ny, int nbd,
		      abcone3d abc, int nop, float* damp);

int
main(int argc, char** argv)
{

    bool verb, fsrf, snap, expl, dabc, cden, adj;
    bool optfd, hybrid, sinc; 
    int jsnap, jdata;

    /* I/O files */
    sf_file file_wav=NULL; /* wavelet */
    sf_file file_vel=NULL; /* velocity */
    sf_file file_den=NULL; /* density */
    sf_file file_wfl=NULL; /* wavefield */
    sf_file file_dat=NULL; /* data */
    sf_file file_src=NULL; /* sources */
    sf_file file_rec=NULL; /* receivers */

    /* cube axes */
    sf_axis at = NULL, az = NULL, ax = NULL, ay = NULL;
    sf_axis as = NULL, ar = NULL;

    int nbd;  /* ABC boundary size */
    int fdorder;  /* finite difference spatial accuracy order */
    int nzpad,nxpad,nypad; /* boundary padded model size */
    int ix,iy,it,is,nx,ny,nz,nt,ns,nr;
    float dx,dy,dz,dt,dt2;
    float* damp=NULL; /* damping profile for hybrid bc */
    float* ws;  /* wavelet */
    float*** vel=NULL;  /* velocity */
    float*** rho=NULL; /* density */
    float*** u0=NULL;  /* wavefield array u@t-1 (u@t+1) */
    float*** u1=NULL;  /* wavefield array u@t */
    float* u_dat=NULL; /* output data */
    float*** ptr_tmp=NULL;   
    pt3d* src3d=NULL;  /* source position */
    pt3d* rec3d=NULL;  /*receiver position*/
    scoef3d cssinc = NULL, crsinc = NULL; 
    lint3d cslint = NULL, crlint = NULL;

    /* FDM structure */
    fdm3d fdm = NULL;
    abcone3d abc = NULL;
    sponge spo = NULL;

    int nbell;

    float* fdcoef_d2;
    float* fdcoef_d1;

    sf_axis acz = NULL, acx = NULL, acy = NULL;
    int nqz, nqx, nqy;
    float oqz, oqx, oqy, dqz, dqx, dqy;

    float** oslice = NULL; /* output 3D wavefield slice-by-slice */
    float*** tmp_array;

    double wall_clock_time_s, wall_clock_time_e;

    const int SECOND_DERIV = 2;
    const int FIRST_DERIV = 1;

    int nop;

#if defined _OPENMP && _DEBUG
    double tic;
    double toc;
#endif

    /* init RSF */
    sf_init(argc,argv);

#ifdef _OPENMP
    omp_init();
    wall_clock_time_s = omp_get_wtime();
#else
    wall_clock_time_s = (double) clock() / CLOCKS_PER_SEC;
#endif

    if (!sf_getbool("verb",&verb))  verb=false; /* Verbosity flag */
    if (!sf_getbool("snap",&snap))  snap=false; /* Wavefield snapshots flag */
    if (!sf_getbool("expl",&expl))  expl=false; /* Multiple sources, one wvlt*/
    if (!sf_getbool("dabc",&dabc))  dabc=false; /* Absorbing BC */
    if (!sf_getbool("cden",&cden))  cden=false; /* Constant density */
    if (!sf_getbool("adj",&adj))    adj=false; /* adjoint flag */

    if (!sf_getbool("free",&fsrf) && !sf_getbool("fsrf",&fsrf)) fsrf=false; /* Free surface flag */

    if (!sf_getint("nbell",&nbell)) nbell=5; /* gaussian for source injection */

    if (!sf_getbool("optfd",&optfd))  optfd=false; /* optimized FD coefficients flag */
    if (!sf_getint("fdorder",&fdorder))  fdorder=4; /* spatial FD order */
    if (!sf_getbool("hybridbc",&hybrid))  hybrid=false;  /* hybrid Absorbing BC */
    if (!sf_getbool("sinc",&sinc)) sinc=false; /* sinc source injection */
  
    /* Initialize variables */
    file_wav = sf_input("in"); /* wavelet */
    file_vel = sf_input("vel"); /* velocity */ 
    file_src = sf_input("sou"); /* sources */
    file_rec = sf_input("rec"); /* receivers */
    file_dat = sf_output("out"); /* data */

    if (snap)  file_wfl = sf_output("wfl"); /* wavefield */
    if (!cden) {
	if (sf_getstring("cden")) {
	    file_den = sf_input ("den"); /* density */
	} else {
	    cden = true;
	    if (verb) sf_warning("No density file provided, running with constant density");
	}
    }
  
    at = sf_iaxa(file_wav,2); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    az = sf_iaxa(file_vel,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(file_vel,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */
    ay = sf_iaxa(file_vel,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /* space */

    as = sf_iaxa(file_src,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
    ar = sf_iaxa(file_rec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */

    nt = sf_n(at); dt = sf_d(at); 
    nz = sf_n(az); dz = sf_d(az); 
    nx = sf_n(ax); dx = sf_d(ax); 
    ny = sf_n(ay); dy = sf_d(ay); 
    ns = sf_n(as);
    nr = sf_n(ar);

    /* other execution parameters */
    if (snap) {
	if (!sf_getint("jsnap",&jsnap))  jsnap=nt;
	/* # of t steps at which to save wavefield */
    }
    if (!sf_getint("jdata",&jdata)) jdata=1;
    /* # of t steps at which to save receiver data */

    /* setup output data header */
    sf_oaxa(file_dat,ar,1);
    sf_setn(at,(nt-1)/jdata+1);
    sf_setd(at,dt*jdata);
    sf_oaxa(file_dat,at,2);

    /* wavefield cut params */
    /* setup output wavefield header */
    if (snap) {
	if (!sf_getint  ("nqz",&nqz)) nqz=sf_n(az); /* Saved wfld window nz */
	if (!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax); /* Saved wfld window nx */
	if (!sf_getint  ("nqy",&nqy)) nqy=sf_n(ay); /* Saved wfld window ny */

	if (!sf_getfloat("oqz",&oqz)) oqz=sf_o(az); /* Saved wfld window oz */
	if (!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax); /* Saved wfld window ox */
	if (!sf_getfloat("oqy",&oqy)) oqy=sf_o(ay); /* Saved wfld window oy */

	if (!sf_getfloat("dqz",&dqz)) dqz=sf_d(az); /* Saved wfld window dz */
	if (!sf_getfloat("dqx",&dqx)) dqx=sf_d(ax); /* Saved wfld window dx */
	if (!sf_getfloat("dqy",&dqy)) dqy=sf_d(ay); /* Saved wfld window dy */
    
	acz = sf_maxa(nqz,oqz,dqz); if (verb) sf_raxa(acz);
	acx = sf_maxa(nqx,oqx,dqx); if (verb) sf_raxa(acx);
	acy = sf_maxa(nqy,oqy,dqy); if (verb) sf_raxa(acy);
	/* check if the imaging window fits in the wavefield domain */
	sf_setn(at,(nt-1)/jsnap+1);
	sf_setd(at,dt*jsnap);
	if (verb) sf_raxa(at);
    
	sf_oaxa(file_wfl,acz,1);
	sf_oaxa(file_wfl,acx,2);
	sf_oaxa(file_wfl,acy,3);
	sf_oaxa(file_wfl,at,4);
    }

    /* 2-2N finite difference coefficient */
    nop = fdorder/2; /* fd half-length stencil */
    if (!sf_getint("nb",&nbd) || nbd<nop)  nbd=nop;
    if (dabc && hybrid && nbd<=nop) nbd = 2*nop;

    /* expand domain for FD operators and ABC */
    fdm = fdutil3d_init(verb,fsrf,az,ax,ay,nbd,1);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if (verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if (verb) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); if (verb) sf_raxa(ay);

    /* Precompute coefficients */
    dt2 = dt*dt;
    nzpad = nz+2*nbd;  nxpad = nx+2*nbd;  nypad = ny+2*nbd;

    fdcoef_d2 = compute_fdcoef(nop,dz,dx,dy,optfd,SECOND_DERIV);
    fdcoef_d1 = compute_fdcoef(nop,dz,dx,dy,optfd,FIRST_DERIV);

    /* Allocate memories */
    if (expl) ws = sf_floatalloc(1);
    else      ws = sf_floatalloc(ns);
    vel = sf_floatalloc3(nzpad,nxpad,nypad);
    if (!cden) rho = sf_floatalloc3(nzpad,nxpad,nypad);
    u_dat = sf_floatalloc(nr);
    src3d = pt3dalloc1(ns);
    rec3d = pt3dalloc1(nr);
    if (snap) oslice = sf_floatalloc2(sf_n(acz),sf_n(acx));

    /* source and receiver position */
    pt3dread1(file_src,src3d,ns,3);  /* read format: (x,y,z) */
    if (sinc) cssinc = sinc3d_make(ns,src3d,fdm);
    else      cslint = lint3d_make(ns,src3d,fdm);

    pt3dread1(file_rec,rec3d,nr,3);  /* read format: (x,y,z) */
    if (sinc) crsinc = sinc3d_make(nr,rec3d,fdm);
    else      crlint = lint3d_make(nr,rec3d,fdm);

    if (!sinc) fdbell3d_init(nbell);

    /* temperary array */
    tmp_array = sf_floatalloc3(nz,nx,ny);

    /* read velocity and pad */
    sf_floatread(tmp_array[0][0],nz*nx*ny,file_vel);
    expand3d(tmp_array,vel,fdm);
    /* read density and pad */
    if (!cden) {
	sf_floatread(tmp_array[0][0],nz*nx*ny,file_den);
	expand3d(tmp_array,rho,fdm);
    }

    free(**tmp_array);  free(*tmp_array);  free(tmp_array);

    /* A1 one-way ABC implicit scheme coefficients  */
    if (dabc) {
	abc = abcone3d_make(nbd,dt,vel,fsrf,fdm);
	if (hybrid)
	    damp = damp_make(nbd-nop); /* compute damping profiles for hybrid bc */
	else
	    spo = sponge_make(fdm->nb);
    }

    /* allocate memory for wavefield variables */
    u0 = sf_floatalloc3(nzpad,nxpad,nypad);
    u1 = sf_floatalloc3(nzpad,nxpad,nypad);

    /* initialize variables */
    memset(u0[0][0],0,sizeof(float)*nzpad*nxpad*nypad);
    memset(u1[0][0],0,sizeof(float)*nzpad*nxpad*nypad);
    memset(u_dat,0,sizeof(float)*nr);

    /* v = (v*dt)^2 */
    for (ix=0;ix<nzpad*nxpad*nypad;ix++)
	*(vel[0][0]+ix) *= *(vel[0][0]+ix)*dt2;
    if (fsrf && !hybrid) {
	for (iy=0; iy<nypad; iy++)
	    for (ix=0; ix<nxpad; ix++)
		memset(vel[iy][ix],0,sizeof(float)*fdm->nb);
    }

    for (it=0; it<nt; it++) {
	if (verb)  sf_warning("it=%d;",it+1);
#if defined _OPENMP && _DEBUG
	tic=omp_get_wtime();
#endif
    
	step_forward(u0,u1,vel,rho,fdcoef_d2,fdcoef_d1,nop,nzpad,nxpad,nypad);
    
	if (adj) { /* backward inject source wavelet */
	    if (expl) {
		sf_seek(file_wav,(off_t)(nt-it-1)*sizeof(float),SEEK_SET);
		sf_floatread(ws,1,file_wav);
		ws[0] *= dt2;
		if (sinc) sinc3d_inject1(u0,ws[0],cssinc);
		else      lint3d_inject1(u0,ws[0],cslint);
	    } else { 
		sf_seek(file_wav,(off_t)(nt-it-1)*ns*sizeof(float),SEEK_SET);
		sf_floatread(ws,ns,file_wav);
		for (is=0; is<ns; is++) ws[is] *= dt2;
		if (sinc) sinc3d_inject(u0,ws,cssinc);
		else      lint3d_inject(u0,ws,cslint);
	    }
	} else { /* forward inject source wavelet */
	    if (expl) {
		sf_floatread(ws,1,file_wav);
		ws[0] *= dt2;
		if (sinc) sinc3d_inject1(u0,ws[0],cssinc);
		else      lint3d_inject1(u0,ws[0],cslint);
	    } else {
		sf_floatread(ws,ns,file_wav);
		for (is=0; is<ns; is++) ws[is] *= dt2;
		if (sinc) sinc3d_inject(u0,ws,cssinc);
		else      lint3d_inject(u0,ws,cslint);
	    }
	}

	/* apply abc */
	if (dabc) {
	    if (hybrid) apply_abc(u0,u1,nz,nx,ny,nbd,abc,nop,damp);
	    else {
		abcone3d_apply(u0,u1,nop,abc,fdm);
		sponge3d_apply(u0,spo,fdm);
		sponge3d_apply(u1,spo,fdm);
	    }
	}

	/* loop over pointers */
	ptr_tmp = u0;  u0 = u1;  u1 = ptr_tmp;

	/* extract snapshot */
	if (snap && it%jsnap==0) {
	    int fy = (floor)((sf_o(acy)-fdm->oypad)/fdm->dy);
	    int jy = floor(sf_d(acy)/fdm->dy);
	    float **ptr_slice;
	    for (iy=0; iy<sf_n(acy); iy++) {
		ptr_slice = u0[fy+iy*jy];
		cut3d_slice(ptr_slice,oslice,fdm,acz,acx);
		sf_floatwrite(oslice[0],sf_n(acz)*sf_n(acx),file_wfl);
	    }
	}

	/* extract receiver data */
	if (sinc) sinc3d_extract(u0,u_dat,crsinc);
	else      lint3d_extract(u0,u_dat,crlint);

	sf_floatwrite(u_dat,nr,file_dat);

#if defined _OPENMP && _DEBUG
	toc=omp_get_wtime(); 
	fprintf(stderr,"%5.2gs",(float)(toc-tic));
#endif
    }
#ifdef _OPENMP
    wall_clock_time_e = omp_get_wtime();
#else
    wall_clock_time_e = (double) clock() / CLOCKS_PER_SEC;
#endif
    if (verb)
	fprintf(stderr,"\nElapsed time: %lf s\n",wall_clock_time_e-wall_clock_time_s);

    free(**u0); free(*u0); free(u0);
    free(**u1); free(*u1); free(u1);
    free(**vel); free(*vel); free(vel);
    free(u_dat);
    free(ws);
    free(fdcoef_d2); free(fdcoef_d1);
    if (snap) { free(*oslice); free(oslice); }
    if(!cden) { free(**rho); free(*rho); free(rho); }
    if (hybrid) free(damp);
    free(src3d); free(rec3d);

    return 0;
}


static void 
step_forward(float*** restrict u0, float*** restrict u1,
             float*** restrict vel, float*** restrict rho,
             float* restrict fdcoef_d2, float* restrict fdcoef_d1,
             int nop, int nzpad, int nxpad, int nypad)
{
    int ix, iy, iz, iop;
    float c0 = fdcoef_d2[0];
    float *cz = &fdcoef_d2[0];
    float *cx = &fdcoef_d2[nop];
    float *cy = &fdcoef_d2[nop+nop];
    float *bz, *bx, *by;
    float drho_dot_du;
    float du_z = 0.f, du_x = 0.f, du_y = 0.f;
    float drho_z = 0.f, drho_x = 0.f, drho_y = 0.f;
    float lap;

    if (rho != NULL) { /* variable density */
	bz = &fdcoef_d1[0];
	bx = &fdcoef_d1[nop];
	by = &fdcoef_d1[nop+nop];
    }

#ifdef __SSE__
#define vType __m128
#define _prefix_loadu_ps _mm_loadu_ps
#define _prefix_storeu_ps _mm_storeu_ps
#define _prefix_add_ps _mm_add_ps
#define _prefix_sub_ps _mm_sub_ps
#define _prefix_mul_ps _mm_mul_ps
#define _prefix_div_ps _mm_div_ps
#define _prefix_set1_ps _mm_set1_ps
#elif defined __AVX__
#define vType __m256
#define _prefix_loadu_ps _mm256_loadu_ps
#define _prefix_storeu_ps _mm256_storeu_ps
#define _prefix_add_ps _mm256_add_ps
#define _prefix_sub_ps _mm256_sub_ps
#define _prefix_mul_ps _mm256_mul_ps
#define _prefix_div_ps _mm256_div_ps
#define _prefix_set1_ps _mm256_set1_ps
#endif

#ifdef _OPENMP
#pragma omp parallel for				\
    schedule(static,1)					\
    shared(nxpad,nypad,nzpad,u0,u1,vel,c0,cx,cy,cz) \
    private(ix,iy,iz,iop,drho_dot_du,du_z,du_x,du_y,drho_z,drho_x,drho_y,lap)
#endif
#if defined (__SSE__) || defined (__AVX__)
    for (iy=nop; iy<nypad-nop; iy++) {
	for (ix=nop; ix<nxpad-nop; ix++) {
#ifdef __SSE__
	    for (iz=nop; iz<nzpad-nop; iz+=4) {
#elif defined __AVX__
		for (iz=nop; iz<nzpad-nop; iz+=8) {
#endif
		    /* load u0 u1 vel from memory to register */
		    vType vec_u0 = _prefix_loadu_ps(&u0[iy][ix][iz]);
		    vType vec_u1 = _prefix_loadu_ps(&u1[iy][ix][iz]);
		    vec_u0 = _prefix_sub_ps(_prefix_mul_ps(_prefix_set1_ps(2.f),vec_u1),vec_u0);

		    vType vec_lap = _prefix_mul_ps(vec_u1,_prefix_set1_ps(c0));
		    for (iop=1; iop<=nop; iop++) {
			/* z axis derivative <1st axis> */
			vec_lap = _prefix_add_ps(vec_lap,
						 _prefix_mul_ps(
						     _prefix_set1_ps(cz[iop]),
						     _prefix_add_ps(
							 _prefix_loadu_ps(&u1[iy][ix][iz+iop]),
							 _prefix_loadu_ps(&u1[iy][ix][iz-iop])
							 )
						     )
			    );
			/* x axis derivative <2nd axis> */
			vec_lap = _prefix_add_ps(vec_lap,
						 _prefix_mul_ps(
						     _prefix_set1_ps(cx[iop]),
						     _prefix_add_ps(
							 _prefix_loadu_ps(&u1[iy][ix+iop][iz]),
							 _prefix_loadu_ps(&u1[iy][ix-iop][iz])
							 )
						     )
			    );
			/* y axis derivative <3rd axis> */
			vec_lap = _prefix_add_ps(vec_lap,
						 _prefix_mul_ps(
						     _prefix_set1_ps(cy[iop]),
						     _prefix_add_ps(
							 _prefix_loadu_ps(&u1[iy+iop][ix][iz]),
							 _prefix_loadu_ps(&u1[iy-iop][ix][iz])
							 )
						     )
			    );
		    }

		    if (rho != NULL) {
			vType vec_du_z = _prefix_set1_ps(0.f);
			vType vec_du_x = _prefix_set1_ps(0.f);
			vType vec_du_y = _prefix_set1_ps(0.f);
			vType vec_drho_z = _prefix_set1_ps(0.f);
			vType vec_drho_x = _prefix_set1_ps(0.f);
			vType vec_drho_y = _prefix_set1_ps(0.f);
			for (iop=1; iop<=nop; iop++) {
			    vec_du_z = _prefix_add_ps(vec_du_z,
						      _prefix_mul_ps(
							  _prefix_set1_ps(bz[iop]),
							  _prefix_sub_ps(
							      _prefix_loadu_ps(&u1[iy][ix][iz+iop]),
							      _prefix_loadu_ps(&u1[iy][ix][iz-iop])
							      )
							  )
				);
			    vec_drho_z = _prefix_add_ps(vec_drho_z,
							_prefix_mul_ps(
							    _prefix_set1_ps(bz[iop]),
							    _prefix_sub_ps(
								_prefix_loadu_ps(&rho[iy][ix][iz+iop]),
								_prefix_loadu_ps(&rho[iy][ix][iz-iop])
								)
							    )
				);
			    vec_du_x = _prefix_add_ps(vec_du_x,
						      _prefix_mul_ps(
							  _prefix_set1_ps(bx[iop]),
							  _prefix_sub_ps(
							      _prefix_loadu_ps(&u1[iy][ix+iop][iz]),
							      _prefix_loadu_ps(&u1[iy][ix-iop][iz])
							      )
							  )
				);
			    vec_drho_x = _prefix_add_ps(vec_drho_x,
							_prefix_mul_ps(
							    _prefix_set1_ps(bx[iop]),
							    _prefix_sub_ps(
								_prefix_loadu_ps(&rho[iy][ix+iop][iz]),
								_prefix_loadu_ps(&rho[iy][ix-iop][iz])
								)
							    )
				);
			    vec_du_y = _prefix_add_ps(vec_du_y,
						      _prefix_mul_ps(
							  _prefix_set1_ps(by[iop]),
							  _prefix_sub_ps(
							      _prefix_loadu_ps(&u1[iy+iop][ix][iz]),
							      _prefix_loadu_ps(&u1[iy-iop][ix][iz])
							      )
							  )
				);
			    vec_drho_y = _prefix_add_ps(vec_drho_y,
							_prefix_mul_ps(
							    _prefix_set1_ps(by[iop]),
							    _prefix_sub_ps(
								_prefix_loadu_ps(&rho[iy+iop][ix][iz]),
								_prefix_loadu_ps(&rho[iy-iop][ix][iz])
								)
							    )
				);
			}
			vec_lap = _prefix_sub_ps(vec_lap,
						 _prefix_div_ps(
						     _prefix_add_ps(
							 _prefix_mul_ps(vec_du_z,vec_drho_z),
							 _prefix_add_ps(
							     _prefix_mul_ps(vec_du_x,vec_drho_x),
							     _prefix_mul_ps(vec_du_y,vec_drho_y)
							     )
							 ),
						     _prefix_loadu_ps(&rho[iy][ix][iz])
						     )
			    );
		    }
		    vec_u0 = _prefix_add_ps(vec_u0,_prefix_mul_ps(_prefix_loadu_ps(&vel[iy][ix][iz]),vec_lap));
		    _prefix_storeu_ps(&u0[iy][ix][iz],vec_u0);
		}
	    }
	}
#undef vType
#undef _prefix_loadu_ps
#undef _prefix_storeu_ps
#undef _prefix_add_ps
#undef _prefix_sub_ps
#undef _prefix_mul_ps
#undef _prefix_div_ps
#undef _prefix_set1_ps

#else
	for (iy=nop; iy<nypad-nop; iy++) {
	    for (ix=nop; ix<nxpad-nop; ix++) {
		for (iz=nop; iz<nzpad-nop; iz++) {
		    lap = u1[iy][ix][iz]*c0;
		    for (iop=1; iop<=nop; iop++) {
			lap += (u1[iy][ix][iz-iop] + u1[iy][ix][iz+iop]) * cz[iop]
			    + (u1[iy][ix-iop][iz] + u1[iy][ix+iop][iz]) * cx[iop]
			    + (u1[iy-iop][ix][iz] + u1[iy+iop][ix][iz]) * cy[iop];
		    }
		    if (rho != NULL) { /* variable density term */
        du_z = du_x = du_y = drho_z = drho_x = drho_y = 0.f;
			for (iop=1; iop<=nop; iop++) {
			    du_z += (u1[iy][ix][iz+iop] - u1[iy][ix][iz-iop]) * bz[iop];
			    du_x += (u1[iy][ix+iop][iz] - u1[iy][ix-iop][iz]) * bx[iop];
			    du_y += (u1[iy+iop][ix][iz] - u1[iy-iop][ix][iz]) * by[iop];
			    drho_z += (rho[iy][ix][iz+iop] - rho[iy][ix][iz-iop]) * bz[iop];
			    drho_x += (rho[iy][ix+iop][iz] - rho[iy][ix-iop][iz]) * bx[iop];
			    drho_y += (rho[iy+iop][ix][iz] - rho[iy-iop][ix][iz]) * by[iop];
			}
			drho_dot_du = (du_z*drho_z + du_x*drho_x + du_y*drho_y)/rho[iy][ix][iz];
			lap -= drho_dot_du;
		    }
		    u0[iy][ix][iz] = 2.*u1[iy][ix][iz] - u0[iy][ix][iz] + vel[iy][ix][iz]*lap;
		}
	    }
	}
#endif
	return;
    }

static float*
    damp_make(int ntransit)
{
    int ib;
    float* damp = NULL;
    float sb, fb;
    if (ntransit>0) damp = sf_floatalloc(ntransit);
    sb = 4.0*ntransit;
    for(ib=0; ib<ntransit; ib++) {
	fb = ib/(sqrt(2.0)*sb);
	damp[ib] = exp(-fb*fb);
    }
    return damp;
}

static float*
    compute_fdcoef(int nop, float dz, float dx, float dy, bool is_optimized, const int d_order)
    /*
      Optimized fd coeffifients from:
      1. Yang Liu. "Globally optimal finite-difference schemes based on least squares" Geophysics 78.4 (2013)
      2. Zhang, Jin-Hai, and Zhen-Xing Yao. "Optimized explicit finite-difference schemes for spatial derivatives using maximum norm." Journal of Computational Physics 250 (2013)

      Conventional fd coefficients formula from:
      3. Chunlei Chu, Paul Stoffa. "Determination of finite-difference weights using scaled binomial windows" Geophysics 77.3 (2012)
    */
{
    int idim, ii;
#define ndim 3
    float d2[ndim];
    float d1[ndim];
    float *d2_fdcoef;
    float *ccc;
    float *d1_fdcoef;
    float *bbb;

    d2[0] = dz*dz; d2[1] = dx*dx; d2[2] = dy*dy;
    d1[0] = dz; d1[1] = dx; d1[2] = dy;

    if (d_order == 2) {
	d2_fdcoef = sf_floatalloc(ndim*nop+1);
	if (is_optimized) ccc= optimal_fdcoef(nop,2);
	else ccc = normal_fdcoef(nop,d_order);
	
	d2_fdcoef[0] = 0.f;
	for (idim=0; idim<ndim; idim++) {
	    for (ii=1; ii<=nop; ii++) {
		d2_fdcoef[idim*nop+ii] = ccc[ii]/d2[idim];
		d2_fdcoef[0] += d2_fdcoef[idim*nop+ii];
	    }
	}
	d2_fdcoef[0] *= - 2.0f;
	return d2_fdcoef;
    } else {
	d1_fdcoef = sf_floatalloc(ndim*nop+1);
	
	if (is_optimized && nop<=6) bbb = optimal_fdcoef(nop,1);
	else bbb = normal_fdcoef(nop,d_order);
	d1_fdcoef[0] = 0.0f;
	for (idim=0; idim<ndim; idim++) {
	    for (ii=1; ii<=nop; ii++) {
		d1_fdcoef[idim*nop+ii] = bbb[ii]/d1[idim];
	    }
	}
	return d1_fdcoef;
    }
}

static int
    factorial(int n)
{
    int i;
    int result = 1;
    for (i=1; i<=n; ++i)
	result *= i;
    return result;
}

static float*
    normal_fdcoef(int nop, const int d_order)
{
    int n;
    float *cc = calloc(nop+1,sizeof(float));
    float *bb = calloc(nop+1,sizeof(float));
    int halfN_fact = factorial(nop);
    halfN_fact *= halfN_fact;
    for (n=1; n<=nop; n++) {
	cc[n] = - 2.f / (n*n) * cos(n*SF_PI) * halfN_fact/factorial(nop+n)/factorial(nop-n); 
	bb[n] = cc[n]*n/2.f;
	}
    if (d_order == 1) return bb;
    else return cc;
}

static float*
    optimal_fdcoef(int nop, const int d_order)
{
    float *cc;
    float *bb;

    if (d_order == 2)  {
	float *opt_c[11];
	float opt_c1[2] = {0.f, 1.f}; 
	float opt_c2[3] = {0.f, 1.369074f, - 0.09266816}; 
	float opt_c3[4] = {0.f, 1.573661f, - 0.1820268f, 0.01728053f}; 
	float opt_c4[5] = {0.f, 1.700010f, - 0.2554615f, 0.04445392f, - 0.004946851f}; 
	float opt_c5[6] = {0.f, 1.782836f, - 0.3124513f, 0.07379487f, - 0.01532122f,
			   0.001954439f}; 
	float opt_c6[7] = {0.f, 1.837023f, - 0.3538895f, 0.09978343f, - 0.02815486f,
			   0.006556587f, - 0.0009405699f}; 
	float opt_c7[8] = {0.f, 1.874503f, - 0.3845794f, 0.1215162f, - 0.04121749f,
			   0.01295522f, - 0.003313813f, 0.0005310053f}; 
	float opt_c8[9] = {0.f, 1.901160f, - 0.4074304f, 0.1390909f, - 0.05318775f,
			   0.02004823f, - 0.006828249f, 0.001895771f, - 0.0003369052f}; 
	float opt_c9[10] = {0.f, 1.919909f, - 0.4240446f, 0.1526043f, - 0.06322328f,
			    0.02676005f, - 0.01080739f, 0.003907747f, - 0.001158024f,
				0.0002240247f}; 
	float opt_c10[11] = {0.f, 1.918204f, - 0.4225858f, 0.1514992f, - 0.06249474f,
			     0.02637196f, - 0.01066631f, 0.003915625f, - 0.001219872f,
			     0.0002863976f, - 0.00003744830f}; 
	opt_c[1] = opt_c1;  opt_c[2] = opt_c2;  opt_c[3] = opt_c3;  opt_c[4] = opt_c4;
	opt_c[5] = opt_c5;  opt_c[6] = opt_c6;  opt_c[7] = opt_c7;  opt_c[8] = opt_c8;
	opt_c[9] = opt_c9;  opt_c[10] = opt_c10;
	cc = sf_floatalloc(nop+1);
	memcpy(cc,opt_c[nop],sizeof(float)*(nop+1));
	return cc;
    } else {
	float *opt_b[7];
	float opt_b1[2] = {0.0f, 0.5f}; 
	float opt_b2[3] = {0.0f, 0.67880327, - 0.08962729}; 
	float opt_b3[4] = {0.0f, 0.77793115, - 0.17388691, 0.02338713}; 
	float opt_b4[5] = {0.0f, 0.84149635, - 0.24532989, 0.06081891, - 0.00839807}; 
	float opt_b5[6] = {0.0f, 0.88414717, - 0.30233648, 0.10275057, - 0.02681517,
			   0.00398089}; 
	float opt_b6[7] = {0.0f, 0.91067892, - 0.34187892, 0.13833962, - 0.04880710,
			   0.01302148, - 0.00199047}; 
	
	opt_b[1] = opt_b1;  opt_b[2] = opt_b2;  opt_b[3] = opt_b3;
	opt_b[4] = opt_b4;  opt_b[5] = opt_b5;  opt_b[6] = opt_b6;
	bb = sf_floatalloc(nop+1);
	memcpy(bb,opt_b[nop],sizeof(float)*(nop+1));
	return bb;
    }
}

static void
apply_abc(float*** restrict uu2, float*** restrict uu1, int nz, int nx, int ny, int nbd,
	  abcone3d abc, int nop, float* restrict damp)
{
    int ix, iy, iz, ib;
    float damp_ib;
    float uu2_bc;
    int nxpad = nx + 2*nbd;
    int nypad = ny + 2*nbd;
    int nzpad = nz + 2*nbd;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(ix,iy,iz,ib,damp_ib,uu2_bc)
#endif
    for (iy=0; iy<nypad; iy++) {
	for (ix=0; ix<nxpad; ix++) {
	    for (ib=nbd-nop; ib<nbd; ib++) {
		iz = nbd-ib-1;
		if (abc->free) {
		    uu2[iy][ix][iz] = 0.0f;
		} else {
		    uu2[iy][ix][iz] = uu1[iy][ix][iz+1] 
			+ (uu1[iy][ix][iz] - uu2[iy][ix][iz+1])*abc->bzl[iy][ix];
		}
		    iz = nzpad-nbd+ib;
		    uu2[iy][ix][iz] = uu1[iy][ix][iz-1] 
                        + (uu1[iy][ix][iz] - uu2[iy][ix][iz-1])*abc->bzh[iy][ix];
	    }
	    if (damp != NULL) {
		for (ib=0; ib<nbd-nop; ib++) {
		    damp_ib = damp[ib];
		    iz = nbd-ib-1;
		    if (abc->free) {
			uu2[iy][ix][iz] = 0.0f;
		    } else {
			uu2_bc = uu1[iy][ix][iz+1] + (uu1[iy][ix][iz] - uu2[iy][ix][iz+1])*abc->bzl[iy][ix];
			uu2[iy][ix][iz] = uu2_bc*(1.f - damp_ib) + uu2[iy][ix][iz]*damp_ib;
		    }
		    iz = nzpad-nbd+ib;
		    uu2_bc = uu1[iy][ix][iz-1] + (uu1[iy][ix][iz] - uu2[iy][ix][iz-1])*abc->bzh[iy][ix];
		    uu2[iy][ix][iz] = uu2_bc*(1.f - damp_ib) + uu2[iy][ix][iz]*damp_ib;
		}
	    }
	}
    }
    
#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(ix,iy,iz,ib,damp_ib,uu2_bc)
#endif
    for (iy=0; iy<nypad; iy++) {
	for (iz=0;iz<nzpad; iz++) {
	    for (ib=nbd-nop; ib<nbd; ib++) {
		ix = nbd-ib-1;
		uu2[iy][ix][iz] = uu1[iy][ix+1][iz] 
		    + (uu1[iy][ix][iz] - uu2[iy][ix+1][iz])*abc->bxl[iy][iz];
		ix = nxpad-nbd+ib;
		uu2[iy][ix][iz] = uu1[iy][ix-1][iz] 
		    + (uu1[iy][ix][iz] - uu2[iy][ix-1][iz])*abc->bxh[iy][iz];
	    }
	    if (damp != NULL) {
		for (ib=0; ib<nbd-nop; ib++) {
		    damp_ib = damp[ib];
		    ix = nbd-ib-1;
		    uu2_bc = uu1[iy][ix+1][iz] + (uu1[iy][ix][iz] - uu2[iy][ix+1][iz])*abc->bxl[iy][iz];
		    uu2[iy][ix][iz] = uu2_bc*(1.f - damp_ib) + uu2[iy][ix][iz]*damp_ib;
		    ix = nxpad-nbd+ib;
		    uu2_bc = uu1[iy][ix-1][iz] + (uu1[iy][ix][iz] - uu2[iy][ix-1][iz])*abc->bxh[iy][iz];
		    uu2[iy][ix][iz] = uu2_bc*(1.f - damp_ib) + uu2[iy][ix][iz]*damp_ib;
		}
	    }
	}
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) private(ix,iy,iz,ib,damp_ib,uu2_bc)
#endif
    for (ix=0; ix<nxpad; ix++) {
	for (iz=0;iz<nzpad; iz++) {
	    for (ib=nbd-nop; ib<nbd; ib++) {
		iy = nbd-ib-1;
		uu2[iy][ix][iz] = uu1[iy+1][ix][iz] 
		    + (uu1[iy][ix][iz] - uu2[iy+1][ix][iz])*abc->byl[ix][iz];
		iy = nypad-nbd+ib;
		uu2[iy][ix][iz] = uu1[iy-1][ix][iz] 
		    + (uu1[iy][ix][iz] - uu2[iy-1][ix][iz])*abc->byh[ix][iz];
	    }
	    if (damp != NULL) {
		for (ib=0; ib<nbd-nop; ib++) {
		    damp_ib = damp[ib];
		    iy = nbd-ib-1;
		    uu2_bc = uu1[iy+1][ix][iz] + (uu1[iy][ix][iz] - uu2[iy+1][ix][iz])*abc->byl[ix][iz];
		    uu2[iy][ix][iz] = uu2_bc*(1.f - damp_ib) + uu2[iy][ix][iz]*damp_ib;
		    iy = nypad-nbd+ib;
		    uu2_bc = uu1[iy-1][ix][iz] + (uu1[iy][ix][iz] - uu2[iy-1][ix][iz])*abc->byh[ix][iz];
		    uu2[iy][ix][iz] = uu2_bc*(1.f - damp_ib) + uu2[iy][ix][iz]*damp_ib;
		}
	    }
	}
    }
    return;
}
