#include <rsf.h>

#ifndef _PREP_H

typedef struct in_para_struct in_para_struct_t;
/*^*/

typedef struct wfl_struct wfl_struct_t;
/*^*/

typedef struct acq_struct acq_struct_t;
/*^*/

typedef struct mod_struct mod_struct_t;
/*^*/

struct in_para_struct{
  bool verb;
  bool fsrf;
  bool dabc;
  bool snap;
  bool adj;
  bool dpt;
  int nb;
  int jsnap;
};
/*^*/

struct wfl_struct{
  // pressure
  float *pc;  // current time
  float *pp;  // previous time
  float *pa;  // aux
  // velocities
  float *v1c; // component 1 current
  float *v2c; // component 2 current
  float *v1p; // component 1 previous
  float *v2p; // component 2 previous
  float *v1a; // component 1 aux
  float *v2a; // component 1 aux
  float *tap1;  // ABC coeffs
  float *tap2;  //
  // buffers to store the data
  float *rdata;
  // buffer for the wavefield snapshot
  float *bwfl;
  // dimensions
  long nabc;  // size of the absorbing boundary
  long modN1;
  long modN2;
  long simN1;
  long simN2;
  float modO1;
  float modO2;
  float simO1;
  float simO2;
  float d1;
  float d2;
  //free surface flag
  bool freesurf;
  //wavefield snapshots
  bool snap;
  int jsnap;
  // pointers to files
  sf_file Fdata;
  sf_file Fwfl;
  // born modeling (scattered data)
  sf_file Fsdata;
  sf_file Fswfl;
  char* pvtmpfilename;
  char* prtmpfilename;
  FILE* Fpvdiv;
  FILE* Fprgrd;
};
/*^*/

struct acq_struct{
  // wavelet and time parameters
  long nt;
  long ntdat;  // number of samples to have an intenger multiple of the undersampling of the wavefield
  long ntsnap;
  float dt;
  float ot;
  float *wav;
  float *dat;
  // acquisition geometry
  long ns;    // NUMBER OF SOURCES TO BE SIMULATED SIMULTANEOUSLY
  long nr;
  float *scoord;
  float *rcoord;
  // interpolation coefficients for source injection and extraction
  float *hicksSou1;
  float *hicksSou2;
  float *hicksRcv1;
  float *hicksRcv2;
};
/*^*/

struct mod_struct{
  long n1;
  long n2;
  float d1;
  float d2;
  float o1;
  float o2;
  float *vmod;
  float *dmod;
  // parameters for modeling
  float *incomp;
  float *buoy;
  // perturnations (for born modeling)
  float *velpert;
  float *denpert;
};
/*^*/

#endif

void init_param(in_para_struct_t* par)
/*< parameter initialization >*/
{

  par->adj=false;
  par->dabc=true;
  par->dpt=false;
  par->fsrf=false;
  par->snap=false;
  par->verb=true;
  par->jsnap=1;
  par->nb=1;
}

void print_param(in_para_struct_t in_para)
/*< display the parameters >*/
{
  sf_warning("verbosity          = %s",((in_para.verb==false)?"no":"yes"));
  sf_warning("free surface       = %s",((in_para.fsrf==false)?"no":"yes"));
  sf_warning("absorbing boundary = %s",((in_para.dabc==false)?"no":"yes"));
  if (in_para.dabc) sf_warning("- sponge thickness = %d",in_para.nb);
  sf_warning("wavefield snapshots= %s",((in_para.snap==false)?"no":"yes"));
  if (in_para.snap) sf_warning("- wavefield time undersampling = %d",in_para.jsnap);
}

static float* build_extended_model_2d(float const *mod, long n1, long n2, int next)
{

  sf_warning(" Extend the 2d models with absorbing boundaries..");

  long n1ext = n1+2*next;
  long n2ext = n2+2*next;
  long nelem = n1ext*n2ext;
  float* extd = sf_floatalloc(nelem);

  // core
  for (long i2=next,j2=0; i2<next+n2; i2++,j2++){
    for (long i1=next,j1=0; i1<next+n1; i1++,j1++){
      float const v = mod[j1+j2*n1];
      extd[i1+i2*n1ext] = v;
    }
  }

  // extend laterally
  for (long i2=0;i2<next;i2++){
    for (long i1=0,j1=next; i1<n1; i1++,j1++){
      extd[j1+i2*n1ext] = extd[j1+next*n1ext];
      extd[j1+(n2ext-1-i2)*n1ext] = extd[j1+(n2ext-1-next)*n1ext];
    }
  }

  // extend up and down
  for (long i2=0;i2<n2ext;i2++){
    for (long i1=0; i1<next; i1++){
      extd[i1+i2*n1ext] = extd[next+i2*n1ext];
      extd[(n1ext-1-i1)+i2*n1ext] = extd[(n1ext-1-next)+i2*n1ext];
    }
  }

  return extd;

}

void prepare_model_2d(mod_struct_t* mod,
                      in_para_struct_t para,
                      sf_axis axvel[2], sf_axis axden[2],
                      sf_file Fvmod, sf_file Fdmod)
/*< Sets up the specified model cube >*/
{

  if ((sf_n(axvel[0])!=sf_n(axden[0])) || (sf_n(axvel[1])!=sf_n(axden[1])))
    sf_error("Inconsistent model dimensions!");

  mod->n1 = sf_n(axvel[0]);
  mod->n2 = sf_n(axvel[1]);
  mod->d1 = sf_d(axvel[0]);
  mod->d2 = sf_d(axvel[1]);
  mod->o1 = sf_o(axvel[0]);
  mod->o2 = sf_o(axvel[1]);

  long nelem = mod->n1*mod->n2;

  mod->vmod = (float*) sf_floatalloc(nelem);
  sf_floatread(mod->vmod,nelem,Fvmod);

  double vave = 0;
  for (int i=0; i<nelem; i++)
    vave += mod->vmod[i];
  vave /= (nelem);
  sf_warning("Velocity Model average value = %g",vave);

  mod->dmod = (float*) sf_floatalloc(nelem);
  sf_floatread(mod->dmod,nelem,Fdmod);

  double dave = 0;
  for (int i=0; i<nelem; i++)
    dave += mod->dmod[i];
  dave /= (nelem);
  sf_warning("Density Model average value = %g",dave);

  // modeling parameters
  long n1 = mod->n1;
  long n2 = mod->n2;
  int nabc = para.nb;
  mod->incomp = build_extended_model_2d(mod->vmod,n1,n2,nabc);
  mod->buoy = build_extended_model_2d(mod->dmod,n1,n2,nabc);

  // compute incompressibility and buoyancy
  long nelemext = (n1+2*nabc)*(n2+2*nabc);
  for (long i=0; i<nelemext; i++){
    float v = mod->incomp[i];
    float d = mod->buoy[i];
    mod->incomp[i] = v*v*d;
    mod->buoy[i] = 1./d;
  }

}

void prepare_born_model_2d(mod_struct_t * const mod,
                           sf_axis axVel[2],
                           sf_file Fvpert,
                           sf_file Frpert)
/*< Prepare the born operator model parameters >*/
{
  sf_warning(" Read the model perturbation files..");

  long n1 = sf_n(axVel[0]);
  long n2 = sf_n(axVel[1]);
  long n12 = n1*n2;

  mod->velpert = sf_floatalloc(n12);
  mod->denpert = sf_floatalloc(n12);

  memset(mod->velpert,0,n12*sizeof(float));
  memset(mod->denpert,0,n12*sizeof(float));

  if (Fvpert)
    sf_floatread(mod->velpert,n12,Fvpert);

  if (Frpert)
    sf_floatread(mod->denpert,n12,Frpert);

}

void make_born_sources_2d(wfl_struct_t * const wfl, mod_struct_t const * mod, acq_struct_t const * acq)
/*< Make the born sources for FWD born modelling>*/
{
  long n1 = mod->n1;
  long n2 = mod->n2;

  long nt = acq->ntdat;
  float dt = acq->dt;

  wfl->Fpvdiv = sf_tempfile(&(wfl->pvtmpfilename),"w+");

  float *snapc = sf_floatalloc(n1*n2);
  float *snapn = sf_floatalloc(n1*n2);

  sf_floatread(snapc,n1*n2,wfl->Fwfl);
  for (long it=0; it<nt-1; it++){

    sf_floatread(snapn,n1*n2,wfl->Fwfl);

    // compute the divergence of the particle velocity from the pressure field
    for (long i=0; i<n1*n2; i++){
      float v = mod->vmod[i];
      float r = mod->dmod[i];
      float k = 1.f/(v*v*r);
      wfl->bwfl[i] = k*(snapn[i] - snapc[i])/dt;
    }
    fwrite(wfl->bwfl,n1*n2,sizeof(float),wfl->Fpvdiv);

    float *tmp = snapc;
    snapc = snapn;
    snapn = tmp;
  }

  memset(wfl->bwfl,0,n1*n2*sizeof(float));
  fwrite(wfl->bwfl,n1*n2,sizeof(float),wfl->Fpvdiv);

  rewind(wfl->Fpvdiv);

  free(snapc);
  free(snapn);

}

void stack_wfl_2d(sf_file Fvpert, sf_file Frpert, wfl_struct_t * const wfl, mod_struct_t const * mod, acq_struct_t const * acq)
/*< project the wavefields in the born model space>*/
{
  long nt = acq->ntdat;
  long n1 = mod->n1;
  long n2 = mod->n2;

  float dt = acq->dt;

  float *srcwfl = sf_floatalloc(n1*n2*nt);
  float *tmp = sf_floatalloc(n1*n2);

  fread(srcwfl,sizeof(float),n1*n2*nt,wfl->Fpvdiv);
  for (long it=0; it<nt; it++){
    float *wp = srcwfl + (nt-1-it)*n1*n2;
    sf_floatread(tmp,n1*n2,wfl->Fswfl);

    for (long i=0; i<n1*n2; i++)
      wp[i] = wp[i]*tmp[i];

  }

  //stack the wavefield
  memset(tmp,0,n1*n2*sizeof(float));
  for (long it=0; it<nt; it++){
    float *wp = srcwfl + it*n1*n2;

    for (long i=0; i<n1*n2; i++){
      float v = mod->vmod[i];
      float r = mod->dmod[i];
      tmp[i]+= 2.f*v*r*wp[i]*dt;
    }
  }

  sf_floatwrite(tmp,n1*n2,Fvpert);

  free(srcwfl);
  free(tmp);

}

void clear_model_2d(mod_struct_t* mod)
/*< Free the model parameter cubes >*/
{
  free(mod->vmod);
  free(mod->dmod);
  free(mod->incomp);
  free(mod->buoy);
}

void clear_born_model_2d(mod_struct_t* mod)
/*< Free the model perturbation parameters cubes for the born operator >*/
{
  free(mod->velpert);
  free(mod->denpert);
}

void prepare_acquisition_2d( acq_struct_t* acq, in_para_struct_t para,
                             sf_axis axsou[2], sf_axis axrec[2], sf_axis axwav[2],
                             sf_file Fsou, sf_file Frec, sf_file Fwav)
/*< Read the acquisition geometry from files >*/
{
  //
  if (sf_n(axsou[0])!=2)
    sf_error("Wrong number of coordinates in the source file!");

  // shot coordinates
  acq->ns = sf_n(axsou[1]);
  acq->scoord = sf_floatalloc(2*acq->ns);
  for (long isou=0; isou<acq->ns; isou++)
    sf_floatread(acq->scoord+2*isou,2,Fsou);

  //
  if (sf_n(axrec[0])!=2)
    sf_error("Wrong number of coordinates in the receiver file!");

  // receiver coordinates
  acq->nr = sf_n(axrec[1]);
  acq->rcoord = sf_floatalloc(2*acq->nr);
  for (long irec=0; irec<acq->nr; irec++)
    sf_floatread(acq->rcoord+2*irec,2,Frec);

  // wavelet parameters
  acq->nt = sf_n(axwav[1]);
  acq->dt = sf_d(axwav[1]);
  acq->ot = sf_o(axwav[1]);
  acq->ntdat = (acq->nt/para.jsnap)*para.jsnap;
  acq->ntsnap= (acq->nt+para.jsnap-1)/para.jsnap;
  sf_warning("Ntsnap = %d\n",acq->ntsnap);
  sf_warning("Ntdat  = %d\n",acq->ntdat);

  long nsouwav = sf_n(axwav[0]);
  long nwavsamp = nsouwav*acq->ntdat;

  if (nsouwav==1){
    sf_warning("Using the same wavelet for all shots!");
  }
  else{
    if (nsouwav==acq->ns){
      sf_warning("Every shot has a different wavelet");
    }
    else{
      sf_error("Inconsistent number of wavelets and shots");
      return;
    }
  }

  acq->wav = sf_floatalloc(nwavsamp);
  memset(acq->wav,0,nwavsamp*sizeof(float));
  for (int isou=0; isou<nsouwav; isou++)
    sf_floatread(acq->wav+isou*acq->ntdat,acq->nt,Fwav);

}

void prepare_scatt_data_2d(acq_struct_t * acq,sf_file Fsdat)
/*< prepare the scattered data for backward extrapolation >*/
{
  long nr = acq->nr;
  long nt = acq->ntdat;

  acq->dat = sf_floatalloc(nr*nt);
  for(long it=0; it<nt; it++){
    float* wp = acq->dat + (nt-1-it)*nr;
    sf_floatread(wp,nr,Fsdat);
  }
}

/*
 * Evaluates the Bessel function (from Dave Hale's KaiserWindow class in the Java Mines TK)
 */
static double ino(double x)
{
  double s = 1.0;
  double ds = 1.0;
  double d = 0.0;
  do {
    d += 2.0;
    ds *= (x*x)/(d*d);
    s += ds;
  } while (ds>s*DBL_EPSILON);
  return s;
}

/*
 * Sinc function
 */
static double sinc(double x){

  if (x==0.0)
    return 1;
  else
    return sin(M_PI*x)/(M_PI*x);
}

/*
 * Kaiser window with fixed scale parameter and radius
 */
static double kwin(double x, double xmax){
  double b = 6.31;
  double xx = x*x;
  double xxmax = xmax*xmax;
  double scale = 1./ino(b);

  return (xx<xxmax)? scale*ino(b*sqrt(1.- xx/xxmax)):0.;

}

void set_sr_interpolation_coeffs(acq_struct_t * const acq, wfl_struct_t const * wfl)
/*< interpolation coefficients for source injection and receiver extraction >*/
{
  long nsous = acq->ns;
  long nrecs = acq->nr;

  float o1 = wfl->simO1;
  float o2 = wfl->simO2;
  float d1 = wfl->d1;
  float d2 = wfl->d2;

  acq->hicksSou1 = sf_floatalloc(8*nsous);
  acq->hicksSou2 = sf_floatalloc(8*nsous);
  acq->hicksRcv1 = sf_floatalloc(8*nrecs);
  acq->hicksRcv2 = sf_floatalloc(8*nrecs);

  for (long isou=0; isou<nsous; isou++){
    float x1s = acq->scoord[2*isou+1]; // z coordinate
    float x2s = acq->scoord[2*isou  ]; // x coordinate
    long ix1s = (x1s-o1)/d1;
    long ix2s = (x2s-o2)/d2;
    float rem1 = (x1s - (ix1s*d1+o1))/d1;
    float rem2 = (x2s - (ix2s*d2+o2))/d2;
    for (int i=-3,ii=0; i<=4; i++,ii++){
      acq->hicksSou1[ii+isou*8] = sinc(i-rem1)*kwin(i-rem1,4.5);
      acq->hicksSou2[ii+isou*8] = sinc(i-rem2)*kwin(i-rem2,4.5);
    }
  }

  for (long irec=0; irec<nrecs; irec++){
    float x1r = acq->rcoord[2*irec+1]; // z coordinate
    float x2r = acq->rcoord[2*irec  ]; // x coordinate
    long ix1r = (x1r-o1)/d1;
    long ix2r = (x2r-o2)/d2;
    float rem1 = (x1r - (ix1r*d1+o1))/d1;
    float rem2 = (x2r - (ix2r*d2+o2))/d2;
    for (int i=-3,ii=0; i<=4; i++,ii++){
      acq->hicksRcv1[ii+irec*8] = sinc(i-rem1)*kwin(i-rem1,4.5);
      acq->hicksRcv2[ii+irec*8] = sinc(i-rem2)*kwin(i-rem2,4.5);
    }
  }

}

void clear_acq_2d(acq_struct_t *acq)
/*< Free the arrays in the acquisition structure >*/
{

  free(acq->hicksRcv1);
  free(acq->hicksRcv2);

  free(acq->hicksSou1);
  free(acq->hicksSou2);

  free(acq->scoord);
  free(acq->rcoord);

  free(acq->wav);
}

void prepare_wfl_2d(wfl_struct_t *wfl,mod_struct_t *mod, sf_file Fdata, sf_file Fwfl, in_para_struct_t para)
/*< Allocate the wavefield structure >*/
{

  wfl->modN1 = mod->n1;
  wfl->modN2 = mod->n2;

  wfl->nabc = para.nb;
  wfl->simN1 = mod->n1 + 2*para.nb;
  wfl->simN2 = mod->n2 + 2*para.nb;

  wfl->d1 = mod->d1;
  wfl->d2 = mod->d2;

  wfl->modO1 = mod->o1;
  wfl->modO2 = mod->o2;
  wfl->simO1 = mod->o1 - para.nb*mod->d1;
  wfl->simO2 = mod->o2 - para.nb*mod->d2;

  long nelem = wfl->simN1*wfl->simN2;
  long wflsize = nelem*sizeof(float);
  wfl->pc = sf_floatalloc(nelem);
  wfl->pp = sf_floatalloc(nelem);
  wfl->pa = sf_floatalloc(nelem);
  memset(wfl->pc,0,wflsize);
  memset(wfl->pp,0,wflsize);
  memset(wfl->pa,0,wflsize);

  wfl->v1c = sf_floatalloc(nelem);
  wfl->v1p = sf_floatalloc(nelem);
  wfl->v1a = sf_floatalloc(nelem);
  memset(wfl->v1c,0,wflsize);
  memset(wfl->v1p,0,wflsize);
  memset(wfl->v1a,0,wflsize);

  wfl->v2c = sf_floatalloc(nelem);
  wfl->v2p = sf_floatalloc(nelem);
  wfl->v2a = sf_floatalloc(nelem);
  memset(wfl->v2c,0,wflsize);
  memset(wfl->v2p,0,wflsize);
  memset(wfl->v2a,0,wflsize);

  // sponge
  wfl->tap1 = sf_floatalloc(wfl->simN1);
  wfl->tap2 = sf_floatalloc(wfl->simN2);

  wfl->bwfl = sf_floatalloc(wfl->modN1*wfl->modN2);

  wfl->freesurf = para.fsrf;

  wfl->Fdata = Fdata;
  wfl->Fwfl  = Fwfl;
  wfl->snap = para.snap;
  wfl->jsnap=para.jsnap;

}

void prepare_born_wfl_2d(wfl_struct_t * const wfl,sf_file Fsdat, sf_file Fswfl)
/*< Set up the extra parameters for the born operator>*/
{
  wfl->Fsdata = Fsdat;
  wfl->Fswfl  = Fswfl;

  wfl->Fpvdiv = NULL;
  wfl->Fprgrd = NULL;
}

void clear_wfl_2d(wfl_struct_t *wfl)
/*< Clear the wavefield structure >*/
{

  free(wfl->pc);
  free(wfl->pp);
  free(wfl->pa);

  free(wfl->v1c);
  free(wfl->v1p);
  free(wfl->v1a);

  free(wfl->v2c);
  free(wfl->v2p);
  free(wfl->v2a);

  free(wfl->tap1);
  free(wfl->tap2);

  free(wfl->bwfl);

}

void clear_born_wfl_2d(wfl_struct_t *wfl)
/*< clear the source for born modeling >*/
{
  if (wfl->Fpvdiv)
    remove(wfl->pvtmpfilename);
}
