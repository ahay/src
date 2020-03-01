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

typedef struct born_setup_struct born_setup_struct_t;
/*^*/

typedef enum adj_enum_t{FWD,ADJ} adj_t;
/*^*/

struct in_para_struct{
  bool verb;
  bool fsrf;
  bool dabc;
  bool snap;
  adj_t adj;
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
  float *v3c; // component 3 current
  float *v1p; // component 1 previous
  float *v2p; // component 2 previous
  float *v3p; // component 3 previous
  float *v1a; // component 1 aux
  float *v2a; // component 2 aux
  float *v3a; // component 3 aux
  float *tap1;  // ABC coeffs
  float *tap2;  //
  float *tap3;
  // buffers to store the data
  float *rdata;
  // buffer for the wavefield snapshot
  float *bwfl;
  // dimensions
  long nabc;  // size of the absorbing boundary
  long modN1;
  long modN2;
  long modN3;
  long simN1;
  long simN2;
  long simN3;
  float modO1;
  float modO2;
  float modO3;
  float simO1;
  float simO2;
  float simO3;
  float d1;
  float d2;
  float d3;
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
  float *hicksSou3;
  float *hicksRcv1;
  float *hicksRcv2;
  float *hicksRcv3;
};
/*^*/

struct mod_struct{
  long n1;
  long n2;
  long n3;
  float d1;
  float d2;
  float d3;
  float o1;
  float o2;
  float o3;
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

struct born_setup_struct{
  bool inputVelPerturbation;
  bool inputDenPerturbation;
  //
  bool outputBackgroundWfl;
  bool outputScatteredWfl;
  //
  bool outputBackgroundData;
  //
  bool outputVelPertImage;
  bool outputDenPertImage;
  //
  // I need this to compute the secondary sources
  char* bckwflfilename;
  FILE* Fbwfl;
  char* sctwflfilename;
  FILE* Fswfl;
  char* pv1wflfilename; // to compute the density perturbation in the adjoint
  char* pv2wflfilename;
  FILE* Fpv1;
  FILE* Fpv2;
};
/*^*/

#define NOP 3 /* derivative operator half-size */
/*^*/

#define IDX2D(i1,i2)((i1) + n1*(i2))
/*^*/

#define IDX3D(i1,i2,i3)( (i1) + n1*( (i2) + n2*(i3) ) )
/*^*/

/* LS coefficients */
#define C1 +1.1989919
/*^*/
#define C2 -0.08024696
/*^*/
#define C3 +0.00855954
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

static float* build_extended_model_3d(float const *mod, long n1, long n2, long n3, int next)
{
  long n1ext = n1+2*next;
  long n2ext = n2+2*next;
  long n3ext = n3+2*next;
  long nelem = n1ext*n2ext*n3ext;
  float* extd = sf_floatalloc(nelem);

  // core
  for (long i3=next,j3=0; i3<next+n3; i3++,j3++){
    for (long i2=next,j2=0; i2<next+n2; i2++,j2++){
      for (long i1=next,j1=0; i1<next+n1; i1++,j1++){
        float const v = mod[j1+n1*(j2+n2*j3)];
        extd[i1+n1ext*(i2+n2ext*i3)] = v;
      }
    }
  }

  // extend left and right
  for (long i3=next; i3<next+n3; i3++){
    for (long i2=0; i2<next; i2++){
      for (long i1=next; i1<next+n1; i1++){
        extd[i1+n1ext*(        i2+n2ext*i3)] = extd[i1+n1ext*(        next+n2ext*i3)];
        extd[i1+n1ext*(n2ext-1-i2+n2ext*i3)] = extd[i1+n1ext*(n2ext-1-next+n2ext*i3)];
      }
    }
  }

  // extend front and back
  for (long i3=0; i3<next; i3++){
    for (long i2=0; i2<n2ext; i2++){
      for (long i1=next; i1<next+n1; i1++){
        extd[i1+n1ext*(i2+n2ext*(        i3))] = extd[i1+n1ext*(i2+n2ext*(        next))];
        extd[i1+n1ext*(i2+n2ext*(n3ext-1-i3))] = extd[i1+n1ext*(i2+n2ext*(n3ext-1-next))];
      }
    }
  }

  // extend up and down
  for (long i3=0; i3<n3ext; i3++){
    for (long i2=0; i2<n2ext; i2++){
      for (long i1=0; i1<next; i1++){
        extd[        i1+n1ext*(i2+n2ext*i3)] = extd[        next+n1ext*(i2+n2ext*i3)];
        extd[n1ext-1-i1+n1ext*(i2+n2ext*i3)] = extd[n1ext-1-next+n1ext*(i2+n2ext*i3)];
      }
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
  if (para.verb) sf_warning("\tVelocity Model average value = %g",vave);

  mod->dmod = (float*) sf_floatalloc(nelem);
  sf_floatread(mod->dmod,nelem,Fdmod);

  double dave = 0;
  for (int i=0; i<nelem; i++)
    dave += mod->dmod[i];
  dave /= (nelem);
  if (para.verb) sf_warning("\tDensity Model average value = %g",dave);

  // modeling parameters
  long n1 = mod->n1;
  long n2 = mod->n2;
  int nabc = para.nb;
  if (para.verb) sf_warning("\tExtend the 2d models with absorbing boundaries..");
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

void prepare_model_3d(mod_struct_t* mod,
                      in_para_struct_t para,
                      sf_axis axvel[3], sf_axis axden[3],
                      sf_file Fvmod, sf_file Fdmod)
/*< Sets up the specified model cube >*/
{
  mod->n1 = sf_n(axvel[0]);
  mod->n2 = sf_n(axvel[1]);
  mod->n3 = sf_n(axvel[2]);
  mod->d1 = sf_d(axvel[0]);
  mod->d2 = sf_d(axvel[1]);
  mod->d3 = sf_d(axvel[2]);
  mod->o1 = sf_o(axvel[0]);
  mod->o2 = sf_o(axvel[1]);
  mod->o3 = sf_o(axvel[2]);

  long nelem = mod->n1*mod->n2*mod->n3;

  mod->vmod = (float*) sf_floatalloc(nelem);
  sf_floatread(mod->vmod,nelem,Fvmod);

  double vave = 0;
  for (int i=0; i<nelem; i++)
    vave += mod->vmod[i];
  vave /= (nelem);
  if (para.verb) sf_warning("\tVelocity Model average value = %g",vave);

  mod->dmod = (float*) sf_floatalloc(nelem);
  sf_floatread(mod->dmod,nelem,Fdmod);

  double dave = 0;
  for (int i=0; i<nelem; i++)
    dave += mod->dmod[i];
  dave /= (nelem);
  if (para.verb) sf_warning("\tDensity Model average value = %g",dave);

  // modeling parameters
  long n1 = mod->n1;
  long n2 = mod->n2;
  long n3 = mod->n3;
  int nabc = para.nb;
  if (para.verb) sf_warning("\tExtend the 2d models with absorbing boundaries..");
  mod->incomp = build_extended_model_3d(mod->vmod,n1,n2,n3,nabc);
  mod->buoy = build_extended_model_3d(mod->dmod,n1,n2,n3,nabc);

  // compute incompressibility and buoyancy
  long nelemext = (n1+2*nabc)*(n2+2*nabc)*(n3+2*nabc);
  for (long i=0; i<nelemext; i++){
    float v = mod->incomp[i];
    float d = mod->buoy[i];
    mod->incomp[i] = v*v*d;
    mod->buoy[i] = 1./d;
  }

}

void prepare_born_model_2d(mod_struct_t * const mod,
                            sf_axis axVel[2], sf_axis axDen[2],
                            sf_file Fvel, sf_file Fden,
                            sf_file Fvpert, sf_file Frpert,
                            in_para_struct_t in_para)
/*< Prepare the born operator model parameters >*/
{
  if (in_para.verb) sf_warning("\tPrepare the background model..");
  prepare_model_2d(mod,in_para,axVel,axDen,Fvel,Fden);

  if (in_para.verb) sf_warning("\tPrepare the perturbation cubes..");

  long n1 = sf_n(axVel[0]);
  long n2 = sf_n(axVel[1]);
  long n12 = n1*n2;

  mod->velpert = sf_floatalloc(n12);
  mod->denpert = sf_floatalloc(n12);

  memset(mod->velpert,0,n12*sizeof(float));
  memset(mod->denpert,0,n12*sizeof(float));

  if (in_para.adj==FWD){
    if (Fvpert) sf_floatread(mod->velpert,n12,Fvpert);
    if (Frpert) sf_floatread(mod->denpert,n12,Frpert);
  }

}

void make_born_velocity_sources_2d(wfl_struct_t * const wfl,
                                   mod_struct_t const * mod,
                                   acq_struct_t const * acq,
                                   born_setup_struct_t *para)
/*< Make the born sources for FWD born modelling>*/
{
  sf_warning("\tPARTICLE VELOCITY secondary sources (pressure gradient)..");

  long n1 = mod->n1;
  long n2 = mod->n2;

  long nt = acq->ntdat;
  float d1 = mod->d1;
  float d2 = mod->d2;

  wfl->Fprgrd = sf_tempfile(&(wfl->prtmpfilename),"w+");

  float *pp = sf_floatalloc(n1*n2);
  float *v1a= sf_floatalloc(n1*n2);
  float *v2a= sf_floatalloc(n1*n2);

  for (long it=0; it<nt; it++){
    // read the background wavefield
    if (para->outputBackgroundWfl)
      sf_floatread(pp,n1*n2,wfl->Fwfl);
    else
      fread(pp,sizeof(float),n1*n2,para->Fbwfl);

    for (long i2=0; i2<n2; i2++){
      for (long i1=2; i1<n1-3; i1++){
        v1a[i1+i2*n1] = (C1*(pp[i1+1+i2*n1] - pp[i1  +i2*n1])+
                         C2*(pp[i1+2+i2*n1] - pp[i1-1+i2*n1])+
                         C3*(pp[i1+3+i2*n1] - pp[i1-2+i2*n1]))/d1;
      }
    }

    for (long i2=2; i2<n2-3; i2++){
      for (long i1=0; i1<n1; i1++){
        v2a[i1+i2*n1] = (C1*(pp[i1+(i2+1)*n1] - pp[i1+i2    *n1])+
                         C2*(pp[i1+(i2+2)*n1] - pp[i1+(i2-1)*n1])+
                         C3*(pp[i1+(i2+3)*n1] - pp[i1+(i2-2)*n1]))/d2;
      }
    }

    fwrite(v1a,sizeof(float),n1*n2,wfl->Fprgrd);
    fwrite(v2a,sizeof(float),n1*n2,wfl->Fprgrd);

  }

  rewind(wfl->Fprgrd);
  // read the background wavefield
  if (para->outputBackgroundWfl)
    sf_seek(wfl->Fwfl,0,SEEK_SET);
  else
    rewind(para->Fbwfl);

  free(pp);
  free(v1a);
  free(v2a);

}

void make_born_pressure_sources_2d(wfl_struct_t * const wfl,
                                   mod_struct_t const * mod,
                                   acq_struct_t const * acq,
                                   born_setup_struct_t * para)
/*< Make the born sources for FWD born modelling>*/
{
  sf_warning("\tPRESSURE secondary sources (particle velocity divergence)..");

  long n1 = mod->n1;
  long n2 = mod->n2;
  long nelem = n1*n2;

  long nt = acq->ntdat;
  float dt = acq->dt;

  wfl->Fpvdiv = sf_tempfile(&(wfl->pvtmpfilename),"w+");

  float *snapc = sf_floatalloc(n1*n2);
  float *snapn = sf_floatalloc(n1*n2);

  if (para->outputBackgroundWfl){
    // read the background wavefield
    sf_floatread(snapc,nelem,wfl->Fwfl);

    for (long it=0; it<nt-1; it++){
      // read the background wavefield
      sf_floatread(snapn,nelem,wfl->Fwfl);

      // compute the divergence of the particle velocity from the pressure field
      for (long i=0; i<nelem; i++){
        float const v = mod->vmod[i];
        float const r = mod->dmod[i];
        float const scale = 1.f/(v*v*r)/dt;
        wfl->bwfl[i] = scale*(snapn[i] - snapc[i]);
      }
      fwrite(wfl->bwfl,nelem,sizeof(float),wfl->Fpvdiv);

      float *tmp = snapc;
      snapc = snapn;
      snapn = tmp;
    }

  }
  else{
    // read the background wavefield
    fread(snapc,sizeof(float),nelem,para->Fbwfl);

    for (long it=0; it<nt-1; it++){
      // read the background wavefield
      fread(snapn,sizeof(float),nelem,para->Fbwfl);

      // compute the divergence of the particle velocity from the pressure field
      for (long i=0; i<nelem; i++){
        float const v = mod->vmod[i];
        float const r = mod->dmod[i];
        float const scale = 1.f/(v*v*r)/dt;
        wfl->bwfl[i] = scale*(snapn[i] - snapc[i]);
      }
      fwrite(wfl->bwfl,nelem,sizeof(float),wfl->Fpvdiv);

      float *tmp = snapc;
      snapc = snapn;
      snapn = tmp;
    }
  }

  memset(wfl->bwfl,0,n1*n2*sizeof(float));
  fwrite(wfl->bwfl,n1*n2,sizeof(float),wfl->Fpvdiv);

  // rewind
  rewind(wfl->Fpvdiv);
  if (para->outputBackgroundWfl)
    sf_seek(wfl->Fwfl,0,SEEK_SET);
  else
    rewind(para->Fbwfl);

  free(snapc);
  free(snapn);

}

void stack_velocity_part_2d(wfl_struct_t * const wfl,
                            mod_struct_t const * mod,
                            acq_struct_t const * acq,
                            born_setup_struct_t *para)
/*< project the wavefields in the born model space >*/
{
  sf_warning("\tPARTICLE VELOCITY component of the density perturbation..");

  long nt = acq->ntdat;
  long n1 = mod->n1;
  long n2 = mod->n2;

  float dt = acq->dt;

  float* v1a = sf_floatalloc(n1*n2);
  float* v2a = sf_floatalloc(n1*n2);
  float* v1r = sf_floatalloc(n1*n2*nt); // back-propagated particle velocities
  float* v2r = sf_floatalloc(n1*n2*nt);
  float *rimg = sf_floatalloc(n1*n2);
  memset(rimg,0,n1*n2*sizeof(float));

  rewind(para->Fpv1);
  rewind(para->Fpv2);

  fread(v1r,sizeof(float),n1*n2*nt,para->Fpv1);
  fread(v2r,sizeof(float),n1*n2*nt,para->Fpv2);

  for (int it=0; it<nt; it++){
    float* w1p = v1r + (nt-1-it)*n1*n2;
    float* w2p = v2r + (nt-1-it)*n1*n2;

    // source side gradient of pressure
    fread(v1a,sizeof(float),n1*n2,wfl->Fprgrd);
    fread(v2a,sizeof(float),n1*n2,wfl->Fprgrd);

    for (long i=0; i<n1*n2; i++){
      v1a[i] *= -1.*w1p[i]; // flipping time flips the sign of the velocity
      v2a[i] *= -1.*w2p[i];
    }

    for (long i=0; i<n1*n2; i++){
      double r = mod->dmod[i];
      double scale = dt/(r*r);
      rimg[i] += (float) scale*(v1a[i]+v2a[i]);
    }

  }

  rewind(para->Fpv1);
  fwrite(rimg,sizeof(float),n1*n2,para->Fpv1);
  rewind(para->Fpv1);

  free(v1a);
  free(v2a);
  free(v1r);
  free(v2r);
  free(rimg);

}

void stack_pressure_part_2d(sf_file Fvpert,
                            sf_file Frpert,
                            wfl_struct_t * const wfl,
                            mod_struct_t const * mod,
                            acq_struct_t const * acq,
                            born_setup_struct_t *para)
/*< project the wavefields in the born model space>*/
{
  sf_warning("\tPRESSURE component of the velocity and density perturbations..");

  long nt = acq->ntdat;
  long n1 = mod->n1;
  long n2 = mod->n2;
  long nelem=n1*n2;

  float dt = acq->dt;

  float *srcwfl = sf_floatalloc(nelem*nt);
  float *tmp = sf_floatalloc(nelem);
  float *vimg = sf_floatalloc(nelem);

  // set
  memset(vimg,0,nelem*sizeof(float));

  fread(srcwfl,sizeof(float),nelem*nt,wfl->Fpvdiv);
  for (long it=0; it<nt; it++){
    float *wp = srcwfl + (nt-1-it)*nelem;

    if (para->outputScatteredWfl)
      sf_floatread(tmp,nelem,wfl->Fswfl);
    else
      fread(tmp,sizeof(float),nelem,para->Fswfl);

    for (long i=0; i<nelem; i++){
      double const v = mod->vmod[i];
      double const r = mod->dmod[i];
      double scale = 2.f*v*r*dt;
      vimg[i] += (float) scale*tmp[i]*wp[i];
    }
  }

  sf_floatwrite(vimg,nelem,Fvpert);

  if (para->outputDenPertImage){

    if (para->outputScatteredWfl)
      sf_seek(wfl->Fswfl,0,SEEK_SET);
    else
      rewind(para->Fswfl);

    float *rimg = sf_floatalloc(nelem);
    fread(rimg,sizeof(float),nelem,para->Fpv1);

    for (long it=0; it<nt; it++){
      float *wp = srcwfl + (nt-1-it)*nelem;

      if (para->outputScatteredWfl)
        sf_floatread(tmp,nelem,wfl->Fswfl);
      else
        fread(tmp,sizeof(float),nelem,para->Fswfl);

      for (long i=0; i<nelem; i++){
        double const v = mod->vmod[i];
        double const scale = v*v*dt;
        rimg[i] += (float) scale*tmp[i]*wp[i];
      }
    }

    sf_floatwrite(rimg,nelem,Frpert);
    free(rimg);
  }

  free(srcwfl);
  free(vimg);
  free(tmp);

}

void clear_model(mod_struct_t* mod)
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

  clear_model(mod);

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
  if (para.verb) sf_warning("\t Number of wavefield snapshots = %d",acq->ntsnap);
  if (para.verb) sf_warning("\t Number of data timesteps simulated = %d",acq->ntdat);

  long nsouwav = sf_n(axwav[0]);
  long nwavsamp = nsouwav*acq->ntdat;

  if (para.verb) {
    if (nsouwav==1){
      sf_warning("\t Using the same wavelet for all shots!");
    }
    else{
      if (nsouwav==acq->ns){
        sf_warning("\t Every shot has a different wavelet");
      }
      else{
        sf_error("Inconsistent number of wavelets and shots");
        return;
      }
    }
  }


  acq->wav = sf_floatalloc(nwavsamp);
  memset(acq->wav,0,nwavsamp*sizeof(float));
  for (int isou=0; isou<nsouwav; isou++)
    sf_floatread(acq->wav+isou*acq->ntdat,acq->nt,Fwav);

}

void prepare_acquisition_3d( acq_struct_t* acq, in_para_struct_t para,
    sf_axis axsou[2], sf_axis axrec[2], sf_axis axwav[2],
    sf_file Fsou, sf_file Frec, sf_file Fwav)
/*< Read the acquisition geometry from files >*/
{
  //
  if (sf_n(axsou[0])!=3)
    sf_error("Wrong number of coordinates in the source file!");

  // shot coordinates
  acq->ns = sf_n(axsou[1]);
  acq->scoord = sf_floatalloc(3*acq->ns);
  for (long isou=0; isou<acq->ns; isou++)
    sf_floatread(acq->scoord+3*isou,3,Fsou);

  //
  if (sf_n(axrec[0])!=3)
    sf_error("Wrong number of coordinates in the receiver file!");

  // receiver coordinates
  acq->nr = sf_n(axrec[1]);
  acq->rcoord = sf_floatalloc(3*acq->nr);
  for (long irec=0; irec<acq->nr; irec++)
    sf_floatread(acq->rcoord+3*irec,3,Frec);

  // wavelet parameters
  acq->nt = sf_n(axwav[1]);
  acq->dt = sf_d(axwav[1]);
  acq->ot = sf_o(axwav[1]);
  acq->ntdat = (acq->nt/para.jsnap)*para.jsnap;
  acq->ntsnap= (acq->nt+para.jsnap-1)/para.jsnap;
  if (para.verb) sf_warning("\t Number of wavefield snapshots = %d",acq->ntsnap);
  if (para.verb) sf_warning("\t Number of data timesteps simulated = %d",acq->ntdat);

  long nsouwav = sf_n(axwav[0]);
  long nwavsamp = nsouwav*acq->ntdat;

  if (para.verb) {
    if (nsouwav==1){
      sf_warning("\t Using the same wavelet for all shots!");
    }
    else{
      if (nsouwav==acq->ns){
        sf_warning("\t Every shot has a different wavelet");
      }
      else{
        sf_error("Inconsistent number of wavelets and shots");
        return;
      }
    }
  }


  acq->wav = sf_floatalloc(nwavsamp);
  memset(acq->wav,0,nwavsamp*sizeof(float));
  for (int isou=0; isou<nsouwav; isou++)
    sf_floatread(acq->wav+isou*acq->ntdat,acq->nt,Fwav);

}

void prepare_born_acquisition_2d(acq_struct_t * const acq,
                                 sf_axis axsou[2], sf_axis axrec[2], sf_axis axwav[2],
                                 sf_file Fsou, sf_file Frec, sf_file Fwav,
                                 sf_file Fsdat,
                                 in_para_struct_t in_para)
/*< preparation of the acquisition for born operator >*/
{
  prepare_acquisition_2d(acq, in_para, axsou, axrec, axwav, Fsou, Frec,Fwav);

  if (in_para.adj){
    // Read the data to backproject
    long nr = acq->nr;
    long nt = acq->ntdat;

    acq->dat = sf_floatalloc(nr*nt);
    for(long it=0; it<nt; it++){
      float* wp = acq->dat + (nt-1-it)*nr;
      sf_floatread(wp,nr,Fsdat);
    }
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

void set_sr_interpolation_coeffs_2d(acq_struct_t * const acq, wfl_struct_t const * wfl)
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

void set_sr_interpolation_coeffs_3d(acq_struct_t * const acq, wfl_struct_t const * wfl)
/*< interpolation coefficients for source injection and receiver extraction >*/
{
  long nsous = acq->ns;
  long nrecs = acq->nr;

  float o1 = wfl->simO1;
  float o2 = wfl->simO2;
  float o3 = wfl->simO3;
  float d1 = wfl->d1;
  float d2 = wfl->d2;
  float d3 = wfl->d3;

  acq->hicksSou1 = sf_floatalloc(8*nsous);
  acq->hicksSou2 = sf_floatalloc(8*nsous);
  acq->hicksSou3 = sf_floatalloc(8*nsous);
  acq->hicksRcv1 = sf_floatalloc(8*nrecs);
  acq->hicksRcv2 = sf_floatalloc(8*nrecs);
  acq->hicksRcv3 = sf_floatalloc(8*nrecs);

  for (long isou=0; isou<nsous; isou++){
    float x1s = acq->scoord[3*isou+2]; // z coordinate
    float x3s = acq->scoord[3*isou+1]; // y coordinate
    float x2s = acq->scoord[3*isou  ]; // x coordinate

    long ix1s = (x1s-o1)/d1;
    long ix2s = (x2s-o2)/d2;
    long ix3s = (x3s-o3)/d3;
    float rem1 = (x1s - (ix1s*d1+o1))/d1;
    float rem2 = (x2s - (ix2s*d2+o2))/d2;
    float rem3 = (x3s - (ix3s*d3+o3))/d3;
    for (int i=-3,ii=0; i<=4; i++,ii++){
      acq->hicksSou1[ii+isou*8] = sinc(i-rem1)*kwin(i-rem1,4.5);
      acq->hicksSou2[ii+isou*8] = sinc(i-rem2)*kwin(i-rem2,4.5);
      acq->hicksSou3[ii+isou*8] = sinc(i-rem3)*kwin(i-rem3,4.5);
    }
  }

  for (long irec=0; irec<nrecs; irec++){
    float x1r = acq->rcoord[3*irec+2]; // z coordinate
    float x3r = acq->rcoord[3*irec+1]; // z coordinate
    float x2r = acq->rcoord[3*irec  ]; // x coordinate

    long ix1r = (x1r-o1)/d1;
    long ix2r = (x2r-o2)/d2;
    long ix3r = (x3r-o3)/d3;
    float rem1 = (x1r - (ix1r*d1+o1))/d1;
    float rem2 = (x2r - (ix2r*d2+o2))/d2;
    float rem3 = (x3r - (ix3r*d3+o3))/d3;
    for (int i=-3,ii=0; i<=4; i++,ii++){
      acq->hicksRcv1[ii+irec*8] = sinc(i-rem1)*kwin(i-rem1,4.5);
      acq->hicksRcv2[ii+irec*8] = sinc(i-rem2)*kwin(i-rem2,4.5);
      acq->hicksRcv3[ii+irec*8] = sinc(i-rem3)*kwin(i-rem3,4.5);
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

void clear_acq_3d(acq_struct_t *acq)
/*< Free the arrays in the acquisition structure >*/
{

  free(acq->hicksRcv1);
  free(acq->hicksRcv2);
  free(acq->hicksRcv3);

  free(acq->hicksSou1);
  free(acq->hicksSou2);
  free(acq->hicksSou3);

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

void prepare_wfl_3d(wfl_struct_t *wfl,mod_struct_t *mod, sf_file Fdata, sf_file Fwfl, in_para_struct_t para)
/*< Allocate the wavefield structure >*/
{
  wfl->modN1 = mod->n1;
  wfl->modN2 = mod->n2;
  wfl->modN3 = mod->n3;

  wfl->nabc = para.nb;
  wfl->simN1 = mod->n1 + 2*para.nb;
  wfl->simN2 = mod->n2 + 2*para.nb;
  wfl->simN3 = mod->n3 + 2*para.nb;

  wfl->d1 = mod->d1;
  wfl->d2 = mod->d2;
  wfl->d3 = mod->d3;

  wfl->modO1 = mod->o1;
  wfl->modO2 = mod->o2;
  wfl->modO3 = mod->o3;
  wfl->simO1 = mod->o1 - para.nb*mod->d1;
  wfl->simO2 = mod->o2 - para.nb*mod->d2;
  wfl->simO3 = mod->o3 - para.nb*mod->d3;

  long nelem = wfl->simN1*wfl->simN2*wfl->simN3;
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

  wfl->v3c = sf_floatalloc(nelem);
  wfl->v3p = sf_floatalloc(nelem);
  wfl->v3a = sf_floatalloc(nelem);
  memset(wfl->v3c,0,wflsize);
  memset(wfl->v3p,0,wflsize);
  memset(wfl->v3a,0,wflsize);

  // sponge
  wfl->tap1 = sf_floatalloc(wfl->simN1);
  wfl->tap2 = sf_floatalloc(wfl->simN2);
  wfl->tap3 = sf_floatalloc(wfl->simN3);

  wfl->bwfl = sf_floatalloc(wfl->modN1*wfl->modN2*wfl->modN3);

  wfl->freesurf = para.fsrf;

  wfl->Fdata = Fdata;
  wfl->Fwfl  = Fwfl;
  wfl->snap = para.snap;
  wfl->jsnap=para.jsnap;
}

void prepare_born_wfl_2d( wfl_struct_t * const wfl, mod_struct_t * mod,
                          sf_file Fdata, sf_file Fwfl,
                          sf_file Fsdat, sf_file Fswfl,
                          in_para_struct_t para)
/*< Set up the extra parameters for the born operator>*/
{
  prepare_wfl_2d(wfl,mod,Fdata,Fwfl,para);

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

void clear_wfl_3d(wfl_struct_t *wfl)
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

  free(wfl->v3c);
  free(wfl->v3p);
  free(wfl->v3a);

  free(wfl->tap1);
  free(wfl->tap2);
  free(wfl->tap3);

  free(wfl->bwfl);
}

void clear_born_wfl_2d(wfl_struct_t *wfl)
/*< clear the source for born modeling >*/
{
  clear_wfl_2d(wfl);

  if (wfl->Fpvdiv)
    remove(wfl->pvtmpfilename);
  if (wfl->Fprgrd)
    remove(wfl->prtmpfilename);
}
