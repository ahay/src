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
  bool adj;
  int nb;
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
  // dimensions
  long n1;
  long n2;
  float d1;
  float d2;
  float o1;
  float o2;
  // pointers to files
  sf_file Fdata;
  sf_file Fwfl;
};
/*^*/

struct acq_struct{
  // wavelet and time parameters
  int nt;
  float dt;
  float ot;
  float *wav;
  // acquisition geometry
  long ns;    // NUMBER OF SOURCES TO BE SIMULATED SIMULTANEOUSLY
  long nr;
  float *scoord;
  float *rcoord;
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
};
/*^*/

#endif

void prepare_model_2d(mod_struct_t* mod,
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
  sf_warning("Velocity Model average value = %g",dave);

  // modeling parameters
  mod->incomp = sf_floatalloc(nelem);
  mod->buoy   = sf_floatalloc(nelem);

  for (long i=0; i<nelem; i++){
    float const v = mod->vmod[i];
    float const r = mod->dmod[i];
    mod->incomp[i] = v*v*r;
    mod->buoy[i]   = 1./r;
  }

}

void clear_model_2d(mod_struct_t* mod)
/*< Free the model parameter cubes >*/
{
  free(mod->vmod);
  free(mod->dmod);
  free(mod->incomp);
  free(mod->buoy);
}

void prepare_acquisition_2d( acq_struct_t* acq,
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
  acq->ot = sf_d(axwav[1]);

  long nsouwav = sf_n(axwav[0]);
  long nwavsamp = nsouwav*acq->nt;

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
  sf_floatread(acq->wav,nwavsamp,Fwav);

}

void clear_acq_2d(acq_struct_t *acq)
/*< Free the arrays in the acquisition structure >*/
{
  free(acq->scoord);
  free(acq->rcoord);

  free(acq->wav);
}

void prepare_wfl_2d(wfl_struct_t *wfl,mod_struct_t *mod, sf_file Fdata, sf_file Fwfl, in_para_struct_t para)
/*< Allocate the wavefield structure >*/
{
  // FIXME: the wavefield model need to be extended for absorbing boundaries
  wfl->n1 = mod->n1;
  wfl->n2 = mod->n2;

  wfl->d1 = mod->d1;
  wfl->d2 = mod->d2;

  wfl->o1 = mod->o1;
  wfl->o2 = mod->o2;

  long nelem = wfl->n1*wfl->n2;
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

  wfl->Fdata = Fdata;
  wfl->Fwfl  = Fwfl;

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

}
