#include <rsf.h>

#ifndef _PREP_H

typedef struct wfl_struct wfl_struct_t;
/*^*/

typedef struct acq_struct acq_struct_t;
/*^*/

typedef struct mod_struct mod_struct_t;
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
};
/*^*/

struct acq_struct{
  // wavelet and time parameters
  int nt;
  float dt;
  float ot;
  float *wav;
  // acquisition geometry
  long ns;
  long  *nr;
  float *xs;
  float *ys;
  float *zs;
  float *xr;
  float *yr;
  float *zr;
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

}

void clear_model_2d(mod_struct_t* mod)
/*< Free the model parameter cubes >*/
{
  free(mod->vmod);
  free(mod->dmod);
}

void prepare_acquisition_2d( acq_struct_t* acq,
                          sf_axis axsou[2], sf_axis axrec[2], sf_axis axwav[2],
                          sf_file Fsou, sf_file Frec, sf_file Fwav)
/*< Read the acquisition geometry from files >*/
{
  // NUMBER OF SHOTS AND POSITIONS
  if (sf_n(axsou[1])!=2)
    sf_error("The source coordinate file is not correct!");

  // shot coordinates
  acq->ns = sf_n(axsou[0]);
  acq->xs = sf_floatalloc(acq->ns);
  acq->zs = sf_floatalloc(acq->ns);
  sf_floatread(acq->xs,acq->ns,Fsou);
  sf_floatread(acq->zs,acq->ns,Fsou);

  // receiver coordinates


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
  free(acq->xs);
  free(acq->zs);

  free(acq->wav);
}

void prepare_wfl_2d(wfl_struct_t *wfl,mod_struct_t *mod)
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
