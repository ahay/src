#include <rsf.h>

#ifndef _PREP_H

typedef struct wfl_struct wfl_struct_t;
/*^*/

typedef struct acq_struct acq_struct_t;
/*^*/

typedef struct mod_struct mod_struct_t;
/*^*/

struct wfl_struct{
  float *pc;
  float *pp;
  float *pa;
  int n1;
  int n2;
  float d1;
  float d2;
  float o1;
  float o2;
};
/*^*/

struct acq_struct{
  int ns;
  float *nr;
  float *xs;
  float *ys;
  float *zs;
  float *xr;
  float *yr;
  float *zr;
};
/*^*/

struct mod_struct{
  int n1;
  int n2;
  float d1;
  float d2;
  float o1;
  float o2;
  float *vmod;
  float *dmod;
};
/*^*/

#endif

void prepare_model_2d(mod_struct_t* mod, sf_axis ax[2],sf_file Fmod, char const * model)
/*< Sets up the specified model cube >*/
{

  if (!strcmp(model,"VELOCITY")){
    int n1 = sf_n(ax[0]);
    int n2 = sf_n(ax[1]);
    mod->vmod = (float*) sf_alloc(n1*n2,sizeof(float));
    sf_floatread(mod->vmod,n1*n2,Fmod);

    double vave = 0;
    for (int i=0; i<n1*n2; i++)
      vave += mod->vmod[i];
    vave /= (n1*n2);
    sf_warning("Velocity Model average value = %g",vave);
  }
  if (!strcmp(model,"DENSITY")){
    int n1 = sf_n(ax[0]);
    int n2 = sf_n(ax[1]);
    mod->dmod = (float*) sf_alloc(n1*n2,sizeof(float));
    sf_floatread(mod->dmod,n1*n2,Fmod);

    double dave = 0;
    for (int i=0; i<n1*n2; i++)
      dave += mod->dmod[i];
    dave /= (n1*n2);
    sf_warning("Velocity Model average value = %g",dave);
  }

}

void clear_model_2d(mod_struct_t* mod)
/*< Free the model parameter cubes >*/
{
  free(mod->vmod);
  free(mod->dmod);
}

void prepare_acquisition(acq_struct_t* acq, sf_axis axsou[2], sf_axis axrec[2], sf_file Fsou, sf_file Frec)
/*< Read the acquisition geometry from files >*/
{
  // NUMBER OF SHOTS AND POSITIONS
  if (sf_n(axsou[1])!=2)
    sf_error("The source coordinate file is not correct!");

  acq->ns = sf_n(axsou[0]);
  acq->xs = sf_floatalloc(acq->ns);
  acq->zs = sf_floatalloc(acq->ns);
  sf_floatread(acq->xs,acq->ns,Fsou);
  sf_floatread(acq->zs,acq->ns,Fsou);

}
