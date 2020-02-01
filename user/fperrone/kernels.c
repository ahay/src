#include <rsf.h>
#include "prep_utils.h"

static void velupd(wfl_struct_t* wfl, mod_struct_t const* mod, acq_struct_t const * acq, bool adjflag)
{

  long n1 = wfl->simN1;
  long n2 = wfl->simN2;

  float *v1c = wfl->v1c;
  float *v2c = wfl->v2c;
  float *v1p = wfl->v1p;
  float *v2p = wfl->v2p;

  float* pp = wfl->pp;

  float* buoy = mod->buoy;
  float dt = acq->dt;
  float d1 = mod->d1;
  float d2 = mod->d2;

  if (!adjflag){

    // 1-component
    for (long i2=0; i2<n2; i2++){
      for (long i1=0; i1<n1-1; i1++){
        float const pd1 = (pp[i1+1+i2*n1] - pp[i1+i2*n1])/d1;
        v1c[i1+i2*n1] = v1p[i1+i2*n1] - buoy[i1+i2*n1]*pd1*dt;
      }
    }

    // 2-component
    for (long i2=0; i2<n2-1; i2++){
      for (long i1=0; i1<n1; i1++){
        float const pd2 = (pp[i1+(i2+1)*n1] - pp[i1+i2*n1])/d2;
        v2c[i1+i2*n1] = v2p[i1+i2*n1] - buoy[i1+i2*n1]*pd2*dt;
      }
    }

  }
  else{

  }

}

static void presupd(wfl_struct_t* wfl, mod_struct_t const* mod, acq_struct_t const* acq, bool adjflag)
{

  long n1 = wfl->simN1;
  long n2 = wfl->simN2;

  float* pc = wfl->pc;
  float* pp = wfl->pp;

  float *v1c = wfl->v1c;
  float *v2c = wfl->v2c;

  float *incomp = mod->incomp;
  float d1 = mod->d1;
  float d2 = mod->d2;
  float dt = acq->dt;

  if (!adjflag){

    for (long i2=1; i2<n2; i2++){
      for (long i1=1; i1<n1; i1++){
        float const v1d1 = -(v1c[i1+i2*n1] - v1c[i1-1+i2*n1])/d1;
        float const v2d2 = -(v2c[i1+i2*n1] - v2c[i1  +(i2-1)*n1])/d2;

        pc[i1+i2*n1] = pp[i1+i2*n1] + incomp[i1+i2*n1]*(v1d1+v2d2)*dt;
      }
    }

  }
  else{

  }

}


static void injectPsource(wfl_struct_t* wfl, mod_struct_t const * mod, acq_struct_t const * acq, long it)
{

  long nsou = acq->ns;
  long nt   = acq->nt;

  long modN1 = mod->n1;

  float modD1 = mod->d1;
  float modD2 = mod->d2;
  float modO1 = mod->o1;
  float modO2 = mod->o2;

  for (long isou=0; isou<nsou; isou++){
    float xs = acq->scoord[isou*2];
    float zs = acq->scoord[isou*2+1];

    int ixs = (xs-modO2)/modD2;
    int izs = (zs-modO1)/modD1;

    wfl->pc[izs + modN1*ixs] += acq->wav[it+isou*nt];
  }

}

static void swapwfl(wfl_struct_t* wfl)
{

  float *tmp = wfl->pc;
  wfl->pc = wfl->pp;
  wfl->pp = tmp;

  tmp = wfl->v1c;
  wfl->v1c = wfl->v1p;
  wfl->v1p = tmp;

  tmp = wfl->v2c;
  wfl->v2c = wfl->v2p;
  wfl->v2p = tmp;

}


void fwdextrap2d(wfl_struct_t * wfl, acq_struct_t const * acq, mod_struct_t const * mod)
/*< extrapolation kernel 2d - forward operator >*/
{
  sf_warning("FORWARD EXTRAPOLATION..");

  int nt = acq->nt;
  long nelem = wfl->simN1*wfl->simN2;

  // loop over time
  for (int it=0; it<nt; it++){
    velupd(wfl,mod,acq,false);
    presupd(wfl,mod,acq,false);
    injectPsource(wfl,mod,acq,it);

    // write the wavefield out
    sf_floatwrite(wfl->pc,nelem,wfl->Fwfl);

    swapwfl(wfl);
  }

}

void adjextrap2d(wfl_struct_t * wfl, acq_struct_t const * acq, mod_struct_t const * mod)
/*< extrapolation kernel 2d - adjoint operator  >*/
{
  sf_warning("ADJOINT EXTRAPOLATION..");

  int nt = acq->nt;
  long nelem = wfl->simN1*wfl->simN2;

  // loop over time
  for (int it=0; it<nt; it++){
    velupd(wfl,mod,acq,true);
    presupd(wfl,mod,acq,true);
    injectPsource(wfl,mod,acq,it);

    // write the wavefield out
    sf_floatwrite(wfl->pc,nelem,wfl->Fwfl);

    swapwfl(wfl);
  }

}


