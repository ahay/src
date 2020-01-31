#include <rsf.h>
#include "prep_utils.h"

static void velupd(wfl_struct_t* wfl, bool adjflag)
{

  long n1 = wfl->n1;
  long n2 = wfl->n2;

  float *v1c = wfl->v1c;
  float *v2c = wfl->v2c;
  float *v1p = wfl->v1p;
  float *v2p = wfl->v2p;

  if (!adjflag){

    for (long i2=0; i2<n2; i2++){
      for (long i1=0; i1<n1; i1++){
        v1c[i1+i2*n1] = v1p[i1+i2*n1];
        v2c[i1+i2*n1] = v2p[i1+i2*n1];
      }
    }

  }
  else{

    for (long i2=0; i2<n2; i2++){
      for (long i1=0; i1<n1; i1++){
        v1c[i1+i2*n1] = v1p[i1+i2*n1];
        v2c[i1+i2*n1] = v2p[i1+i2*n1];
      }
    }

  }

}

static void presupd(wfl_struct_t* wfl, bool adjflag)
{

  if (!adjflag){

  }
  else{

  }

}

static void swapwfl(wfl_struct_t* wfl)
{

}

static void resetwfl(wfl_struct_t* wfl)
{

}


void fwdextrap2d(wfl_struct_t * wfl, acq_struct_t const * acq, mod_struct_t const * mod)
/*< extrapolation kernel 2d - forward operator >*/
{
  sf_warning("FORWARD EXTRAPOLATION..");

  long nshots = acq->ns;
  float* xs = acq->xs;
  float* zs = acq->zs;
  sf_warning("Number of shots to model : %d",nshots);

  int nt = acq->nt;

  // loop over shots
  for (long ishot=0; ishot<acq->ns; ishot++){
    sf_warning("Shot %d",ishot+1);
    sf_warning("xs[%d]=%g",ishot+1,xs[ishot]);
    sf_warning("zs[%d]=%g",ishot+1,zs[ishot]);

    // loop over time
    for (int it=0; it<nt; it++){
      velupd(wfl,false);
      presupd(wfl,false);
      swapwfl(wfl);
    }

    resetwfl(wfl);
  }

}

void adjextrap2d(wfl_struct_t * wfl, acq_struct_t const * acq, mod_struct_t const * mod)
/*< extrapolation kernel 2d - adjoint operator  >*/
{
  sf_warning("ADJOINT EXTRAPOLATION..");

  long nshots = acq->ns;
  float* xs = acq->xs;
  float* zs = acq->zs;
  sf_warning("Number of shots to model : %d",nshots);

  int nt = acq->nt;

  // loop over shots
  for (int ishot=0; ishot<acq->ns; ishot++){
    sf_warning("Shot %d",ishot+1);
    sf_warning("xs[%d]=%g",ishot+1,xs[ishot]);
    sf_warning("zs[%d]=%g",ishot+1,zs[ishot]);

    // loop over time
    for (int it=0; it<nt; it++){
      velupd(wfl,true);
      presupd(wfl,true);
      swapwfl(wfl);
    }

    resetwfl(wfl);
  }

}


