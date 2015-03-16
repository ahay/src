/*
model.c
Igor Terentyev.
*/
/*============================================================================*/

#include "model.h"

/*----------------------------------------------------------------------------*/

void next_step(TIMESTEPINDEX * ts) {
  if (ts->iv == ts->niv-1) {
    ts->iv=0;
    ts->it++;
  }
  else {
    ts->iv++;
  }
}

int less_than(TIMESTEPINDEX t1, TIMESTEPINDEX t2) {
  if ((t1.it<t2.it) || ((t1.it == t2.it) && (t1.iv<t2.iv))) return 1;
  return 0;
}

ireal get_time(TIMESTEPINDEX ts) {
  if (ts.niv>0)  
    return ((ireal)(ts.it)+((ireal)(ts.iv))/((ireal)(ts.niv)))*ts.dt;
  else return REAL_ZERO;
}

/*int im_setnull(IMODEL *model)*/
int im_construct(IMODEL * model) {

  model->specs=NULL;

  /* WWS 22.02.08 removed zero of ndim, added null constr for grid */  
  init_default_grid(&(model->g));

  model->nnei = 0;
  /*  sten_setnull(&(model->sten)); */

  rd_a_setnull(&(model->ld_a));
  rd_a_setnull(&(model->ld_c));
  rd_a_setnull(&(model->ld_p));
  
  model->ld_s = model->ld_r = NULL;

  model->tsind.it=0;
  model->tsind.iv=0;
  model->tsind.niv=0;
  model->tsind.dt=REAL_ZERO;
  model->tsind.rhs=REAL_ZERO;

  IASN(model->nls,IPNT_0);
  IASN(model->nrs,IPNT_0);

  model->active.resize(0);

  model->sfg.resize(RDOM_MAX_NARR);
  for (int i=0;i<RDOM_MAX_NARR; i++) init_default_grid(&(model->sfg[i]));

  return 0;
}

/*----------------------------------------------------------------------------*/

int im_destroy(IMODEL *model, void (*destr)(void **))
{
  rd_a_destroy(&(model->ld_a));
  if (model->ld_s) userfree_(model->ld_s);
  destr(&(model->specs));
  /* since fd object allocates, should also free */
  /* if (model->specs) userfree_(model->specs); */
  return im_construct(model);
}
/*----------------------------------------------------------------------------*/

int im_setndim(IMODEL *model)
{
    int nnei, i, err;
    RDOM *p;

    err = gen_3n1(model->g.dim, &nnei);
    if ( err ) return err;

    p = (RDOM*)usermalloc_(nnei * 2 * sizeof(RDOM));
    if ( p == NULL ) return E_ALLOC;
    
    for ( i = 0; i < nnei * 2; ++ i ) rd_a_setnull(p + i);

    if ( model->ld_s != NULL ) userfree_(model->ld_s);
    model->ld_s = p;
    model->ld_r = model->ld_s + nnei;
    model->nnei = nnei;
    
    return 0;
}
/*----------------------------------------------------------------------------*/
