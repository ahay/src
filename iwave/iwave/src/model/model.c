/* 
model.c
Igor Terentyev.
*/
/*============================================================================*/

#include "utils.h"
#include "fdfun.h"

/*----------------------------------------------------------------------------*/

void init_step(TIMESTEPINDEX * ts, int it) {
  ts->it=it;
  ts->iv=0;
}

void next_step(TIMESTEPINDEX * ts, TIMESTEPINFO tsinfo) {
  if (ts->iv == tsinfo.npair-1) {
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
  return ts.it*ts.dt;
}

/*int im_setnull(IMODEL *model)*/
int im_construct(IMODEL * model) {

  model->specs=NULL;

  /* WWS 22.02.08 removed zero of ndim, added null constr for grid */  
  init_default_grid(&(model->g));

  model->nnei = 0;
  //  sten_setnull(&(model->sten));

  model->tsinfo.narr = 0;
  rd_a_setnull(&(model->ld_a));
  rd_a_setnull(&(model->ld_c));
  rd_a_setnull(&(model->ld_p));
  
  model->ld_s = model->ld_r = NULL;

  init_step(&(model->tsind),INT_MAX);

  return 0;
}

/*----------------------------------------------------------------------------*/

int im_destroy(IMODEL *model)
{
  rd_a_destroy(&(model->ld_a));
  free(model->ld_s);
  fd_modeldest(model);
  free(model->specs);
  return im_construct(model);
}
/*----------------------------------------------------------------------------*/

int im_setndim(IMODEL *model)
{
    int nnei, i, err;
    RDOM *p;

    err = gen_3n1(model->g.dim, &nnei);
    if ( err ) return err;

    p = (RDOM*)malloc(nnei * 2 * sizeof(RDOM));
    if ( p == NULL ) return E_ALLOC;
    
    for ( i = 0; i < nnei * 2; ++ i ) rd_a_setnull(p + i);

    if ( model->ld_s != NULL ) free(model->ld_s);
    model->ld_s = p;
    model->ld_r = model->ld_s + nnei;
    model->nnei = nnei;
    
    return 0;
}
/*----------------------------------------------------------------------------*/
