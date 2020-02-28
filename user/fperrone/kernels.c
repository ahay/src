#include <rsf.h>
#include "prep_utils.h"

static void velupd(wfl_struct_t* wfl, mod_struct_t const* mod, acq_struct_t const * acq, adj_t adjflag)
{

  long const n1 = wfl->simN1;
  long const n2 = wfl->simN2;

  float * const v1c = wfl->v1c;
  float * const v2c = wfl->v2c;
  float const * v1p = wfl->v1p;
  float const * v2p = wfl->v2p;

  float * const v1a = wfl->v1a;
  float * const v2a = wfl->v2a;

  float const * pp = wfl->pp;
  float * const pa = wfl->pa;

  float const * tap1 = wfl->tap1;
  float const * tap2 = wfl->tap2;

  float const * buoy = mod->buoy;
  float const * incomp = mod->incomp;
  float const dt = acq->dt;
  float const d1 = mod->d1;
  float const d2 = mod->d2;

  float const dtd1 = dt/d1;
  float const dtd2 = dt/d2;

  long i2start=NOP;
  long i2end = n2-NOP;
  long i1start=NOP;
  long i1end = n1-NOP;

  switch (adjflag){
  case FWD:

    for (long i2=i2start; i2<i2end; i2++){
      for (long i1=i1start,idx=IDX2D(i1,i2); i1<i1end; i1++,idx++){

        v1a[idx] = (C1*(pp[i1+1+i2*n1] - pp[i1  +i2*n1])+
                    C2*(pp[i1+2+i2*n1] - pp[i1-1+i2*n1])+
                    C3*(pp[i1+3+i2*n1] - pp[i1-2+i2*n1]))*dtd1;
      }
    }

    for (long i2=i2start; i2<i2end; i2++){
      for (long i1=i1start,idx=IDX2D(i1,i2); i1<i1end; i1++,idx++){

        v2a[idx] = (C1*(pp[i1+(i2+1)*n1] - pp[i1+ i2   *n1])+
                    C2*(pp[i1+(i2+2)*n1] - pp[i1+(i2-1)*n1])+
                    C3*(pp[i1+(i2+3)*n1] - pp[i1+(i2-2)*n1]))*dtd2;

      }
    }


    // 1-component
    for (long i2=i2start; i2<i2end; i2++){
      float const spox = tap2[i2];
      for (long i1=i1start,idx=IDX2D(i1,i2); i1<i1end; i1++,idx++){
        float const spo = spox*tap1[i1];
        float const b =buoy[idx];
        v1c[idx] = spo*(v1p[idx] - b*v1a[idx]);
      }
    }

    // 2-component
    for (long i2=i2start; i2<i2end; i2++){
      float const spox = tap2[i2];
      for (long i1=i1start,idx=IDX2D(i1,i2); i1<i1end; i1++,idx++){
        float const spo = spox*tap1[i1];
        float const b =buoy[idx];
        v2c[idx] = spo*(v2p[idx] - b*v2a[idx]);
      }
    }

    break;
  case ADJ:
    // ===============================================================
    // 2nd order in time
    for (int i2=i2start; i2<i2end; i2++){
      float const spo2 = tap2[i2];
      for (int i1=i1start, idx=IDX2D(i1,i2  ); i1<i1end; i1++,idx++){
        float const spo = tap1[i1]* spo2;
        float k = incomp[idx]*dt;
        pa[idx] = spo*k*pp[idx];
      }
    }

    for (long i2=i2start; i2<i2end; i2++){
      for (long i1=i1start, idx=IDX2D(i1,i2  ); i1<i1end; i1++,idx++){
        v1a[idx] = (C1*(pa[i1+1+i2*n1] - pa[i1  +i2*n1])+
                    C2*(pa[i1+2+i2*n1] - pa[i1-1+i2*n1])+
                    C3*(pa[i1+3+i2*n1] - pa[i1-2+i2*n1]))/d1;
      }
    }

    for (long i2=i2start; i2<i2end; i2++){
      for (long i1=i1start, idx=IDX2D(i1,i2  ); i1<i1end; i1++,idx++){
        v2a[idx] = (C1*(pa[i1+(i2+1)*n1] - pa[i1+i2    *n1])+
                    C2*(pa[i1+(i2+2)*n1] - pa[i1+(i2-1)*n1])+
                    C3*(pa[i1+(i2+3)*n1] - pa[i1+(i2-2)*n1]))/d2;
      }
    }

    // 1-component
    for (long i2=i2start; i2<i2end; i2++){
      float const spox = tap2[i2];
      for (long i1=i1start, idx=IDX2D(i1,i2  ); i1<i1end; i1++,idx++){
        float const spo = spox*tap1[i1];
        v1c[idx] = spo*v1p[idx] - v1a[idx];
      }
    }

    // 2-component
    for (long i2=i2start; i2<i2end; i2++){
      float const spox = tap2[i2];
      for (long i1=i1start, idx=IDX2D(i1,i2  ); i1<i1end; i1++,idx++){
        float const spo = spox*tap1[i1];
        v2c[idx] = spo*v2p[idx] - v2a[idx];
      }
    }

    break;
  }

}

static void presupd(wfl_struct_t* wfl, mod_struct_t const* mod, acq_struct_t const* acq, adj_t adjflag)
{

  long const n1 = wfl->simN1;
  long const n2 = wfl->simN2;

  float * const pc = wfl->pc;
  float const * pp = wfl->pp;

  float const *v1c = wfl->v1c;
  float const *v2c = wfl->v2c;
  float * const v1p = wfl->v1p;
  float * const v2p = wfl->v2p;

  float * const v1a = wfl->v1a;
  float * const v2a = wfl->v2a;

  float const * tap1 = wfl->tap1;
  float const * tap2 = wfl->tap2;

  float const *buoy   = mod->buoy;
  float const *incomp = mod->incomp;
  float const d1 = mod->d1;
  float const d2 = mod->d2;
  float const dt = acq->dt;

  float const dtd1 = dt/d1;
  float const dtd2 = dt/d2;

  long i2start=NOP;
  long i2end = n2-NOP;
  long i1start=NOP;
  long i1end = n1-NOP;

  switch (adjflag){
  case FWD:

    for (long i2=i2start; i2<i2end; i2++){
      for (long i1=i1start,idx=IDX2D(i1,i2); i1<i1end; i1++,idx++){
            v1a[idx] = -(C1*(v1c[i1  +i2*n1] - v1c[i1-1+i2*n1])+
                         C2*(v1c[i1+1+i2*n1] - v1c[i1-2+i2*n1])+
                         C3*(v1c[i1+2+i2*n1] - v1c[i1-3+i2*n1]))*dtd1;
      }
    }

    for (long i2=i2start; i2<i2end; i2++){
      for (long i1=i1start,idx=IDX2D(i1,i2); i1<i1end; i1++,idx++){
        v2a[idx] = -(C1*(v2c[i1+(i2  )*n1] - v2c[i1  +(i2-1)*n1])+
                     C2*(v2c[i1+(i2+1)*n1] - v2c[i1  +(i2-2)*n1])+
                     C3*(v2c[i1+(i2+2)*n1] - v2c[i1  +(i2-3)*n1]))*dtd2;
      }
    }

    for (long i2=i2start; i2<i2end; i2++){
      float const spox = tap2[i2];
      for (long i1=i2start,idx=IDX2D(i1,i2); i1<i1end; i1++,idx++){
        float const spo = spox*tap1[i1];
        float const k = incomp[idx];
        pc[idx] = spo*(pp[idx] + k*(v1a[idx]+v2a[idx]));
      }
    }

    break;
  case ADJ:
    // ===============================================================
    // 2nd order in time
    for (long i2=i2start; i2<i2end; i2++){
      float const spo2 = tap2[i2];
      for (long i1=i1start,idx=IDX2D(i1,i2); i1<i1end; i1++,idx++){
        float const spo = spo2*tap1[i1];
        float irho = buoy[idx]*dt;
        v2p[idx] = spo*irho*v2c[idx];
        v1p[idx] = spo*irho*v1c[idx];
      }
    }

    for (long i2=i2start; i2<i2end; i2++){
      for (long i1=i1start,idx=IDX2D(i1,i2); i1<i1end; i1++,idx++){
        v1a[idx] = -(C1*(v1p[i1  +i2*n1] - v1p[i1-1+i2*n1])+
                     C2*(v1p[i1+1+i2*n1] - v1p[i1-2+i2*n1])+
                     C3*(v1p[i1+2+i2*n1] - v1p[i1-3+i2*n1]))/d1;
      }
    }

    for (long i2=i2start; i2<i2end; i2++){
      for (long i1=i1start,idx=IDX2D(i1,i2); i1<i1end; i1++,idx++){
        v2a[idx] = -(C1*(v2p[i1+(i2  )*n1] - v2p[i1  +(i2-1)*n1])+
                     C2*(v2p[i1+(i2+1)*n1] - v2p[i1  +(i2-2)*n1])+
                     C3*(v2p[i1+(i2+2)*n1] - v2p[i1  +(i2-3)*n1]))/d2;
      }
    }

    for (long i2=i2start; i2<i2end; i2++){
      float const spox = tap2[i2];
      for (long i1=i1start,idx=IDX2D(i1,i2); i1<i1end; i1++,idx++){
        float const spo = spox*tap1[i1];
        pc[idx] = spo*pp[idx] + (v1a[idx]+v2a[idx]);
      }
    }

    break;
  }

}


static void injectPsource(wfl_struct_t* wfl,
                          mod_struct_t const * mod,
                          acq_struct_t const * acq,
                          long it)
{

  long nsou = acq->ns;

  long N1 = wfl->simN1;

  float modD1 = mod->d1;
  float modD2 = mod->d2;
  float o1 = wfl->simO1;
  float o2 = wfl->simO2;

  float dt = acq->dt;
  float scale = dt/(modD1*modD2);

  float * wf = wfl->pc;

  for (long isou=0; isou<nsou; isou++){
    float xs = acq->scoord[isou*2];
    float zs = acq->scoord[isou*2+1];

    int ixs = (xs-o2)/modD2;
    int izs = (zs-o1)/modD1;
    float force = acq->wav[isou + nsou*it]*scale;
    long idx = izs + N1*ixs;

    for (int j=-3,jh=0; j<=4; j++,jh++){
      const float hicks2 = acq->hicksSou2[jh+isou*8];
      for (int i=-3,ih=0; i<=4; i++,ih++){
        const float hc = acq->hicksSou1[ih+isou*8]*hicks2;
        wf[idx + i + N1*j] += hc*force;
      }
    }

  }

}

static void injectPdata(wfl_struct_t* wfl, mod_struct_t const * mod, acq_struct_t const * acq, long it){

  long nrec = acq->nr;

  long N1 = wfl->simN1;

  float modD1 = mod->d1;
  float modD2 = mod->d2;
  float o1 = wfl->simO1;
  float o2 = wfl->simO2;

  for (long irec=0; irec<nrec; irec++){
    float xr = acq->rcoord[irec*2];
    float zr = acq->rcoord[irec*2+1];

    int ixr = (xr-o2)/modD2;
    int izr = (zr-o1)/modD1;
    long idx = izr + N1*ixr;

    float force = acq->dat[irec + nrec*it];

    for (int j=-3,jh=0; j<=4; j++,jh++){
      const float hicks2 = acq->hicksRcv2[jh+irec*8];
      for (int i=-3,ih=0; i<=4; i++,ih++){
        const float hc = acq->hicksRcv1[ih+irec*8]*hicks2;
        wfl->pc[idx + i + N1*j] += hc*force;
      }
    }

  }

}

static void injectBornVelocitySource(wfl_struct_t * const wfl, mod_struct_t const *mod, acq_struct_t const * acq, long it){
  long modN1 = wfl->modN1;
  long modN2 = wfl->modN2;
  long n12 = modN1*modN2;
  long nb = wfl->nabc;

  float dt = acq->dt;

  fread(wfl->v1a,sizeof(float),n12,wfl->Fprgrd);
  fread(wfl->v2a,sizeof(float),n12,wfl->Fprgrd);

  for (long i2=0; i2<modN2; i2++){
    for (long i1=0; i1<modN1; i1++){
      long simIdx = (i1+nb) + (i2+nb)*wfl->simN1;
      long modIdx = i1 + i2*modN1;

      float dr = mod->denpert[modIdx];
      float rh = mod->dmod[modIdx];
      float rp = dr*dt/(rh*rh);

      wfl->v1c[simIdx] += rp*wfl->v1a[modIdx];
      wfl->v2c[simIdx] += rp*wfl->v2a[modIdx];
    }
  }

}

static void injectBornPressureSource(wfl_struct_t * const wfl, mod_struct_t const *mod, acq_struct_t const * acq, long it){

  long modN1 = wfl->modN1;
  long modN2 = wfl->modN2;
  long n12 = modN1*modN2;
  long nb = wfl->nabc;

  float dt = acq->dt;

  fread(wfl->bwfl,n12,sizeof(float),wfl->Fpvdiv);

  for (long i2=0; i2<modN2; i2++){
    for (long i1=0; i1<modN1; i1++){
      long simIdx = (i1+nb) + (i2+nb)*wfl->simN1;
      long modIdx = i1 + i2*modN1;

      float const dv = mod->velpert[modIdx];
      float const dr = mod->denpert[modIdx];
      float const vc = mod->vmod[modIdx];
      float const rh = mod->dmod[modIdx];
      float const vp= 2.f*vc*dv*rh;
      float const rp= vc*vc*dr;

      wfl->pc[simIdx] += (vp+rp)*dt*wfl->bwfl[modIdx];

    }
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


static void extract_vel_wfl_2d(wfl_struct_t* wfl, int comp)
{
  long modN1 = wfl->modN1;
  long modN2 = wfl->modN2;
  long nabc  = wfl->nabc;
  long simN1 = wfl->simN1;

  // copy the write chunk of wavefield
  switch (comp){
  case 1:
    for (long i2=0; i2<modN2; i2++)
      memcpy(wfl->bwfl+i2*modN1,wfl->v1c+(nabc+(i2+nabc)*simN1),modN1*sizeof(float));
    break;
  case 2:
    for (long i2=0; i2<modN2; i2++)
      memcpy(wfl->bwfl+i2*modN1,wfl->v2c+(nabc+(i2+nabc)*simN1),modN1*sizeof(float));
    break;
  }

}

static void extract_pres_wfl_2d(wfl_struct_t * const wfl){
  long modN1 = wfl->modN1;
  long modN2 = wfl->modN2;
  long nabc  = wfl->nabc;
  long simN1 = wfl->simN1;

  // copy the write chunk of wavefield
  for (long i2=0; i2<modN2; i2++)
    memcpy(wfl->bwfl+i2*modN1,wfl->pc+(nabc+(i2+nabc)*simN1),modN1*sizeof(float));

}

static void extract_dat_2d(wfl_struct_t* wfl,acq_struct_t const * acq){

  long nr = acq->nr;
  long n1 = wfl->simN1;
  float o1 = wfl->simO1;
  float o2 = wfl->simO2;
  float d1 = wfl->d1;
  float d2 = wfl->d2;

  wfl->rdata = sf_floatalloc(nr);
  for (long ir=0; ir<nr; ir++){
    float xr = acq->rcoord[2*ir];
    float zr = acq->rcoord[2*ir+1];
    long ixr = (xr - o2)/d2;
    long izr = (zr - o1)/d1;
    long idx = izr + n1*ixr;
    float rv = 0.;
    for (int j=-3,jh=0; j<=4; j++,jh++){
      const float hicks2 = acq->hicksRcv2[jh+ir*8];
      for (int i=-3,ih=0; i<=4; i++,ih++){
        const float hc = acq->hicksRcv1[ih+ir*8]*hicks2;
        rv += hc*wfl->pc[idx + i +j*n1];
      }
    }
    wfl->rdata[ir] = rv;
  }

  sf_floatwrite(wfl->rdata,nr,wfl->Fdata);

  free(wfl->rdata);

}

static void extract_scat_dat_2d(wfl_struct_t * const wfl,acq_struct_t const *acq){

  long nr = acq->nr;
  long n1 = wfl->simN1;
  float o1 = wfl->simO1;
  float o2 = wfl->simO2;
  float d1 = wfl->d1;
  float d2 = wfl->d2;

  wfl->rdata = sf_floatalloc(nr);
  for (long ir=0; ir<nr; ir++){
    float xr = acq->rcoord[2*ir];
    float zr = acq->rcoord[2*ir+1];
    long ixr = (xr - o2)/d2;
    long izr = (zr - o1)/d1;
    long idx = izr + n1*ixr;
    float rv = 0.;
    for (int j=-3,jh=0; j<=4; j++,jh++){
      const float hicks2 = acq->hicksRcv2[jh+ir*8];
      for (int i=-3,ih=0; i<=4; i++,ih++){
        const float hc = acq->hicksRcv1[ih+ir*8]*hicks2;
        rv += hc*wfl->pc[idx + i +j*n1];
      }
    }
    wfl->rdata[ir] = rv;
  }

  sf_floatwrite(wfl->rdata,nr,wfl->Fsdata);

  free(wfl->rdata);

}

static void applyFreeSurfaceBC(wfl_struct_t *wfl){

  long const nb = wfl->nabc;
  long const n1 = wfl->simN1;
  long const n2 = wfl->simN2;

  for (long i2=0; i2<n2; i2++)
    for (long i1=0; i1<nb+1; i1++)
      wfl->pc[i1+i2*n1]=0.;

}

void fwdextrap2d(wfl_struct_t * wfl, acq_struct_t const * acq, mod_struct_t const * mod)
/*< extrapolation kernel 2d - forward operator >*/
{
  sf_warning("FORWARD EXTRAPOLATION..");

  long modN1=wfl->modN1;
  long modN2=wfl->modN2;
  long nelem=modN1*modN2;
  long nt = acq->ntdat;

  // loop over time
  for (long it=0; it<nt; it++){
    bool save = ((wfl->snap) && !(it%wfl->jsnap));
    if ((it+1)%100==0)
      sf_warning("Time step %4d of %4d (%2.1f %%)",it, nt, ((100.*it)/nt));

    velupd(wfl,mod,acq,FWD);
    presupd(wfl,mod,acq,FWD);
    injectPsource(wfl,mod,acq,it);

    if (wfl->freesurf)
      applyFreeSurfaceBC(wfl);

    // write the wavefield out
    if (save){
      extract_pres_wfl_2d(wfl);
      sf_floatwrite(wfl->bwfl,nelem,wfl->Fwfl);
    }
    // extract the data at the receiver locations
    extract_dat_2d(wfl,acq);

    swapwfl(wfl);
  }

}

void bornbckwfl2d(wfl_struct_t * wfl, acq_struct_t const * acq,  mod_struct_t const * mod, born_setup_struct_t para)
/*< Born background wavefield extrapolation >*/
{
  long modN1 =wfl->modN1;
  long modN2 =wfl->modN2;
  int nt = acq->ntdat;
  long nelem = modN1*modN2;
  bool saveData= para.outputBackgroundData;

  // loop over time
  for (int it=0; it<nt; it++){

    velupd(wfl,mod,acq,FWD);
    presupd(wfl,mod,acq,FWD);
    injectPsource(wfl,mod,acq,it);

    if (wfl->freesurf)
      applyFreeSurfaceBC(wfl);

    // write the wavefield out
    extract_pres_wfl_2d(wfl);

    if (para.outputBackgroundWfl)
      sf_floatwrite(wfl->bwfl,nelem,wfl->Fwfl);
    else
      fwrite(wfl->bwfl,sizeof(float),nelem,para.Fbwfl);

    // extract the data at the receiver locations
    if (saveData) extract_dat_2d(wfl,acq);

    swapwfl(wfl);
  }

}

void bornfwdextrap2d(wfl_struct_t * wfl, acq_struct_t const * acq, mod_struct_t const * mod, born_setup_struct_t para)
/*< kernel for Born forward extrapolation >*/
{
  long modN1 = wfl->modN1;
  long modN2 = wfl->modN2;
  long nelem = modN1*modN2;
  int nt = acq->ntdat;

  // loop over time
  for (int it=0; it<nt; it++){
    bool save = (wfl->Fswfl);

    velupd(wfl,mod,acq,FWD);
    if (para.inputDenPerturbation)
      injectBornVelocitySource(wfl,mod,acq,it);
    presupd(wfl,mod,acq,FWD);
    if (para.inputVelPerturbation)
      injectBornPressureSource(wfl,mod,acq,it);

    if (wfl->freesurf)
      applyFreeSurfaceBC(wfl);

    // write the wavefield out
    if (save){
      extract_pres_wfl_2d(wfl);
      sf_floatwrite(wfl->bwfl,nelem,wfl->Fswfl);
    }

    // extract the data at the receiver locations
    extract_scat_dat_2d(wfl,acq);

    swapwfl(wfl);
  }
}

void adjextrap2d(wfl_struct_t * wfl, acq_struct_t const * acq, mod_struct_t const * mod)
/*< extrapolation kernel 2d - adjoint operator  >*/
{
  sf_warning("ADJOINT EXTRAPOLATION..");

  long modN1=wfl->modN1;
  long modN2=wfl->modN2;
  long nelem=modN1*modN2;
  long nt = acq->ntdat;

  // loop over time
  for (long it=0; it<nt; it++){
    bool save = ((wfl->snap) && !(it%wfl->jsnap));
    if ((it+1)%100==0)
      sf_warning("Time step %4d of %4d (%2.1f %%)",it, nt, ((100.*it)/nt));

    velupd(wfl,mod,acq,ADJ);
    presupd(wfl,mod,acq,ADJ);
    injectPsource(wfl,mod,acq,it);

    if (wfl->freesurf)
      applyFreeSurfaceBC(wfl);

    // write the wavefield out
    if (save) {
      extract_pres_wfl_2d(wfl);
      sf_floatwrite(wfl->bwfl,nelem,wfl->Fwfl);
    }

    swapwfl(wfl);
  }

}

void bornadjextrap2d(wfl_struct_t * wfl,
                     acq_struct_t const * acq,
                     mod_struct_t const * mod,
                     born_setup_struct_t *para)
/*< kernel for Born forward extrapolation >*/
{
  long modN1 = wfl->modN1;
  long modN2 = wfl->modN2;
  long nelem = modN1*modN2;
  int nt = acq->ntdat;
  bool save = para->outputScatteredWfl;
  // loop over time
  for (int it=0; it<nt; it++){

    velupd(wfl,mod,acq,ADJ);
    if (para->outputDenPertImage){
      extract_vel_wfl_2d(wfl,1);
      fwrite(wfl->bwfl,sizeof(float),nelem,para->Fpv1);
      extract_vel_wfl_2d(wfl,2);
      fwrite(wfl->bwfl,sizeof(float),nelem,para->Fpv2);
    }
    presupd(wfl,mod,acq,ADJ);
    injectPdata(wfl,mod,acq,it);

    if (wfl->freesurf)
      applyFreeSurfaceBC(wfl);

    // write the wavefield out
    extract_pres_wfl_2d(wfl);
    if (save)
      sf_floatwrite(wfl->bwfl,nelem,wfl->Fswfl);
    else
      fwrite(wfl->bwfl,sizeof(float),nelem,para->Fswfl);

    swapwfl(wfl);
  }

}

void setupABC(wfl_struct_t* wfl)
/*< Setup of the coefficients for the absorbing boundary >*/
{

  float tapbase = 0.92;
  float alpha = sqrt(-log(tapbase));

  for (long i=0; i<wfl->simN1; i++)
    wfl->tap1[i] = 1.;
  for (long i=0; i<wfl->simN2; i++)
    wfl->tap2[i] = 1.;

  for (int i=0; i<NOP; i++){
    float tap = exp(-alpha*alpha);
    wfl->tap1[i] = tap;
    wfl->tap2[i] = tap;
    wfl->tap1[wfl->simN1-1-i] = tap;
    wfl->tap2[wfl->simN2-1-i] = tap;
  }

  float taplen1 = wfl->nabc-NOP;
  float taplen2 = wfl->nabc-NOP;
  for (int i=0,j=NOP; i<taplen1; i++,j++){
    float arg = alpha*(taplen1-i)/taplen1;
    float tap = exp(-arg*arg);
    wfl->tap1[j] = tap;
    wfl->tap1[wfl->simN1-1-j] = tap;
  }
  for (int i=0,j=NOP; i<taplen2; i++,j++){
    float arg = alpha*(taplen2-i)/taplen2;
    float tap = exp(-arg*arg);
    wfl->tap2[j] = tap;
    wfl->tap2[wfl->simN2-1-j] = tap;
  }

}

void reset_wfl(wfl_struct_t* wfl)
/*< reset the wavefields to zero >*/
{

  long n1=wfl->simN1;
  long n2=wfl->simN2;
  memset(wfl->v1c,0,n1*n2*sizeof(float));
  memset(wfl->v2c,0,n1*n2*sizeof(float));
  memset(wfl->v1p,0,n1*n2*sizeof(float));
  memset(wfl->v2p,0,n1*n2*sizeof(float));

  memset(wfl->pc,0,n1*n2*sizeof(float));
  memset(wfl->pp,0,n1*n2*sizeof(float));

}

