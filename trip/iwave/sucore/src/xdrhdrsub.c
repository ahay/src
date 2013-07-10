/*
 * xdrhdrsub.c - subroutine for segy header R/W
 * THIS FILE IS GENERATED AUTOMATICALLY - 
 * see the makefile in this directory
 */

#ifdef SUXDR
#include "su_xdr.h"
int xdrhdrsub(XDR *segy_xdr, segy *trace)
{
 int i1, i2;
 int status;
 cwp_Bool outgoing;
 unsigned int u1, u2, utemp;
#define HHTOU(H1,H2,U) u1=H1;u1&=65535;u1<<=16;\
                       u2=H2;u2&=65535;U=u1|u2;
#define UTOHH(U,H1,H2) \
        u2=(U);i2=u2&65535;H2=(i2>32767)?(i2-65536):i2;\
        u1=(U)>>16;i1=u1&65535;H1=(i1>32767)?(i1-65536):i1;

#define HUTOU(H1,U2,U) u1=H1;u1&=65535;u1<<=16;\
                       u2=U2;u2&=65535;U=u1|u2;
#define UTOHU(U,H1,U2) \
        u2=(U);U2=u2&65535;\
        u1=(U)>>16;i1=u1&65535;H1=(i1>32767)?(i1-65536):i1;

#define UHTOU(U1,H2,U) u1=U1;u1&=65535;u1<<=16;\
                       u2=H2;u2&=65535;U=u1|u2;
#define UTOUH(U,U1,H2) \
        u2=(U);i2=u2&65535;H2=(i2>32767)?(i2-65536):i2;\
        u1=(U)>>16;U1=u1&65535;

#define UUTOU(U1,U2,U) u1=U1;u1&=65535;u1<<=16;\
                       u2=U2;u2&=65535;U=u1|u2;
#define UTOUU(U,U1,U2) \
        u2=(U);U2=u2&65535;\
        u1=(U)>>16;U1=u1&65535;


 status=TRUE;
 outgoing = (segy_xdr->x_op == XDR_ENCODE)?cwp_true:cwp_false;
 if(FALSE == xdr_int(segy_xdr,(&trace->tracl))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->tracr))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->fldr))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->tracf))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->ep))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->cdp))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->cdpt))) return(FALSE);
 if(outgoing) {
    HHTOU(trace->trid,trace->nvs,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->trid,trace->nvs)
    }
 if(outgoing) {
    HHTOU(trace->nhs,trace->duse,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->nhs,trace->duse)
    }
 if(FALSE == xdr_int(segy_xdr,(&trace->offset))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->gelev))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->selev))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->sdepth))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->gdel))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->sdel))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->swdep))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->gwdep))) return(FALSE);
 if(outgoing) {
    HHTOU(trace->scalel,trace->scalco,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->scalel,trace->scalco)
    }
 if(FALSE == xdr_int(segy_xdr,(&trace->sx))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->sy))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->gx))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->gy))) return(FALSE);
 if(outgoing) {
    HHTOU(trace->counit,trace->wevel,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->counit,trace->wevel)
    }
 if(outgoing) {
    HHTOU(trace->swevel,trace->sut,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->swevel,trace->sut)
    }
 if(outgoing) {
    HHTOU(trace->gut,trace->sstat,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->gut,trace->sstat)
    }
 if(outgoing) {
    HHTOU(trace->gstat,trace->tstat,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->gstat,trace->tstat)
    }
 if(outgoing) {
    HHTOU(trace->laga,trace->lagb,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->laga,trace->lagb)
    }
 if(outgoing) {
    HHTOU(trace->delrt,trace->muts,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->delrt,trace->muts)
    }
 if(outgoing) {
    HUTOU(trace->mute,trace->ns,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHU(utemp,trace->mute,trace->ns)
    }
 if(outgoing) {
    UHTOU(trace->dt,trace->gain,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOUH(utemp,trace->dt,trace->gain)
    }
 if(outgoing) {
    HHTOU(trace->igc,trace->igi,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->igc,trace->igi)
    }
 if(outgoing) {
    HHTOU(trace->corr,trace->sfs,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->corr,trace->sfs)
    }
 if(outgoing) {
    HHTOU(trace->sfe,trace->slen,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->sfe,trace->slen)
    }
 if(outgoing) {
    HHTOU(trace->styp,trace->stas,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->styp,trace->stas)
    }
 if(outgoing) {
    HHTOU(trace->stae,trace->tatyp,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->stae,trace->tatyp)
    }
 if(outgoing) {
    HHTOU(trace->afilf,trace->afils,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->afilf,trace->afils)
    }
 if(outgoing) {
    HHTOU(trace->nofilf,trace->nofils,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->nofilf,trace->nofils)
    }
 if(outgoing) {
    HHTOU(trace->lcf,trace->hcf,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->lcf,trace->hcf)
    }
 if(outgoing) {
    HHTOU(trace->lcs,trace->hcs,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->lcs,trace->hcs)
    }
 if(outgoing) {
    HHTOU(trace->year,trace->day,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->year,trace->day)
    }
 if(outgoing) {
    HHTOU(trace->hour,trace->minute,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->hour,trace->minute)
    }
 if(outgoing) {
    HHTOU(trace->sec,trace->timbas,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->sec,trace->timbas)
    }
 if(outgoing) {
    HHTOU(trace->trwf,trace->grnors,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->trwf,trace->grnors)
    }
 if(outgoing) {
    HHTOU(trace->grnofr,trace->grnlof,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->grnofr,trace->grnlof)
    }
 if(outgoing) {
    HHTOU(trace->gaps,trace->otrav,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->gaps,trace->otrav)
    }
 if(FALSE == xdr_float(segy_xdr,&(trace->d1))) return(FALSE);
 if(FALSE == xdr_float(segy_xdr,&(trace->f1))) return(FALSE);
 if(FALSE == xdr_float(segy_xdr,&(trace->d2))) return(FALSE);
 if(FALSE == xdr_float(segy_xdr,&(trace->f2))) return(FALSE);
 if(FALSE == xdr_float(segy_xdr,&(trace->ungpow))) return(FALSE);
 if(FALSE == xdr_float(segy_xdr,&(trace->unscale))) return(FALSE);
 if(FALSE == xdr_int(segy_xdr,(&trace->ntr))) return(FALSE);
 if(outgoing) {
    HHTOU(trace->mark,trace->shortpad,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->mark,trace->shortpad)
    }
 if(outgoing) {
    HHTOU(trace->unass[0],trace->unass[1],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(trace->unass[2],trace->unass[3],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(trace->unass[4],trace->unass[5],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(trace->unass[6],trace->unass[7],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(trace->unass[8],trace->unass[9],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(trace->unass[10],trace->unass[11],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(trace->unass[12],trace->unass[13],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->unass[0],trace->unass[1])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->unass[2],trace->unass[3])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->unass[4],trace->unass[5])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->unass[6],trace->unass[7])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->unass[8],trace->unass[9])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->unass[10],trace->unass[11])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,trace->unass[12],trace->unass[13])
    }

 return(status);
}
#else
void xdrhdrsub(){return;}
#endif
