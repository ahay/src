/*
 * xdrbhdrsub.c - subroutine for segy header R/W
 * THIS FILE IS GENERATED AUTOMATICALLY - 
 * see the makefile in this directory
 */

#ifdef SUXDR
#include "su_xdr.h"
int xdrbhdrsub(XDR *segy_xdr, bhed *binhdr)
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
 if(FALSE == xdr_int(segy_xdr,(&binhdr->jobid))) status = FALSE;
 if(FALSE == xdr_int(segy_xdr,(&binhdr->lino))) status = FALSE;
 if(FALSE == xdr_int(segy_xdr,(&binhdr->reno))) status = FALSE;
 if(outgoing) {
    HHTOU(binhdr->ntrpr,binhdr->nart,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->ntrpr,binhdr->nart)
    }
 if(outgoing) {
    HHTOU(binhdr->hdt,binhdr->dto,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hdt,binhdr->dto)
    }
 if(outgoing) {
    HHTOU(binhdr->hns,binhdr->nso,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hns,binhdr->nso)
    }
 if(outgoing) {
    HHTOU(binhdr->format,binhdr->fold,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->format,binhdr->fold)
    }
 if(outgoing) {
    HHTOU(binhdr->tsort,binhdr->vscode,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->tsort,binhdr->vscode)
    }
 if(outgoing) {
    HHTOU(binhdr->hsfs,binhdr->hsfe,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hsfs,binhdr->hsfe)
    }
 if(outgoing) {
    HHTOU(binhdr->hslen,binhdr->hstyp,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hslen,binhdr->hstyp)
    }
 if(outgoing) {
    HHTOU(binhdr->schn,binhdr->hstas,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->schn,binhdr->hstas)
    }
 if(outgoing) {
    HHTOU(binhdr->hstae,binhdr->htatyp,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hstae,binhdr->htatyp)
    }
 if(outgoing) {
    HHTOU(binhdr->hcorr,binhdr->bgrcv,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hcorr,binhdr->bgrcv)
    }
 if(outgoing) {
    HHTOU(binhdr->rcvm,binhdr->mfeet,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->rcvm,binhdr->mfeet)
    }
 if(outgoing) {
    HHTOU(binhdr->polyt,binhdr->vpol,utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->polyt,binhdr->vpol)
    }
 if(outgoing) {
    HHTOU(binhdr->hunass[0],binhdr->hunass[1],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[2],binhdr->hunass[3],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[4],binhdr->hunass[5],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[6],binhdr->hunass[7],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[8],binhdr->hunass[9],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[10],binhdr->hunass[11],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[12],binhdr->hunass[13],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[14],binhdr->hunass[15],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[16],binhdr->hunass[17],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[18],binhdr->hunass[19],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[20],binhdr->hunass[21],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[22],binhdr->hunass[23],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[24],binhdr->hunass[25],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[26],binhdr->hunass[27],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[28],binhdr->hunass[29],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[30],binhdr->hunass[31],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[32],binhdr->hunass[33],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[34],binhdr->hunass[35],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[36],binhdr->hunass[37],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[38],binhdr->hunass[39],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[40],binhdr->hunass[41],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[42],binhdr->hunass[43],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[44],binhdr->hunass[45],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[46],binhdr->hunass[47],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[48],binhdr->hunass[49],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[50],binhdr->hunass[51],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[52],binhdr->hunass[53],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[54],binhdr->hunass[55],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[56],binhdr->hunass[57],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[58],binhdr->hunass[59],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[60],binhdr->hunass[61],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[62],binhdr->hunass[63],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[64],binhdr->hunass[65],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[66],binhdr->hunass[67],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[68],binhdr->hunass[69],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[70],binhdr->hunass[71],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[72],binhdr->hunass[73],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[74],binhdr->hunass[75],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[76],binhdr->hunass[77],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[78],binhdr->hunass[79],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[80],binhdr->hunass[81],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[82],binhdr->hunass[83],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[84],binhdr->hunass[85],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[86],binhdr->hunass[87],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[88],binhdr->hunass[89],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[90],binhdr->hunass[91],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[92],binhdr->hunass[93],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[94],binhdr->hunass[95],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[96],binhdr->hunass[97],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[98],binhdr->hunass[99],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[100],binhdr->hunass[101],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[102],binhdr->hunass[103],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[104],binhdr->hunass[105],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[106],binhdr->hunass[107],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[108],binhdr->hunass[109],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[110],binhdr->hunass[111],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[112],binhdr->hunass[113],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[114],binhdr->hunass[115],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[116],binhdr->hunass[117],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[118],binhdr->hunass[119],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[120],binhdr->hunass[121],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[122],binhdr->hunass[123],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[124],binhdr->hunass[125],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[126],binhdr->hunass[127],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[128],binhdr->hunass[129],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[130],binhdr->hunass[131],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[132],binhdr->hunass[133],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[134],binhdr->hunass[135],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[136],binhdr->hunass[137],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[138],binhdr->hunass[139],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[140],binhdr->hunass[141],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[142],binhdr->hunass[143],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[144],binhdr->hunass[145],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[146],binhdr->hunass[147],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[148],binhdr->hunass[149],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[150],binhdr->hunass[151],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[152],binhdr->hunass[153],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[154],binhdr->hunass[155],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[156],binhdr->hunass[157],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[158],binhdr->hunass[159],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[160],binhdr->hunass[161],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[162],binhdr->hunass[163],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[164],binhdr->hunass[165],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[166],binhdr->hunass[167],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    HHTOU(binhdr->hunass[168],binhdr->hunass[169],utemp)
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
 } else {
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[0],binhdr->hunass[1])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[2],binhdr->hunass[3])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[4],binhdr->hunass[5])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[6],binhdr->hunass[7])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[8],binhdr->hunass[9])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[10],binhdr->hunass[11])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[12],binhdr->hunass[13])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[14],binhdr->hunass[15])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[16],binhdr->hunass[17])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[18],binhdr->hunass[19])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[20],binhdr->hunass[21])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[22],binhdr->hunass[23])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[24],binhdr->hunass[25])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[26],binhdr->hunass[27])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[28],binhdr->hunass[29])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[30],binhdr->hunass[31])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[32],binhdr->hunass[33])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[34],binhdr->hunass[35])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[36],binhdr->hunass[37])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[38],binhdr->hunass[39])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[40],binhdr->hunass[41])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[42],binhdr->hunass[43])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[44],binhdr->hunass[45])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[46],binhdr->hunass[47])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[48],binhdr->hunass[49])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[50],binhdr->hunass[51])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[52],binhdr->hunass[53])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[54],binhdr->hunass[55])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[56],binhdr->hunass[57])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[58],binhdr->hunass[59])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[60],binhdr->hunass[61])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[62],binhdr->hunass[63])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[64],binhdr->hunass[65])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[66],binhdr->hunass[67])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[68],binhdr->hunass[69])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[70],binhdr->hunass[71])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[72],binhdr->hunass[73])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[74],binhdr->hunass[75])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[76],binhdr->hunass[77])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[78],binhdr->hunass[79])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[80],binhdr->hunass[81])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[82],binhdr->hunass[83])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[84],binhdr->hunass[85])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[86],binhdr->hunass[87])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[88],binhdr->hunass[89])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[90],binhdr->hunass[91])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[92],binhdr->hunass[93])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[94],binhdr->hunass[95])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[96],binhdr->hunass[97])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[98],binhdr->hunass[99])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[100],binhdr->hunass[101])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[102],binhdr->hunass[103])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[104],binhdr->hunass[105])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[106],binhdr->hunass[107])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[108],binhdr->hunass[109])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[110],binhdr->hunass[111])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[112],binhdr->hunass[113])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[114],binhdr->hunass[115])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[116],binhdr->hunass[117])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[118],binhdr->hunass[119])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[120],binhdr->hunass[121])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[122],binhdr->hunass[123])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[124],binhdr->hunass[125])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[126],binhdr->hunass[127])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[128],binhdr->hunass[129])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[130],binhdr->hunass[131])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[132],binhdr->hunass[133])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[134],binhdr->hunass[135])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[136],binhdr->hunass[137])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[138],binhdr->hunass[139])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[140],binhdr->hunass[141])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[142],binhdr->hunass[143])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[144],binhdr->hunass[145])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[146],binhdr->hunass[147])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[148],binhdr->hunass[149])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[150],binhdr->hunass[151])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[152],binhdr->hunass[153])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[154],binhdr->hunass[155])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[156],binhdr->hunass[157])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[158],binhdr->hunass[159])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[160],binhdr->hunass[161])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[162],binhdr->hunass[163])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[164],binhdr->hunass[165])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[166],binhdr->hunass[167])
    if(FALSE == xdr_u_int(segy_xdr,&utemp))status= FALSE;
    UTOHH(utemp,binhdr->hunass[168],binhdr->hunass[169])
    }

 return(status);
}
#else
void xdrbhdrsub(){return;}
#endif
