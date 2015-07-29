/* HWT subroutines */
/*
  Copyright (C) 2006 Colorado School of Mines
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.hh>
/*^*/

#include "hwt2d.hh"

static sf_axa az,ax;
static sf_axa at;
static float **vv;
/*------------------------------------------------------------*/

void hwt2d_init(float** vv_in    /* velocity */,
		sf_axis     az_in    /* z axis   */,
		sf_axis     ax_in    /* x axis   */,
		sf_axis     at_in    /* t axis   */,
		sf_axis     ag_in    /* g axis   */)
/*< initialize hwt2d >*/
{
    az = sf_nod(az_in);
    ax = sf_nod(ax_in);
    at = sf_nod(at_in);
//    ag = sf_nod(ag_in);

    vv = vv_in;
}

/*------------------------------------------------------------*/
double hwt2d_getv(pt2d P) 
/*< get velocity from 2-D cube >*/
{
    float  z, x;
    float rz,rx;
    int   jz,jx;
    float fz,fx;
    float v;

    z = SF_MIN(SF_MAX(az.o,P.z),az.o+(az.n-2)*az.d);
    x = SF_MIN(SF_MAX(ax.o,P.x),ax.o+(ax.n-2)*ax.d);

    rz = (z-az.o)/az.d;
    jz = (int)(rz);
    fz =       rz-jz;

    rx = (x-ax.o)/ax.d;
    jx = (int)(rx);
    fx =       rx-jx;

    v = vv[ jx  ][ jz  ] * (1-fz)*(1-fx) +
	vv[ jx  ][ jz+1] * (  fz)*(1-fx) +
	vv[ jx+1][ jz  ] * (1-fz)*(  fx) +
	vv[ jx+1][ jz+1] * (  fz)*(  fx);

    return(v);
}

/*------------------------------------------------------------*/

bool hwt2d_cusp(pt2d Qo, pt2d Pm, pt2d Po, pt2d Pp)
/*< find cusp points >*/
{
  int  sjm,sjp;

  sjm = SF_SIG( JAC2d(Po,Qo,Pm) );
  sjp = SF_SIG( JAC2d(Po,Qo,Pp) );

  if (sjm*sjp==1)
    return true;
  else
    return false;
}

/*------------------------------------------------------------*/

pt2d hwt2d_raytr(pt2d Qo, pt2d Po)
/*< ray tracing >*/
{
    pt2d Pm,Pp, Ro;
    float lo;
    double sina,cosa,tana;
    int    ss;

    lo = Po.v * at.d;

    if(Po.z!=Qo.z) {
	tana = (Po.x-Qo.x) / (Po.z-Qo.z);	
	ss   = SF_SIG(tana);

	tana = SF_ABS(tana);	
	sina = tana/sqrt(1+tana*tana);
	cosa =    1/sqrt(1+tana*tana);
    } else {
	ss = 1;

	cosa = 0;
	sina = 1;
    }

    Pm.x=Po.x-lo*cosa; Pm.z=Po.z+ss*lo*sina; Pm.v=hwt2d_getv(Pm);
    Pp.x=Po.x+lo*cosa; Pp.z=Po.z-ss*lo*sina; Pp.v=hwt2d_getv(Pp);

    Ro=hwt2d_step(Qo,Pm,Po,Pp);
    return(Ro);
}

/*------------------------------------------------------------*/

pt2d hwt2d_wfttr(pt2d Qo, pt2d Pm, pt2d Po, pt2d Pp)
/*< wavefront tracing >*/
{
    pt2d Ro;
    Ro=hwt2d_step(Qo,Pm,Po,Pp);
    return(Ro);
}

/*------------------------------------------------------------*/

pt2d hwt2d_step(pt2d Qo, pt2d Pm, pt2d Po, pt2d Pp)
/*< one HWT time step >*/
{
    pt2d   Ro,Sm,Sp;
    double a,b,c;
    float  lo;
    double ka, kb;

    a= Pm.x-Pp.x;
    b= Pm.z-Pp.z;
    c=(Pm.v-Pp.v)*at.d;

    lo = Po.v * at.d; /* circle radius at Po */

    ka =                            c /(a*a+b*b);
    kb = SF_SIG(a) * sqrt(a*a+b*b-c*c)/(a*a+b*b); 

    Sm.x = Po.x - lo*( ka*a - kb*b );
    Sm.z = Po.z - lo*( ka*b + kb*a );
    Sm.v = 0;

    Sp.x = Po.x - lo*( ka*a + kb*b );
    Sp.z = Po.z - lo*( ka*b - kb*a );
    Sp.v = 0;

    Ro = DST2d(Sm,Qo) > DST2d(Sp,Qo) ? Sm : Sp;
    Ro.v=hwt2d_getv(Ro);
    return(Ro);
}

/*------------------------------------------------------------*/

pt2d hwt2d_orth( pt2d Pm, pt2d Po, pt2d Pp)
/*< orthogonal step forward  >*/
{
    pt2d   Ro;
    float  lo;
    double sina,cosa,tana;
    int    ss;

/*    lo = 0.5 * DST2d(Pm,Pp);*/
    lo = Po.v * at.d;

    if(Pp.x!=Pm.x) {
	tana = (Pp.z-Pm.z) / (Pp.x-Pm.x);
	ss   = SF_SIG(tana);
	
	tana = SF_ABS(tana);	
	sina = tana/sqrt(1+tana*tana);
	cosa =    1/sqrt(1+tana*tana);
    } else {
	ss = 1;

	cosa = 0;
	sina = 1;
    }

    Ro.x = Po.x - ss * lo * sina;
    Ro.z = Po.z +      lo * cosa;
    Ro.v = 0;
    Ro.v=hwt2d_getv(Ro);
    
    return(Ro);
}
