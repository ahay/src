/* Velocity analysis.

Inverse of sfvelmod
*/
/*
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

#include <math.h>
#include <rsf.h>
#include "fint1.h"
#include <assert.h>

//static float v, x, h, v1, s;

static float Ax, Bx2, Ch2;/* half offset!*/

static float crsrk(float t, int it) 
{
    return sqrtf( (t + Ax)*(t + Ax) + t*Bx2 + t*Ch2);//+v*v
}
static float crsrk_dbg(float t, int it) 
{
    return sqrtf(t*t + Ch2);//+v*v
}
/*
static float hyperb(float t, int it) 
{
    return hypotf(t,v);
}
*/
 /******* MF
%
% Knip = a(1); Kn = a(2); B = a(3);
% m [Nm x 1]
% x [Nx x 1]
function [dT, t0, num_invalid] = f_T_mf(a,m,x)

global dT_exact V_0 V_0_grad is_V_0_grad 

Nm = size (m,1);
Nx = size (x,1);
dT = zeros (Nm,Nx);

num_invalid = 0;
invalid_value = 10;

if is_V_0_grad    V_midpoint_array = V_0_grad;else    V_midpoint_array = V_0 * ones (Nm,1); end
% V slope  s.t. m
if (size(a,1) > 3) 
    V_midpoint_array = V_0 + a(4) * m; 
end;
% V0
if (size(a,1) > 4) 
    V_midpoint_array = a(5) + a(4) * m; 
end

for im=1:Nm

    Vm = V_midpoint_array(im);
    
    for ix=1:Nx
         
         [dT(im,ix), invalid] = f_dT_mf_R(a, m(im), x(ix), Vm);
         
         num_invalid = num_invalid + invalid;
         if (invalid)
             dT(im,ix) = 15*dT_exact(im,ix); %invalid_value; BUGBUG: temporary solution assume no-error at invalid place
         end
        %dT(im,ix)=f_dT_mf(a, m(im), x(ix));
    end
end
t0 = 2 / V_midpoint_array(floor(Nm/2+0.5)) * (1/a(1)) / cos (a(3));%
num_invalid
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MF validity: noV_0_grad = 3500 + m;t evry S-R pair is valid
%
function [is_valid] = f_mf_is_valid (a,m,h)
Rnip = a(1); Rn = a(2); B = a(3);
sb = sin(B);
is_valid = 1;
xg = m + h;
xs = m - h;
% the validity condition:
if (abs (B) > 1e-6)
    if (B < -1e-6)
        x_left = Rnip / sb;
        if (xs < x_left || xg < x_left)
            is_valid = 0;
        end
    else
        x_right = Rnip / (sb);
        if (xs > x_right || xg > x_right)
            is_valid = 0;
        end
    end
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T _mfV_0_grad = 3500 + m;
% a [3 x 1] = { Knip, Kn, betha }
% 
function [dT] = f_dT_mf(a,m,h)% without V0 and t0
global V0
if (~f_mf_is_valid (a,m,h))
    dT = 0;
    return
end
%Knip = a(1); Kn = a(2); B = a(3);
Knip = 1/a(1); Kn = 1/a(2); B = a(3);
sb = sin(B);

xg = m + h;
xs = m - h;

mone = xg - xs;
mahane = (xg+xs-Knip * sb * xg * xs);

if (abs (mahane) > 1e-4)
    g = mone / mahane; %% (xg-xs) / (xg+xs-Knip * sb * xg * xs);
else
    if (abs (mone - mahane) < 1e-6)
        g = 1;
    else
        g = 1e10;
    end;
end;

% g = 0
Ks = -Kn;    
Kg = Kn;

% otherwise
if (abs (g) > 1e-6)
    G = 1 / g;
    % Source
    if (abs (G-1) < 1e-6)
        Ks = 0; % TODO: Rs = infty
    else
        Ks = (Knip - G*Kn) / (1 - G);
    end
    % Geophone
    if (abs (G+1) < 1e-6)
        Kg = 0; % TODO: Rg = infty
    else
        Kg = (Knip + G*Kn) / (1 + G);
    end
end

Ts = (sqrt (1 + 2*Ks * sb * xs + (Ks*xs)*(Ks*xs)) - 1) / Ks;

Tg = (sqrt (1 + 2*Kg * sb * xg + (Kg*xg)*(Kg*xg)) - 1) / Ks;

dT = (Ts + Tg)/V0;

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% d t2 \ da1 [da2,da3] = m [mm, hh]
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MF in Rnre Rcee
function [dT, invalid] = f_dT_mf_R(a,m,h, V_midpoint)% without V0 and t0

if (~f_mf_is_valid (a,m,h))
    invalid = 1;
    dT = 0;
    return
end
invalid = 0;

Rnip = a(1); 
Rn = a(2); 
betha = a(3);
%%%%%%%%%%%%%%%%%%%%
xg = m + h;
xs = m - h;

radCre = Rnip; 
radCee = Rn; 
%%%%%%%%%
delx = xs; 
dely = xg;
%%%%%%%%%
s_angle = 2 * sin (betha) * delx * dely;
per_trace_idlet_y = 2 * dely * sin (betha);
per_trace_idlet_x = 2 * delx * sin (betha);
%%%%%%%%%
s = delx + dely - s_angle/radCre;
if (abs (s) > 1e-10)
    sigma = (delx - dely) / s;
else
    sigma = 1e10;
end
%%%%%%%%%
if ( abs (sigma*radCee + radCre) > 1e-10)
    rplus = radCee / (sigma*radCee+radCre)*(1+sigma)*radCre;
else
    rplus = 1e10;    
end
if ( abs (-sigma*radCee + radCre) > 1e-10)
    rminus = -radCee / (-sigma*radCee+radCre)*(1-sigma)*radCre;
else
    rminus = 1e10;
end
%%%%%%%%%
delty = sqrt (rminus*rminus + rminus*per_trace_idlet_y + dely*dely) - abs (rminus);
if (rminus < 0)
    delty = -delty;
end

deltx = sqrt (rplus*rplus - rplus*per_trace_idlet_x + delx*delx) - abs (rplus);
if (rplus < 0)
    deltx = -deltx;
end

dT = (deltx - delty) / V_midpoint;

return

***********/
int main(int argc, char* argv[])
{
    fint1 crs/*nmo*/;
    bool sembl, half, dsembl, asembl, weight, dbgt0;
    int it,ih,ix,nt,nh,nx,ib,ie,nb,i, nw, mute, *mask=NULL;
    float amp, amp2, dt, dh, t0, h0, smax, num, den, dy, str, sh=0., sh2=0.;
    float *trace=NULL, *** *stack=NULL, *** *stack2=NULL, *** *stackh=NULL, *hh=NULL;
    /* char *time=NULL, *space=NULL, *unit=NULL; */
    //size_t len;
    sf_file cmp=NULL, scan=NULL, offset=NULL, msk=NULL;
    mapfunc crsfunc/*nmofunc*/;

    float a10, a20, a30, da1, da2, da3, iNth, x0, dx, x_central, V0, iV02;
    int   na1, na2, na3, ia1, ia2, ia3, Nta1a2a3;

    sf_init (argc,argv);
    cmp = sf_input("in");
    scan = sf_output("out");

    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(cmp,"n3",&nx)) sf_error("No n3= in input");
    // nx = sf_leftsize(cmp,2);

    if (!sf_getbool("dbgt0",&dbgt0)) dbgt0=false;
    /* if y compute with Rnip=Velocity */

    if (!sf_getbool("semblance",&sembl)) sembl=false;
    /* if y, compute semblance; if n, stack */
    if (sembl || !sf_getbool("diffsemblance",&dsembl)) dsembl=false;
    /* if y, compute differential semblance */
    if (sembl || dsembl || !sf_getbool("avosemblance",&asembl)) asembl=false;
    /* if y, compute AVO-friendly semblance */
    if (!sf_getint("nb",&nb)) nb=2;
    /* semblance averaging */
    if (!sf_getbool("weight",&weight)) weight=true;
    /* if y, apply pseudo-unitary weighting */

    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");


    if (!sf_histfloat(cmp,"o3",&x0)) sf_error("No o3= in input");
    if (!sf_histfloat(cmp,"d3",&dx)) sf_error("No d3= in input");

    x_central = x0 + (nx-1)/2.f*dx;


    if (!sf_getbool("half",&half)) half=false;
    /* if y, the second axis is half-offset instead of full offset */

    //CDPtype=1;
    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	hh = sf_floatalloc(nh);

	h0 = dh = 0.;
    } else {
	if (!sf_histfloat(cmp,"o2",&h0)) sf_error("No o2= in input");
	if (!sf_histfloat(cmp,"d2",&dh)) sf_error("No d2= in input");
	
	sf_putfloat(scan,"h0",h0);
	sf_putfloat(scan,"dh",dh);
	sf_putint(scan,"nh",nh);

	if (!sf_histfloat(cmp,"d3",&dy)) sf_error("No d3 in input");
	/* if (sf_histfloat(cmp,"d3",&dy)) {
	    CDPtype=half? 0.5+dh/dy: 0.5+0.5*dh/dy;
	    if (0 == CDPtype) CDPtype=1;
	    if (1 != CDPtype) {
		sf_histint(cmp,"CDPtype",&CDPtype);
		sf_warning("CDPtype=%d",CDPtype);
	    }
	}*/

	offset = NULL;
	hh = NULL;
    }

    if (NULL != sf_getstring("mask")) {
	msk = sf_input("mask");
	mask = sf_intalloc(nh);
    } else {
	msk = NULL;
	mask = NULL;
    }
    
    /* a1, a2, a3 CRS params: a3=sinBetha, a2=Rn, a1=Rnip */
    if (!sf_getfloat("a10",&a10) && !sf_histfloat(cmp,"a10",&a10)) 
	sf_error("Need a10=");
    if (!sf_getfloat("da1",&da1) && !sf_histfloat(cmp,"da1",&da1)) 
	sf_error("Need da1=");
    if (!sf_getint("na1",&na1) && !sf_histint(cmp,"na1",&na1)) 
	sf_error("Need na1=");

    if (!sf_getfloat("a20",&a20) && !sf_histfloat(cmp,"a20",&a20)) 
	sf_error("Need a20=");
    if (!sf_getfloat("da2",&da2) && !sf_histfloat(cmp,"da2",&da2)) 
	sf_error("Need da2=");
    if (!sf_getint("na2",&na2) && !sf_histint(cmp,"na2",&na2)) 
	sf_error("Need na2=");

    if (!sf_getfloat("a30",&a30) && !sf_histfloat(cmp,"a30",&a30)) 
	sf_error("Need a30=");
    if (!sf_getfloat("da3",&da3) && !sf_histfloat(cmp,"da3",&da3)) 
	sf_error("Need da3=");
    if (!sf_getint("na3",&na3) && !sf_histint(cmp,"na3",&na3)) 
	sf_error("Need na3=");


    /* SCAN */
    sf_putfloat(scan,"o2",a10);
    sf_putfloat(scan,"d2",da1);
    sf_putint(scan,"n2",na1);

    sf_putfloat(scan,"o3",a20);
    sf_putfloat(scan,"d3",da2);
    sf_putint(scan,"n3",na2);

    sf_putfloat(scan,"o4",a30);
    sf_putfloat(scan,"d4",da3);
    sf_putint(scan,"n4",na3);


    if (!sf_getfloat("smax",&smax)) smax=2.0;
    /* maximum heterogeneity */
    /*if (!sf_getint("ns",&ns)) ns=1;
     number of heterogeneity scans  
    ds = ns>1? (smax-1.0)/(ns-1): 0.;

    
    if (ns > 1) {
	sf_putfloat(scan,"o3",1.0);
	sf_putfloat(scan,"d3",ds);
	sf_putint(scan,"n3",ns);

	sf_shiftdim(cmp, scan, 3);
	} */


    /*if (!sf_getbool("slowness",&slow)) slow=false;
     if y, use slowness instead of velocity */
    
    sf_putstring(scan,"label2","a1-Rnip");
    sf_putstring(scan,"label3","a2-sin-betha");
    sf_putstring(scan,"label4","a3-Rn");

    sf_putstring(scan,"unit1","s");
    sf_putstring(scan,"unit2","km");
    sf_putstring(scan,"unit3","km");
    sf_putstring(scan,"unit4","km");

    if (dbgt0)
	crsfunc = crsrk_dbg; /*  hyperb; */
    else
	crsfunc = crsrk; /*  hyperb; */

    /* stack =  sf_floatalloc3(nt,nv,ns); 
    stack2 = (sembl || dsembl || asembl)? sf_floatalloc3(nt,nv,ns) : NULL;
    stackh = asembl? sf_floatalloc3(nt,nv,ns) : NULL; */
    stack = sf_floatalloc4(nt, na1, na2, na3);
    stack2 = (sembl || dsembl || asembl)? sf_floatalloc4(nt, na1, na2, na3) : NULL;
    stackh = asembl? sf_floatalloc4(nt, na1, na2, na3) : NULL; 

    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    if (!sf_getint("mute",&mute)) mute=12;
    /* mute zone */

    if (!sf_getfloat("str",&str)) str=0.5;
    /* maximum stretch allowed */

    if (!sf_getfloat("V0",&V0)) sf_error("V0=? float");
    /* V0 1.5..3 */
    iV02 = 2.f/V0;
    assert(V0 > 0.f);

    if (!sf_getbool("weight",&weight)) weight=true;
    /* if y, apply pseudo-unitary weighting */

    trace = sf_floatalloc(nt);
    crs/*nmo*/ = fint1_init(nw,nt,mute);

    Nta1a2a3 = nt*na1*na2*na3;
    iNth = 1.f/(nt*nh);

	for (it=0; it < Nta1a2a3/*nt*nv*ns*/; it++) {
	    stack[0][0][0][it] = 0.;
	    if (sembl || asembl) stack2[0][0][0][it] = 0.;
	    if (asembl) stackh[0][0][0][it] = 0.;
	}

    for (ix=0; ix < nx; ix++) {

	const float 	x = x0 + ix * dx - x_central,
			x2 = x*x;

	sf_warning("cmp %d of %d x=%f",ix+1,nx,x);
	
	if (NULL != offset) 
	    sf_floatread(hh,nh,offset);
	if (NULL != msk) 
	    sf_intread(mask,nh,msk);

	if (asembl) sh = sh2 = 0.;

	for (ih=0; ih < nh; ih++) {

	    const float h = (NULL != offset)? hh[ih]: 
		0.5 * (h0 + ih * dh); // + (dh/CDPtype)*(ix%CDPtype);

	    assert(!half); //if (!half) h *= 0.5; // h is half offset !

	    const float h2 = h*h;

	    sf_floatread(trace,nt,cmp); 

	    if (NULL != msk && 0==mask[ih]) continue;

	    if (asembl) {
		sh  += h;
		sh2 += h*h;
	    }

	    for (it=0; it < nt; it++) {
		trace[it] *= iNth; // /= nt*nh;
	    }
	    
	    fint1_set(crs/*nmo*/,trace);

	/* Betha */
		for (ia2=0; ia2 < na2; ia2++) {
		    const float     sb = a20 + ia2 * da2,
				    cb2 = 1.f - sb*sb,
				    K = iV02 * cb2;

			Ax = x*sb * iV02;//2 sb / V0

		    assert(fabs(sb) < 1.f);

		/* Rn */
		for (ia3=0; ia3 < na3; ia3++) {
		    const float Rn = a30 + ia3 * da3;

			Bx2 = x2*K / Rn;

		for (ia1=0; ia1 < na1; ia1++) {
		    const float Rnip = a10 + ia1 * da1;

			Ch2 = h2*K / Rnip;

		    stretch(crs/*nmo*/,crsfunc/*nmofunc*/,nt,dt,t0,nt,dt,t0,trace,str);

		    for (it=0; it < nt; it++) {

			amp = weight? fabsf(K*h/Rnip)*trace[it]: trace[it];

			if (dsembl) {

			    if (ih > 0) {
				// amp2 = amp - stack2[is][iv][it];
				amp2 = amp - stack2[ia3][ia2][ia1][it];
				//stack[is][iv][it] += amp2*amp2;
				stack[ia3][ia2][ia1][it] += amp2*amp2;
			    }
			    stack[ia3][ia2][ia1][it] = amp;

			} else {
			    if (sembl || asembl) stack2[ia3][ia2][ia1][it]  += amp*amp;
			    if (asembl) stackh[ia3][ia2][ia1][it] += amp*h;

			    stack[ia3][ia2][ia1][it] += amp;
			}
		    }

		} /* a3 */
		} /* a2 */
		} /* a1 */

		/* } s */
	} /* h */
    } /* x */
    
	if (sembl || asembl) {
	    /*for (is=0; is < ns; is++) {
	      for (iv=0; iv < nv; iv++) {*/
	    for (ia3=0; ia3 < na3; ia3++) {
	    for (ia2=0; ia2 < na2; ia2++) {
	    for (ia1=0; ia1 < na1; ia1++) {

		    for (it=0; it < nt; it++) {
			ib = it-nb;
			ie = it+nb+1;
			if (ib < 0) ib=0;
			if (ie > nt) ie=nt;
			num = 0.;
			den = 0.;
			for (i=ib; i < ie; i++) {
			    num += stack[ia3][ia2][ia1][i]*stack[ia3][ia2][ia1][i];
			    den += stack2[ia3][ia2][ia1][i];

			    /* (h2 s^2 - 2 h s sh + N sh^2)/((-h^2 + h2 N) s2) */
			    if (asembl) 
				num += (nh*stackh[ia3][ia2][ia1][i]*stackh[ia3][ia2][ia1][i] - 
					2.*sh*stack[ia3][ia2][ia1][i]*stackh[ia3][ia2][ia1][i])/sh2;
			}
			if (asembl) {
			    den *= (nh-sh*sh/sh2);
			} else {
			    den *= nh;
			}
			trace[it] = (den > 0.)? num/den: 0.;
		    }
		    sf_floatwrite(trace,nt,scan);

		} /* a1 */
		} /* a2 */
	    } /* a3 */
	} else {
	    sf_floatwrite (stack[0][0][0],Nta1a2a3/*nt*nv*ns*/,scan);
	}

    exit(0);
}
