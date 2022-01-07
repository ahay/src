/* Helper functions for Gauss-Seidel */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
#include <rsf.h>

/* Phase space grid */
static int nz,nx,na,iperiod;
static float oz,ox,oa;
static float dz,dx,da;
static float dzi,dxi,dai;
static int order = 1, lorder = 0;
static const char code[] = "xatzlid";

void gsray_init(int nz1, int nx1, int na1,
                float oz1, float ox1, float oa1,
                float dz1, float dx1, float da1)
/*< initialize >*/
{
    nz = nz1; oz = oz1; dz = dz1;
    nx = nx1; ox = ox1; dx = dx1;
    na = na1; oa = oa1; da = da1;
    iperiod = 2*nx + 2*nz - 4;

    dzi = -1.0/dz;
    dxi = -1.0/dx;
    dai = -1.0/da;
}

void gsray_order(int ord, int lord)
/*< set order of upwind >*/
{
    order = ord;
    lorder  = lord;
}

static float sf_gsray_shift_index (float f, int shift) {
    f -= shift;
    if (f < 0.0)
        f += iperiod;
    return f;
}

static float sf_gsray_unshift_index (float f, int shift) {
    f += shift;
    if (f > iperiod)
        f -= iperiod;
    return f;
}

#define IEPS 0.01
static float sf_gsray_choose_index (float fz, float fx, float fa,
                                    float fzc, float fxc, float fac, float dd) {
    int shift;
    float ts1, ts2, ts3, ts4;
    float fz1, fx1, fa1;

    /* Left upper corner */
    ts1 = (fz*fzc + fx*fxc + fa*fac)/dd;
    /* Right upper corner */
    shift = nx - 1;
    fz1 =  sf_gsray_shift_index (fz, shift);
    fx1 =  sf_gsray_shift_index (fx, shift);
    fa1 =  sf_gsray_shift_index (fa, shift);
    ts2 = (fz1*fzc + fx1*fxc + fa1*fac)/dd;
    ts2 = sf_gsray_unshift_index (ts2, shift);
    /* Right lower corner */
    shift += (nz - 1);
    fz1 =  sf_gsray_shift_index (fz, shift);
    fx1 =  sf_gsray_shift_index (fx, shift);
    fa1 =  sf_gsray_shift_index (fa, shift);
    ts3 = (fz1*fzc + fx1*fxc + fa1*fac)/dd;
    ts3 = sf_gsray_unshift_index (ts3, shift);
    /* Left lower corner */
    shift += (nx - 1);
    fz1 =  sf_gsray_shift_index (fz, shift);
    fx1 =  sf_gsray_shift_index (fx, shift);
    fa1 =  sf_gsray_shift_index (fa, shift);
    ts4 = (fz1*fzc + fx1*fxc + fa1*fac)/dd;
    ts4 = sf_gsray_unshift_index (ts4, shift);
    /* Choose proper solution */
    if (fabsf (ts1 - ts2) < IEPS && fabsf (ts1 - ts3) < IEPS)
        return ts1;
    else if (fabsf (ts1 - ts2) < IEPS && fabsf (ts2 - ts4) < IEPS)
        return ts2;
    else if (fabsf (ts2 - ts3) < IEPS && fabsf (ts3 - ts4) < IEPS)
        return ts3;
    else if (fabsf (ts1 - ts4) < IEPS && fabsf (ts4 - ts3) < IEPS)
        return ts4;
    else
        sf_error ("Index ambiguity");

    /*    if (fabsf (ts1 - ts2) < IEPS)
        return 0.5*(ts1 + ts2);
    else if (fabsf (ts1 - ts3) < IEPS)
        return 0.5*(ts1 + ts3);
    else if (fabsf (ts1 - ts4) < IEPS)
        return 0.5*(ts1 + ts4);
    else if (fabsf (ts2 - ts3) < IEPS)
        return 0.5*(ts2 + ts3);
    else if (fabsf (ts2 - ts4) < IEPS)
        return 0.5*(ts2 + ts4);
    else if (fabsf (ts3 - ts4) < IEPS)
        return 0.5*(ts3 + ts4);
    else
        sf_error ("Index ambiguity");*/
    return ts1;
}

float gs_update(float ***t, 
                float Hp, float Hq, float da_dsigma,
                int iz, int ix, int ia, 
                float dt_dsigma, float dl_dsigma, float dxaz_dsigma,
                int iq /* what to compute */, int iter,
                float zn, float zf, float xn, float xf, float an, float af)
/*< Gauss-Seidel update >*/
{
    float ts, fz, fx, fa, dd, omega = 1.0;

    float fzmax, fzmin;
    float fxmax, fxmin;
    float famax, famin;

    /* Downwind scheme coefficients  */
    fz = dzi*Hp;  /* cs*ss; */
    fx = dxi*Hq;  /* sn*ss; */
    /* fa = dai*(cs*ssx - sn*ssz);*/
    fa = dai*da_dsigma;

    fzmin = SF_MAX(-fz,0.);
    fzmax = SF_MAX(fz,0.);

    fxmin = SF_MAX(-fx,0.);
    fxmax = SF_MAX(fx,0.);

    famin = SF_MAX(-fa,0.);
    famax = SF_MAX(fa,0.);

    /* Diagonal term */
    dd = (fxmax + fxmin + fzmax + fzmin + famax + famin);

    ts = 0.0;

    int izm2, izm1, izp1, izp2;
    int ixm2, ixm1, ixp1, ixp2;
    int iam2, iam1, iap1, iap2;
    izp1 = iz + 1; izp2 = iz + 2;
    izm1 = iz - 1; izm2 = iz - 2;
    ixp1 = ix + 1; ixp2 = ix + 2;
    ixm1 = ix - 1; ixm2 = ix - 2;
    iap1 = ia + 1; iap2 = ia + 2;
    iam1 = ia - 1; iam2 = ia - 2;
    /* z-direction */
    if (iz == 0) {
        izm1 = 0; izm2 = 0;
    } else if (iz == 1) {
        izm2 = 0;
    } else if (iz == (nz - 1)) {
        izp1 = nz - 1; izp2 = nz - 1;
    } else if (iz == (nz - 2)) {
        izp2 = nz - 1;
    }
    /* x-direction */
    if (ix == 0) {
        ixm1 = 0; ixm2 = 0;
    } else if (ix == 1) {
        ixm2 = 0;
    } else if (ix == (nx - 1)) {
        ixp1 = nx - 1; ixp2 = nx - 1;
    } else if (ix == (nx - 2)) {
        ixp2 = nx - 1;
    }
    /* a-direction */
    if (ia == 0) {
        iam1 = na - 1; iam2 = na - 2;
    } else if (ia == 1) {
        iam2 = na - 1;
    } else if (ia == (na - 1)) {
        iap1 = 0; iap2 = 1;
    } else if (ia == (na - 2)) {
        iap2 = 0;
    }

    if (4 == order && iter >= lorder) { /* QUICK, finite-volume formulation */
        /* Z-direction, near face */
        dd = 0.5*zn*dzi;
        if (zn > 0.0) {
            dd += -0.125*zn*dzi;
            ts += zn*dzi*(0.75*t[ia][ix][izm1] - 0.125*t[ia][ix][izm2]);
        } else if (zn < 0.0) {
            dd += 0.25*zn*dzi;
            ts += zn*dzi*(0.375*t[ia][ix][izm1] - 0.125*t[ia][ix][izp1]);
        }
        /* Z-direction, far face */
        dd += -0.5*zf*dzi;
        if (zf > 0.0) {
            dd += -0.25*zf*dzi;
            ts += zf*dzi*(-0.375*t[ia][ix][izp1] + 0.125*t[ia][ix][izm1]);
        } else if (zf < 0.0) {
            dd += 0.125*zf*dzi;
            ts += zf*dzi*(-0.75*t[ia][ix][izp1] + 0.125*t[ia][ix][izp2]);
        }
        /* X-direction, near face */
        dd += 0.5*xn*dxi;
        if (xn > 0.0) {
            dd += -0.125*xn*dxi;
            ts += xn*dxi*(0.75*t[ia][ixm1][iz] - 0.125*t[ia][ixm2][iz]);
        } else if (xn < 0.0) {
            dd += 0.25*xn*dxi;
            ts += xn*dxi*(0.375*t[ia][ixm1][iz] - 0.125*t[ia][ixp1][iz]);
        }
        /* X-direction, far face */
        dd += -0.5*xf*dxi;
        if (xf > 0.0) {
            dd += -0.25*xf*dxi;
            ts += xf*dxi*(-0.375*t[ia][ixp1][iz] + 0.125*t[ia][ixm1][iz]);
        } else if (xf < 0.0) {
            dd += 0.125*xf*dxi;
            ts += xf*dzi*(-0.75*t[ia][ixp1][iz] + 0.125*t[ia][ixp2][iz]);
        }
        /* Angle-direction, near face */
        dd += 0.5*an*dai;
        if (an > 0.0) {
            dd += -0.125*an*dai;
            ts += an*dai*(0.75*t[iam1][ix][iz] - 0.125*t[iam2][ix][iz]);
        } else if (an < 0.0) {
            dd += 0.25*an*dai;
            ts += an*dai*(0.375*t[iam1][ix][iz] - 0.125*t[iap1][ix][iz]);
        }
        /* Angle-direction, far face */
        dd += -0.5*af*dai;
        if (af > 0.0) {
            dd += -0.25*af*dai;
            ts += af*dai*(-0.375*t[iap1][ix][iz] + 0.125*t[iam1][ix][iz]);
        } else if (af < 0.0) {
            dd += 0.125*af*dai;
            ts += af*dai*(-0.75*t[iap1][ix][iz] + 0.125*t[iap2][ix][iz]);
        }
        ts = -ts;
        omega = 0.35;      
//  } else if (4 == order) { /* QUICK */
/*      dd *= 3.0/8.0;
        if (fzmax != 0.0) {
            ts += 1.0/8.0*fzmax*(-t[ia][ix][izp2] + 7.0*t[ia][ix][izp1] - 3.0*t[ia][ix][izm1]);
        } else if (fzmin != 0.0) {
            ts += 1.0/8.0*fzmin*(-t[ia][ix][izm2] + 7.0*t[ia][ix][izm1] - 3.0*t[ia][ix][izp1]);
        }
        if (fxmax != 0.0) {
            ts += 1.0/8.0*fxmax*(-t[ia][ixp2][iz] + 7.0*t[ia][ixp1][iz] - 3.0*t[ia][ixm1][iz]);
        } else if (fxmin != 0.0) {
            ts += 1.0/8.0*fxmin*(-t[ia][ixm2][iz] + 7.0*t[ia][ixm1][iz] - 3.0*t[ia][ixp1][iz]);
        }
        if (famax != 0.0) {
            ts += 1.0/8.0*famax*(-t[iap2][ix][iz] + 7.0*t[iap1][ix][iz] - 3.0*t[iam1][ix][iz]);
        } else if (famin != 0.0) {
            ts += 1.0/8.0*famin*(-t[iam2][ix][iz] + 7.0*t[iam1][ix][iz] - 3.0*t[iap1][ix][iz]);
        }
        omega = 0.5;*/
    } else if (3 == order) { /* Third-order upwind */
        dd *= 0.5;
        if (fzmax != 0.0) {
            ts += 1.0/6.0*fzmax*(-t[ia][ix][izp2] + 6.0*t[ia][ix][izp1] - 2.0*t[ia][ix][izm1]);
        } else if (fzmin != 0.0) {
            ts += 1.0/6.0*fzmin*(-t[ia][ix][izm2] + 6.0*t[ia][ix][izm1] - 2.0*t[ia][ix][izp1]);
        }
        if (fxmax != 0.0) {
            ts += 1.0/6.0*fxmax*(-t[ia][ixp2][iz] + 6.0*t[ia][ixp1][iz] - 2.0*t[ia][ixm1][iz]);
        } else if (fxmin != 0.0) {
            ts += 1.0/6.0*fxmin*(-t[ia][ixm2][iz] + 6.0*t[ia][ixm1][iz] - 2.0*t[ia][ixp1][iz]);
        }
        if (famax != 0.0) {
            ts += 1.0/6.0*famax*(-t[iap2][ix][iz] + 6.0*t[iap1][ix][iz] - 2.0*t[iam1][ix][iz]);
        } else if (famin != 0.0) {
            ts += 1.0/6.0*famin*(-t[iam2][ix][iz] + 6.0*t[iam1][ix][iz] - 2.0*t[iap1][ix][iz]);
        }
        omega = 0.5;
    } else if (2 == order) { /* Second-order upwind */
        dd *= 1.5;
        ts += 0.5*(fzmax*(-t[ia][ix][izp2] + 4.0*t[ia][ix][izp1]) +
                   fzmin*(-t[ia][ix][izm2] + 4.0*t[ia][ix][izm1]) +
                   fxmax*(-t[ia][ixp2][iz] + 4.0*t[ia][ixp1][iz]) +
                   fxmin*(-t[ia][ixm2][iz] + 4.0*t[ia][ixm1][iz]) +
                   famax*(-t[iap2][ix][iz] + 4.0*t[iap1][ix][iz]) +
                   famin*(-t[iam2][ix][iz] + 4.0*t[iam1][ix][iz]));
        omega = 1.0;
    } else { /* First-order upwind */
        fz = fx = fa = 0.0;
        if (fzmax != 0.0)
            fz = t[ia][ix][izp1];
        else if (fzmin != 0.0)
            fz = t[ia][ix][izm1];
        if (fxmax != 0.0)
            fx = t[ia][ixp1][iz];
        else if (fxmin != 0.0)
            fx = t[ia][ixm1][iz];
        if (famax != 0.0)
            fa = t[iap1][ix][iz];
        else if (famin != 0.0)
            fa = t[iam1][ix][iz];
        /* Special case - exit index, check periodicity */
        if ('i' == code[iq]) {
            return sf_gsray_choose_index (fz, fx, fa, fzmin + fzmax,
                                          fxmin + fxmax, famin + famax, dd);
        }
        ts += ((fzmin + fzmax)*fz +
               (fxmin + fxmax)*fx +
               (famin + famax)*fa);
        omega = 1.0;
/*
        ts += (fzmax*t[ia][ix][izp1] + fzmin*t[ia][ix][izm1] +
               fxmax*t[ia][ixp1][iz] + fxmin*t[ia][ixm1][iz] +
               famax*t[iap1][ix][iz] + famin*t[iam1][ix][iz]);*/
    }
    switch (code[iq]) {
	case 't':
	    ts += dt_dsigma; /* ss2; ss*ss */
	    break;
	case 'l':
	    ts += dl_dsigma; /* ss in anisotropic wave p. */
	    break;
	default:
	    ts += dxaz_dsigma;
	    break;
    }
    /* G-S + SOR */
    return (1.0 - omega)*t[ia][ix][iz] + omega*ts/dd;
}

char gs_color(char ***c, 
              int iz, int ix, int ia, 
              float ss, float ssz, float ssx, /* slowness and slowness gradient */
              float cs, float sn /* cosine and sine of the phase angle */)
/*< Gauss-Seidel "color" >*/
{
    char cz;

    float fz;

    float fzmax, fzmin; 

    /* Downwind scheme coefficients  */
    fz = dzi*cs*ss;
    
    fzmin = SF_MAX(-fz,0.);
    fzmax = SF_MAX(fz,0.);
    
    /* z-direction */
    if (fzmax > 0. && iz < nz-1) {
        cz = c[ia][ix][iz+1];
    } else if (fzmin > 0. && iz > 0) {
        cz = c[ia][ix][iz-1];
    } else {
        cz = '0';
    }

    return cz;
}

void boundary_mat(float ***t, int iz, int ix, int ia, int iq)
/*< Initialize outgoing rays for boundary conditions >*/
{
    switch(code[iq]) {
        case 'x':
            t[ia][ix][iz] = ox + ix*dx;
            break;
        case 'a':
            /* convert to degrees */
            t[ia][ix][iz] = (oa + ia*da)*180./SF_PI;
            break;
        case 't':
        case 'l':
            t[ia][ix][iz] = 0;
            break;
        case 'z':
            t[ia][ix][iz] = oz + iz*dz;
            break;
        case 'd':
            t[ia][ix][iz] = 0;
            break;
    }

    return;
}

static float gs_hcircle_dist(float s, int iz, int ix, int ia, int iq) {
    double R = 0.5*(nz + 1)*dz;
    double a = oa + ia*da;
    double x0, tn, xi, yi, len;
    if (a < 0.)
        a = a + SF_PI;
    else
        a = SF_PI - a;
    tn = tan(a);
    /* Cricle-centered x,y coordinates from here */
    x0 = iz*dz - R;
    /* Distance from the center of the circle (middle of real z) */
    if (a > SF_PI/2.0)
        xi = -(sqrt((tn*tn + 1)*R*R - tn*tn*x0*x0) - tn*tn*x0)/(tn*tn + 1);
    else
        xi = (sqrt((tn*tn + 1)*R*R - tn*tn*x0*x0) + tn*tn*x0)/(tn*tn + 1);
    if (xi*xi < R*R)
        yi = sqrt(R*R - xi*xi);
    else
        yi = 0.0;
    yi *= dx/dz; /* Distance in real x */
    len = sqrt((xi - x0)*(xi - x0) + yi*yi); /* Ray length */

    switch(code[iq]) {
        case 'x': {
            if (ix == (nx - 1))
                return (ox + nx*dx + yi);
            else
                return ox - yi*dx/dz;
            }
        case 'a':
            return (oa + ia*da);
        case 't':
            return len*s;
        case 'l':
            return len;
        case 'z':
            return (oz + R + xi);
        case 'd':
            return 0.0;
    }

    return 0.0;
}

void gs_init_bc(float ***t, float **s, bool sph, int iz, int ix, int ia, int iq)
/*< Initialize outgoing rays for boundary conditions >*/
{
    /* On left or right side, calculate intersection with imaginary
       half-spheres */
    if (sph && (ix == 0 || ix == (nx - 1))) {
        t[ia][ix][iz] = gs_hcircle_dist(s[ix][iz], iz, ix, ia, iq);
        return;
    }

    switch(code[iq]) {
        case 'x':
            t[ia][ix][iz] = ox + ix*dx;
            break;
        case 'a':
            /* convert to degrees */
            t[ia][ix][iz] = (oa + ia*da)*180./SF_PI;
            break;
        case 't':
        case 'l':
            t[ia][ix][iz] = 0;
            break;
        case 'z':
            t[ia][ix][iz] = oz + iz*dz;
            break;
        case 'd':
            t[ia][ix][iz] = 0;
            break;
        case 'i': {
            if (0 == iz) /* Top */
                t[ia][ix][iz] = ix;
            else if ((nz - 1) == iz) /* Bottom */
                t[ia][ix][iz] = 2*nx + nz - 3 - ix;
            else if (0 == ix) /* Left */
                t[ia][ix][iz] = 2*nx + 2*nz - 4 - iz;
            else if ((nx - 1) == ix) /* Right */
                t[ia][ix][iz] = nx - 1 + iz;
            break;
        }
    }

    return;
}

void gs_get_face_cf(float **s, float **sz, float **sx, int iz, int ix, int ia,
                    float *zn, float *zf, float *xn, float *xf, float *an, float *af)
/*< Get coefficients for the faces of the control volume (relevant for order=4 only) >*/
{
    float aa, aap, aan, ss, ssz, ssx, cs, sn;
    int izm1, izp1, ixm1, ixp1, iam1, iap1;
    izp1 = iz + 1; izm1 = iz - 1;
    ixp1 = ix + 1; ixm1 = ix - 1;
    iap1 = ia + 1; iam1 = ia - 1;
    /* z-direction */
    if (iz == 0) {
        izm1 = 0;
    } else if (iz == (nz - 1)) {
        izp1 = nz - 1;
    }
    /* x-direction */
    if (ix == 0) {
        ixm1 = 0;
    } else if (ix == (nx - 1)) {
        ixp1 = nx - 1;
    }
    /* a-direction */
    if (ia == 0) {
        iam1 = na - 1;
    } else if (ia == (na - 1)) {
        iap1 = 0;
    }
    aa = oa + ia*da;
    aap = oa + iam1*da;
    aan = oa + iap1*da;
    cs = cos(aa); sn = sin(aa);
    ss = s[ix][iz];
    ssz = sz[ix][iz];
    ssx = sx[ix][iz];
    /* Average cell face values */
    *zn = 0.5*cs*(s[ix][izm1] + ss);
    *zf = 0.5*cs*(s[ix][izp1] + ss);
    *xn = 0.5*sn*(s[ixm1][iz] + ss);
    *xf = 0.5*sn*(s[ixp1][iz] + ss);
    *an = 0.5*ssx*(cos(aap) + cs) - 0.5*ssz*(sin(aap) + sn);
    *af = 0.5*ssx*(cos(aan) + cs) - 0.5*ssz*(sin(aan) + sn);
}

