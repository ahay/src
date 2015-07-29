

/*
  Copyright (C) 2010 Colorado School of Mines
  
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

#include "weiutl.h"
/*^*/

#ifndef _weiutl_h

/*------------------------------------------------------------*/
typedef struct weicub *weicub3d;
/*^*/

struct weicub{
    bool    verb;
    sf_axis  amx;
    sf_axis  amy;
    sf_axis  alx;
    sf_axis  aly;
    sf_axis   az;
    sf_axis   aw;
    sf_axis   af;
    sf_axis   ae;
    sf_axis   ac;
    sf_axis  ahx;
    sf_axis  ahy;
    sf_axis  ahz;
    sf_axis  aht;
    int   ompnth; /* number of threads */
    int      nxy;
    int  pmx,pmy; /* padding */
    int  tmx,tmy; /*  taper  */
    float    eps;
    float  dsmax;
    float  dtmax;
    bool    dflg;
    float   deps;
};
/*^*/
/*------------------------------------------------------------*/
typedef struct weif2d *weifft2d;
/*^*/

struct weif2d{
    int            n1,n2;
    kiss_fft_cfg   forw1; /*   FFT on axis 1 */
    kiss_fft_cfg   invs1;
    kiss_fft_cfg   forw2; /*   FFT on axis 2 */
    kiss_fft_cfg   invs2;
    kiss_fft_cpx  *ctmp1; /* temp array */
    kiss_fft_cpx  *ctmp2; /* temp array */
    float       fftscale; /* FFT scale  */
};
/*^*/

/*------------------------------------------------------------*/
typedef struct weissr *weissr3d;
/*^*/

struct weissr{
    sf_axis     bxx,byy;
    float         **kk;
    int        *lx,*ly;
    sf_complex   ***pk;
    sf_complex   ***wk;
    float        ***wt;
    float        ****wt1;
    float       dsmax2;
    weifft2d      *f2d;
};
/*^*/

/*------------------------------------------------------------*/
typedef struct weitap *weitap3d;
/*^*/

struct weitap{
    int     nt1,   nt2;
    int      n1,    n2;
    bool     b1,    b2;
    float *tap1, *tap2;
};
/*^*/

/*------------------------------------------------------------*/
typedef struct weislo *weislo3d;
/*^*/

struct weislo{
    sf_file       F; /* file */
    int         *nr; /* number of references */
    float      **sm; /* ref slo squared */
    int       nrmax;
    float     dsmax;
    float      ***s;
};
/*^*/


/*------------------------------------------------------------*/
typedef struct weiopf *weiop3f;
/*^*/

struct weiopf{
    sf_complex ***wfld; /*  generic wavefield (nx,ny,nth) */
    sf_complex ***swfl; /*   source wavefield (nx,ny,nth) */
    sf_complex ***rwfl; /* receiver wavefield (nx,ny,nth) */
    /* CIC */
    float     ****icic; /* image (nx,ny,nz) */
    /* EIC */
    float    *****ieic; /*     image (nhx,nhy,nhz,nht,nc) */
    float     ****ihic; /*     image (nhx,nhy,    nht,nc) */
    /* tmp arrays */
    sf_complex ***sone; /*    source wavefield (nx,ny,nz) */
    sf_complex ***rone; /*  receiver wavefield (nx,ny,nz) */
};
/*^*/

/*------------------------------------------------------------*/
typedef struct weiop *weiop3d;
/*^*/

struct weiop{
    sf_complex ***wfld; /*  generic wavefield (nx,ny,nth) */
    sf_complex ***swfl; /*   source wavefield (nx,ny,nth) */
    sf_complex ***rwfl; /* receiver wavefield (nx,ny,nth) */
    sf_complex ***ahic; /*     adjoint source (nx,ny,nth) */
    /* CIC */
    float      ***icic; /* image (nx,ny,nz) */
    /* EIC */
    float    *****ieic; /*     image (nhx,nhy,nhz,nht,nc) */
    float     ****ihic; /*     image (nhx,nhy,    nht,nc) */
    /* tmp arrays */
    sf_complex ***sone; /*    source wavefield (nx,ny,nz) */
    sf_complex ***rone; /*  receiver wavefield (nx,ny,nz) */
};
/*^*/

/*------------------------------------------------------------*/
typedef struct weico *weico3d;
/*^*/

struct weico{
    pt3d          *cc; /* CIP coordinates */
    bool        *ccin; /* in-grid flag */
    float cxmin,cxmax;
    float cymin,cymax;
    float czmin,czmax;
    int      **mcxall; /* indices */
    int      **pcxall;
    int      **mcyall;
    int      **pcyall;
    int      **mczall;
    int      **pczall;
    int       *iczall;
    sf_complex   **tt; /* time delay */
    int        *ncofz;
    int       **icofz;
};
/*^*/

#endif
