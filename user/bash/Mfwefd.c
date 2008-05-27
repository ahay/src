/* Full wave equation 3D finite-difference modeling */
/*
  Copyright (C) 2008 The University of Texas at Austin
  
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

#define C0  1.125 /*  9/8  */	
#define C1  0.04166666666666666667 /* 1/24 */

/* centered FD derivative stencils */
#define DZA(a,i3,i2,i1,h) (C0 * (a[i3][i2][i1 + 1] - a[i3][i2][i1]) -  \
                           C1 * (a[i3][i2][i1 + 2] - a[i3][i2][i1 - 1])) * h
#define DXA(a,i3,i2,i1,h) (C0 * (a[i3][i2 + 1][i1] - a[i3][i2][i1]) -  \
                           C1 * (a[i3][i2 + 2][i1] - a[i3][i2 - 1][i1])) * h
#define DYA(a,i3,i2,i1,h) (C0 * (a[i3 + 1][i2][i1] - a[i3][i2][i1]) -  \
                           C1 * (a[i3 + 2][i2][i1] - a[i3 - 1][i2][i1])) * h
#define DZB(a,i3,i2,i1,h) (C0 * (a[i3][i2][i1] - a[i3][i2][i1 - 1]) -  \
                           C1 * (a[i3][i2][i1 + 1] - a[i3][i2][i1 - 2])) * h
#define DXB(a,i3,i2,i1,h) (C0 * (a[i3][i2][i1] - a[i3][i2 - 1][i1]) -  \
                           C1 * (a[i3][i2 + 1][i1] - a[i3][i2 - 2][i1])) * h
#define DYB(a,i3,i2,i1,h) (C0 * (a[i3][i2][i1] - a[i3 - 1][i2][i1]) -  \
                           C1 * (a[i3 + 1][i2][i1] - a[i3 - 2][i2][i1])) * h


int main (int argc, char* argv[]) {
    /* Input */
    sf_file Fwav = NULL; /* wavelet */
    sf_file Fsou = NULL; /* sources */
    sf_file Frec = NULL; /* receivers */
    sf_file Fvp = NULL; /* compressional velocities */
    sf_file Fvs = NULL; /* shear velocities */
    sf_file Fden = NULL; /* densities */
    /* Output */
    sf_file Fdat = NULL; /* recorded data */
    sf_file Fwfl = NULL; /* computed wavefield */

    /* cube axes */
    sf_axis at, az, ax, ay, as, ar;
    int nt, nx, ny, nz, ns, nr;
    int it, iz, ix, iy;
    float dt, dx, dy, dz;
    float buox, buoy, buoz;
    float muxy, muxz, muyz;

    float ***lam = NULL; /* Lame constant - Lambda */
    float ***mu = NULL; /* Lame constant - Mu */
    float ***buo = NULL; /* Buoyancy - 1 over density */

    float **srcs = NULL; /* Sources, 2D array of x,y,z triplets */
    float **recs = NULL; /* Receivers, 2D array of x,y,z triplets */

    float **wav = NULL; /* Wavelet, 2D array of triplets of the force components */

    float ***vx = NULL;
    float ***vy = NULL;
    float ***vz = NULL;
    float ***txx = NULL;
    float ***tyy = NULL;
    float ***tzz = NULL;
    float ***txy = NULL;
    float ***txz = NULL;
    float ***tyz = NULL;

    bool verb, fsrf;

    sf_init(argc,argv);

    if (!sf_getbool ("verb", &verb)) verb = false; /* verbosity flag */
    if (!sf_getbool ("free", &fsrf)) fsrf = false; /* free surface flag */

    Fwav = sf_input ("in"); /* wavelet */
//    Fsou = sf_input ("sou"); /* sources */
//    Frec = sf_input ("rec"); /* receivers */
    Fvp = sf_input ("vp"); /* compressional velocities */
    Fvs = sf_input ("vs"); /* shear velocities */
    Fden = sf_input ("den"); /* densities */

    Fwfl = sf_output("out"); /* computed wavefield */
    Fdat = sf_output("dat"); /* recorded data */

    /* axes */
    at = sf_iaxa (Fwav, 2); sf_setlabel (at, "t"); if (verb) sf_raxa (at); /* time */
//    as = sf_iaxa (Fsou, 2); sf_setlabel (as, "s"); if (verb) sf_raxa (as); /* sources */
//    ar = sf_iaxa (Frec, 2); sf_setlabel (ar, "r"); if (verb) sf_raxa (ar); /* receivers */
    az = sf_iaxa (Fvp, 1); sf_setlabel (az, "z"); if (verb) sf_raxa (az); /* depth - z */
    ax = sf_iaxa (Fvp, 2); sf_setlabel (ax, "x"); if (verb) sf_raxa (ax); /* space - x */
    ay = sf_iaxa (Fvp, 3); sf_setlabel (ax, "y"); if (verb) sf_raxa (ay); /* space - y */

    nt = sf_n (at); dt = sf_d (at);
//    ns = sf_n (as);
//    nr = sf_n (ar);
    nz = sf_n (az); dz = sf_d (az);
    nx = sf_n (ax); dx = sf_d (ax);
    ny = sf_n (ay); dy = sf_d (ay);

    /* Read wavelet */
    wav = sf_floatalloc2 (3, nt);
    sf_floatread (wav[0], 3 * nt, Fwav);

    /* Read sources */
//    srcs = sf_floatalloc2 (3, ns);
//    sf_floatread (srcs[0], 3 * ns, Fsou);

    /* Read receivers */
//    recs = sf_floatalloc2 (3, nr);
//    sf_floatread (recs[0], 3 * nr, Frec);

    /* Read velocities and densities */
    lam = sf_floatalloc3 (nz, nx, ny);
    sf_floatread (lam[0][0], nx * ny * nz, Fvp);

    mu = sf_floatalloc3 (nz, nx, ny);
    sf_floatread (mu[0][0], nx * ny * nz, Fvs);

    buo = sf_floatalloc3 (nz, nx, ny);
    sf_floatread (buo[0][0], nx * ny * nz, Fden);

    /* Particle velocity fields + stress fields */
    vx = sf_floatalloc3 (nz, nx, ny);
    vy = sf_floatalloc3 (nz, nx, ny);
    vz = sf_floatalloc3 (nz, nx, ny);
    txx = sf_floatalloc3 (nz, nx, ny);
    tyy = sf_floatalloc3 (nz, nx, ny);
    tzz = sf_floatalloc3 (nz, nx, ny);
    txy = sf_floatalloc3 (nz, nx, ny);
    txz = sf_floatalloc3 (nz, nx, ny);
    tyz = sf_floatalloc3 (nz, nx, ny);

    sf_oaxa (Fwfl, az, 1);
    sf_oaxa (Fwfl, ax, 2);
    sf_oaxa (Fwfl, ay, 3);

    /* Recalculate them into Lame constants and buoyancy + initialize wavefield arrays */
    for (iy = 0; iy < ny; iy++) {
        for (ix = 0; ix < ny; ix++) {
            for (iz = 0; iz < nz; iz++) {
                mu[iy][ix][iz] = mu[iy][ix][iz] * mu[iy][ix][iz] * buo[iy][ix][iz];
                lam[iy][ix][iz] = lam[iy][ix][iz] * lam[iy][ix][iz] * buo[iy][ix][iz] - 2.0 * mu[iy][ix][iz];
                buo[iy][ix][iz] = 1.0 / buo[iy][ix][iz];
                vx[iy][ix][iz] = 0.0;
                vy[iy][ix][iz] = 0.0;
                vz[iy][ix][iz] = 0.0;
                txx[iy][ix][iz] = 0.0;
                txy[iy][ix][iz] = 0.0;
                txz[iy][ix][iz] = 0.0;
                txy[iy][ix][iz] = 0.0;
                txz[iy][ix][iz] = 0.0;
                tyz[iy][ix][iz] = 0.0;
            }
        }
    }
    /* Main loop over time */
    for (it = 0; it < nt; it++) {
        if (verb)
            sf_warning ("Calculating time stamp: %f", it * dt);
        for (iy = 2; iy < (ny - 2); iy++) {
            for (ix = 2; ix < (nx - 2); ix++) {
                for (iz = 2; iz < (nz - 2); iz++) {
                    buox = 0.5 * (buo[iy][ix][iz] + buo[iy][ix - 1][iz]);
                    buoy = 0.5 * (buo[iy][ix][iz] + buo[iy - 1][ix][iz]);
                    buoz = 0.5 * (buo[iy][ix][iz] + buo[iy][ix][iz - 1]);
                    vx[iy][ix][iz] = vx[iy][ix][iz] + dt * buox *
                                     (DXB (txx, iy, ix, iz, 1.0 / dx) + 
                                      DYB (txy, iy, ix, iz, 1.0 / dy) + 
                                      DZB (txz, iy, ix, iz, 1.0 / dz));
                    vy[iy][ix][iz] = vy[iy][ix][iz] + dt * buoy *
                                     (DXB (txy, iy, ix, iz, 1.0 / dx) + 
                                      DYB (tyy, iy, ix, iz, 1.0 / dy) + 
                                      DZB (tyz, iy, ix, iz, 1.0 / dz));
                    vz[iy][ix][iz] = vz[iy][ix][iz] + dt * buoz *
                                     (DXB (txz, iy, ix, iz + 1, 1.0 / dx) + 
                                      DYB (tyz, iy, ix, iz + 1, 1.0 / dy) + 
                                      DZB (tzz, iy, ix, iz + 1, 1.0 / dz));
                }
            }
        }
        /* Temporary - inject one source */
        buoy = 0.5 * (buo[ny / 2][nx / 2][100] + buo[ny / 2 - 1][nx / 2][100]);
        buox = 0.5 * (buo[ny / 2][nx / 2][100] + buo[ny / 2][nx / 2 - 1][100]);
        buoz = 0.5 * (buo[ny / 2][nx / 2][100] + buo[ny / 2][nx / 2][99]);
        vx[ny / 2][nx / 2][100] += dt * buox * wav[it][0];
        vy[ny / 2][nx / 2][100] += dt * buoy * wav[it][1];
        vz[ny / 2][nx / 2][100] += dt * buoz * wav[it][2];
sf_warning ("Source: %f %f %f", vy[ny / 2][nx / 2][100], vx[ny / 2][nx / 2][100], vz[ny / 2][nx / 2][100]);
        for (iy = 2; iy < (ny - 2); iy++) {
            for (ix = 2; ix < (nx - 2); ix++) {
                for (iz = 2; iz < (nz - 2); iz++) {
                    muxy = 1.0 / (0.25 * (1.0 / mu[iy][ix][iz] + 1.0 / mu[iy][ix + 1][iz]) + 
                                          1.0 / mu[iy + 1][ix][iz] + 1.0 / mu[iy + 1][ix + 1][iz]);
                    muxz = 1.0 / (0.25 * (1.0 / mu[iy][ix][iz] + 1.0 / mu[iy][ix + 1][iz]) + 
                                          1.0 / mu[iy][ix][iz + 1] + 1.0 / mu[iy][ix + 1][iz + 1]);
                    muyz = 1.0 / (0.25 * (1.0 / mu[iy][ix][iz] + 1.0 / mu[iy + 1][ix][iz]) + 
                                          1.0 / mu[iy][ix][iz + 1] + 1.0 / mu[iy + 1][ix][iz + 1]);
                    txx[iy][ix][iz] = txx[iy][ix][iz] + dt * ((lam[iy][ix][iz] + 2 * mu[iy][ix][iz]) *
                                      DXB (vx, iy + 1, ix + 1, iz + 1, 1.0 / dx) + lam[iy][ix][iz] * 
                                      (DYB (vy, iy + 1, ix + 1, iz + 1, 1.0 / dy) + DZB (vz, iy + 1, ix + 1, iz + 1, 1.0 / dz)));
                    tyy[iy][ix][iz] = tyy[iy][ix][iz] + dt * ((lam[iy][ix][iz] + 2 * mu[iy][ix][iz]) *
                                      DYB (vy, iy + 1, ix + 1, iz + 1, 1.0 / dy) + lam[iy][ix][iz] * 
                                      (DXB (vx, iy, ix, iz, 1.0 / dx) + DZB (vz, iy, ix, iz, 1.0 / dz)));
                    tzz[iy][ix][iz] = tzz[iy][ix][iz] + dt * ((lam[iy][ix][iz] + 2 * mu[iy][ix][iz]) *
                                      DZB (vz, iy + 1, ix + 1, iz + 1, 1.0 / dz) + lam[iy][ix][iz] * 
                                      (DXB (vx, iy + 1, ix + 1, iz + 1, 1.0 / dx) + DYB (vy, iy + 1, ix + 1, iz + 1, 1.0 / dy)));
                    txy[iy][ix][iz] = txy[iy][ix][iz] + dt * muxy *
                                      (DYB (vx, iy + 1, ix + 1, iz + 1, 1.0 / dy) + DXB (vy, iy + 1, ix + 1, iz + 1, 1.0 / dx));
                    txz[iy][ix][iz] = txz[iy][ix][iz] + dt * muxz *
                                      (DZB (vx, iy + 1, ix + 1, iz + 1, 1.0 / dz) + DXB (vz, iy + 1, ix + 1, iz + 1, 1.0 / dx));
                    tyz[iy][ix][iz] = tyz[iy][ix][iz] + dt * muyz *
                                      (DZB (vy, iy + 1, ix + 1, iz + 1, 1.0 / dz) + DYB (vz, iy + 1, ix + 1, iz + 1, 1.0 / dy));
                }
            }
        }
        if (50 == it) {
            for (iy = 0; iy < ny ; iy++) {
                for (ix = 0; ix < nx; ix++) {
                    for (iz = 0; iz < nz; iz++) {
                        vx[iy][ix][iz] = 1.0 / 3.0 * (txx[iy][ix][iz] + tyy[iy][ix][iz] + tzz[iy][ix][iz]);
                    }
                }
            }
            sf_floatwrite (vx[0][0], nx * ny * nz, Fwfl);
            exit (0);
        }
    }
}
