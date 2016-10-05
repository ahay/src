/* generating acquisition geometry file for sfmpicfftrtm  */
/*
  Copyright (C) 2016 University of Texas at Austin
  
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

void geogen(int **geo, int nz, int nx, int ny, int sou_z, int sou_ox, int sou_oy, int sou_jx, int sou_jy, int sou_nx, int sou_ny, int rec_z, int rec_nx, int rec_ny, int npad, int noff, int roll);

int main(int argc, char* argv[])
{
    int nz,nx,ny,sou_z,sou_ox,sou_oy,sou_jx,sou_jy,sou_nx,sou_ny,rec_z,rec_nx,rec_ny,npad,noff,roll;
    int dim1, dim2;
    sf_axis ad1=NULL, ad2=NULL;
    sf_file Fgeo=NULL;
    int **geo=NULL;

    sf_init(argc,argv);

    if (!sf_getint("nz",&nz)) sf_error("Need nz="); /* dimension in z */
    if (!sf_getint("nx",&nx)) sf_error("Need nx="); /* dimension in x */
    if (!sf_getint("ny",&ny)) sf_error("Need ny="); /* dimension in y */
    if (!sf_getint("sou_z", &sou_z )) sf_error("Need sou_z=" ); /* source position in depth      */
    if (!sf_getint("sou_ox",&sou_ox)) sf_error("Need sou_ox="); /* source starting location in x */
    if (!sf_getint("sou_oy",&sou_oy)) sf_error("Need sou_oy="); /* source starting location in y */
    if (!sf_getint("sou_nx",&sou_nx)) sf_error("Need sou_nx="); /* number of sources in x        */
    if (!sf_getint("sou_ny",&sou_ny)) sf_error("Need sou_ny="); /* number of sources in y        */
    if (!sf_getint("sou_jx",&sou_jx)) sou_jx = (nx-sou_ox)/(sou_nx-1); /* source interval in x          */
    if (!sf_getint("sou_jy",&sou_jy)) sou_jy = (ny-sou_oy)/(sou_ny-1); /* source interval in y          */
    if (!sf_getint("rec_z", &rec_z )) sf_error("Need rec_z=" ); /* receiver position in depth */
    if (!sf_getint("rec_nx",&rec_nx)) sf_error("Need rec_nx="); /* number of receivers in x   */
    if (!sf_getint("rec_ny",&rec_ny)) sf_error("Need rec_ny="); /* number of receivers in y   */
    if (!sf_getint("npad",&npad)) sf_error("Need npad="); /* computational domain padding */
    if (!sf_getint("noff",&noff)) sf_error("Need noff="); /* near offset */
    if (!sf_getint("roll",&roll)) sf_error("Need roll="); /* acquisition pattern: 0-> fixed-spread, 1-> towed-streamer to the negative */

    /* double check dimension */
    if (sou_nx > (nx-sou_ox)/sou_jx+1) {
        sou_nx = (nx-sou_ox)/sou_jx+1;
        sf_warning("Setting sou_nx to %d",sou_nx);
    }
    if (sou_ny > 1 && sou_ny > (ny-sou_oy)/sou_jy+1) {
        sou_ny = (ny-sou_oy)/sou_jy+1;
        sf_warning("Setting sou_ny to %d",sou_ny);
    }

    /* do the work */
    dim1 = 14;
    dim2 = sou_nx*sou_ny;

    ad1 = sf_maxa(dim1,0,1); sf_setlabel(ad1,"acqpar"); sf_raxa(ad1);
    ad2 = sf_maxa(dim2,0,1); sf_setlabel(ad2,"shot");   sf_raxa(ad2);

    Fgeo = sf_output("out");
    sf_settype(Fgeo,SF_INT);
    sf_oaxa(Fgeo,ad1,1);
    sf_oaxa(Fgeo,ad2,2);

    geo = sf_intalloc2(dim1,dim2);
    geogen(geo,nz,nx,ny,sou_z,sou_ox,sou_oy,sou_jx,sou_jy,sou_nx,sou_ny,rec_z,rec_nx,rec_ny,npad,noff,roll);

    sf_intwrite(geo[0],dim1*dim2,Fgeo);

    exit(0);
}

void geogen(int **geo,
            int nz, int nx, int ny,
            int sou_z, int sou_ox, int sou_oy, int sou_jx, int sou_jy, int sou_nx, int sou_ny,
            int rec_z, int rec_nx, int rec_ny,
            int npad, int noff, int roll)
/*< generate geometry file >*/
{
    int mod_oz,mod_ox,mod_oy,mod_nz,mod_nx,mod_ny,sou_x,sou_y,rec_ox,rec_oy,rec_nx_new,rec_ny_new;
    int iy,ix;

    if (roll==0) { /* fixed-spread acquisition */
        if (ny > 1) {
            mod_oz = 0;
            mod_nz = nz;
            rec_ny_new = rec_ny;
            rec_nx_new = rec_nx;
            for (iy=0; iy<sou_ny; iy++) {
                sou_y = sou_oy + iy*sou_jy;
                if (sou_y-rec_ny/2 < 0)
                    rec_oy = 0;
                else if (sou_y+rec_ny/2 > ny-1)
                    rec_oy = (ny-1)-rec_ny+1; 
                else
                    rec_oy = sou_y-rec_ny/2;
                if (rec_oy-npad < 0)
                    mod_oy = 0;
                else
                    mod_oy = rec_oy-npad;
                if (rec_oy+(rec_ny_new-1)+npad > ny-1)
                    mod_ny = (ny-1)-mod_oy+1;
                else
                    mod_ny = rec_oy+(rec_ny_new-1)+npad-mod_oy+1;
                for (ix=0; ix<sou_nx; ix++) {
                    sou_x = sou_ox + ix*sou_jx;
                    if (sou_x-rec_nx/2 < 0)
                        rec_ox = 0;
                    else if (sou_x+rec_nx/2 > nx-1)
                        rec_ox = (nx-1)-rec_nx+1;
                    else
                        rec_ox = sou_x-rec_nx/2;
                    if (rec_ox-npad < 0)
                        mod_ox = 0;
                    else
                        mod_ox = rec_ox-npad;
                    if (rec_ox+(rec_nx_new-1)+npad > nx-1)
                        mod_nx = (nx-1)-mod_ox+1;
                    else
                        mod_nx = rec_ox+(rec_nx_new-1)+npad-mod_ox+1;
                    geo[iy*sou_nx+ix][ 0] = mod_oz; geo[iy*sou_nx+ix][ 1] = mod_ox; geo[iy*sou_nx+ix][ 2] = mod_oy;
                    geo[iy*sou_nx+ix][ 3] = mod_nz; geo[iy*sou_nx+ix][ 4] = mod_nx; geo[iy*sou_nx+ix][ 5] = mod_ny;
                    geo[iy*sou_nx+ix][ 6] = sou_z ; geo[iy*sou_nx+ix][ 7] = sou_x ; geo[iy*sou_nx+ix][ 8] = sou_y ;
                    geo[iy*sou_nx+ix][ 9] = rec_z ; geo[iy*sou_nx+ix][10] = rec_ox; geo[iy*sou_nx+ix][11] = rec_oy;
                    geo[iy*sou_nx+ix][12] = rec_nx_new; 
                    geo[iy*sou_nx+ix][13] = rec_ny_new;
                }
            }
        } else {
            mod_oz =  0;
            mod_nz = nz;
            mod_oy =  0;
            mod_ny =  1;
            sou_y  =  0;
            rec_oy =  0;
            rec_nx_new = rec_nx;
            rec_ny_new = 1;
            for (ix=0; ix<sou_nx; ix++) {
                sou_x = sou_ox + ix*sou_jx;
                if (sou_x-rec_nx/2 < 0)
                    rec_ox = 0;
                else if (sou_x+rec_nx/2 > nx-1)
                    rec_ox = (nx-1)-rec_nx+1;
                else
                    rec_ox = sou_x-rec_nx/2;
                if (rec_ox-npad < 0)
                    mod_ox = 0;
                else
                    mod_ox = rec_ox-npad;
                if (rec_ox+(rec_nx_new-1)+npad > nx-1)
                    mod_nx = (nx-1)-mod_ox+1;
                else
                    mod_nx = rec_ox+(rec_nx_new-1)+npad-mod_ox+1;
                geo[ix][ 0] = mod_oz; geo[ix][ 1] = mod_ox; geo[ix][ 2] = mod_oy;
                geo[ix][ 3] = mod_nz; geo[ix][ 4] = mod_nx; geo[ix][ 5] = mod_ny;
                geo[ix][ 6] = sou_z ; geo[ix][ 7] = sou_x ; geo[ix][ 8] = sou_y ;
                geo[ix][ 9] = rec_z ; geo[ix][10] = rec_ox; geo[ix][11] = rec_oy;
                geo[ix][12] = rec_nx_new; 
                geo[ix][13] = rec_ny_new;
            }
        }
    } else { /* towed marine streamer to negative */
        if (ny > 1) {
            mod_oz = 0;
            mod_nz = nz;
            rec_ny_new = rec_ny;
            rec_nx_new = rec_nx;
            for (iy=0; iy<sou_ny; iy++) {
                sou_y = sou_oy + iy*sou_jy;
                if (sou_y-rec_ny/2 < 0)
                    rec_oy = 0;
                else if (sou_y+rec_ny/2 > ny-1)
                    rec_oy = (ny-1)-rec_ny+1;
                else
                    rec_oy = sou_y-rec_ny/2;
                if (rec_oy-npad < 0)
                    mod_oy = 0;
                else
                    mod_oy = rec_oy-npad;
                if (rec_oy+(rec_ny_new-1)+npad > ny-1)
                    mod_ny = (ny-1)-mod_oy+1;
                else
                    mod_ny = rec_oy+(rec_ny_new-1)+npad-mod_oy+1;
                for (ix=0; ix<sou_nx; ix++) {
                    sou_x = sou_ox + ix*sou_jx;
                    if (sou_x-noff-(rec_nx-1) < 0) {
                        rec_ox = 0;
                        rec_nx_new = sou_x - noff;
                    } else {
                        rec_ox = sou_x-noff-(rec_nx-1);
                        rec_nx_new = rec_nx;
                    }
                    if (rec_ox-npad/5 < 0)
                        mod_ox = 0;
                    else
                        mod_ox = rec_ox-npad/5;
                    if (rec_ox+(rec_nx_new-1)+npad > nx-1)
                        mod_nx = (nx-1)-mod_ox+1;
                    else
                        mod_nx = rec_ox+(rec_nx_new-1)+npad-mod_ox+1;
                    geo[iy*sou_nx+ix][ 0] = mod_oz; geo[iy*sou_nx+ix][ 1] = mod_ox; geo[iy*sou_nx+ix][ 2] = mod_oy;
                    geo[iy*sou_nx+ix][ 3] = mod_nz; geo[iy*sou_nx+ix][ 4] = mod_nx; geo[iy*sou_nx+ix][ 5] = mod_ny;
                    geo[iy*sou_nx+ix][ 6] = sou_z ; geo[iy*sou_nx+ix][ 7] = sou_x ; geo[iy*sou_nx+ix][ 8] = sou_y ;
                    geo[iy*sou_nx+ix][ 9] = rec_z ; geo[iy*sou_nx+ix][10] = rec_ox; geo[iy*sou_nx+ix][11] = rec_oy;
                    geo[iy*sou_nx+ix][12] = rec_nx_new; 
                    geo[iy*sou_nx+ix][13] = rec_ny_new;
                }
            }
        } else {
            mod_oz =  0;
            mod_nz = nz;
            mod_oy =  0;
            mod_ny =  1;
            sou_y  =  0;
            rec_oy =  0;
            rec_ny_new = 1;
            for (ix=0; ix<sou_nx; ix++) {
                sou_x = sou_ox + ix*sou_jx;
                if (sou_x-noff-(rec_nx-1) < 0) {
                    rec_ox = 0;
                    rec_nx_new = sou_x - noff;
                } else {
                    rec_ox = sou_x-noff-(rec_nx-1);
                    rec_nx_new = rec_nx;
                }
                if (rec_ox-npad/5 < 0)
                    mod_ox = 0;
                else
                    mod_ox = rec_ox-npad/5;
                if (rec_ox+(rec_nx_new-1)+npad > nx-1)
                    mod_nx = (nx-1)-mod_ox+1;
                else
                    mod_nx = rec_ox+(rec_nx_new-1)+npad-mod_ox+1;
                geo[ix][ 0] = mod_oz; geo[ix][ 1] = mod_ox; geo[ix][ 2] = mod_oy;
                geo[ix][ 3] = mod_nz; geo[ix][ 4] = mod_nx; geo[ix][ 5] = mod_ny;
                geo[ix][ 6] = sou_z ; geo[ix][ 7] = sou_x ; geo[ix][ 8] = sou_y ;
                geo[ix][ 9] = rec_z ; geo[ix][10] = rec_ox; geo[ix][11] = rec_oy;
                geo[ix][12] = rec_nx_new; 
                geo[ix][13] = rec_ny_new;
            }
        }
    }

}
