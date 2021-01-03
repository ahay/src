/* 2-D Fourth-order Optimized Finite-difference wave extrapolation */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
#include <math.h>
#include <limits.h>

int main(int argc, char* argv[]) 
{
    int tmp, nz, nt, iz, it, ik;
    int LEN;
    int *s1;
    float dt, dz;
    sf_complex *nxt, *old, *cur, *wav;
    sf_complex **coef, plap;
    sf_file out, ini, G;
    int mode;
    bool cpxexp;

    sf_init(argc,argv);
    out = sf_output("out");
    ini = sf_input("in");   /* initial condition wavlet*/
    G = sf_input("G");   /* coefficient matrix */

    if (SF_COMPLEX != sf_gettype(ini)) sf_error("Need complex input");
    if (!sf_histint(ini,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histfloat(ini,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("mode",&mode)) mode=0;
    if (!sf_getbool("cpxexp",&cpxexp)) cpxexp=true;
    if (!sf_histint(G,"n1",&tmp)) sf_error("No n1= in input");
    if (!sf_histint(G,"n2",&LEN)) sf_error("No n2= in input");
    if (nz != tmp) sf_error("Mismatch!");

    sf_putint(out,"n1",nz);
    sf_putfloat(out,"d1",dz);
    sf_putint(out,"n2",nt);
    sf_putfloat(out,"d2",dt);
    sf_putfloat(out,"o2",0.0); 
    sf_settype(out,SF_COMPLEX);

    sf_warning("nt=%d",nt);
    sf_warning("dt=%g",dt);
    wav =  sf_complexalloc(nz);
    sf_complexread(wav,nz,ini);

    old  = sf_complexalloc(nz);
    cur  = sf_complexalloc(nz);
    nxt  = sf_complexalloc(nz);

    coef = sf_complexalloc2(nz,LEN);
    
    sf_complexread(coef[0],nz*LEN,G);

    s1 = sf_intalloc(LEN);

    if (cpxexp) {
        //for (iz=0; iz<LEN/2; iz++) s1[iz]= iz - (LEN/2-1);
        //for (iz=LEN/2; iz<LEN; iz++) s1[iz]= iz - LEN/2;
        //for (iz=0; iz<LEN; iz++) s1[iz]= (iz%2==0) ? iz/2 : -iz/2;
        for (iz=0; iz<LEN; iz++) s1[iz]= (iz%2==0) ? iz/2-(LEN/2-1) : -(iz/2-(LEN/2-1));
    } else {
        for (iz=0; iz<LEN; iz++) s1[iz]= iz;
    }
    for (iz=0; iz<LEN; iz++) sf_warning("s1=%d",s1[iz]); 
    
    for (iz=0; iz < nz; iz++) {
        cur[iz] = wav[iz];
        old[iz] = sf_cmplx(0,0);
        nxt[iz] = sf_cmplx(0,0);
    }

    /* propagation in time */
    for (it=0; it < nt; it++) {
        sf_warning("it=%d;",it);

        sf_complexwrite(cur,nz,out);

        //for (iz=LEN2-1; iz < nz-(LEN2-1); iz++) {  
        for (iz=0; iz < nz; iz++) {  

            // apply lowrank FD coefficients
            plap = sf_cmplx(0,0);
            for (ik=0; ik < LEN; ik++) {
#ifdef SF_HAS_COMPLEX_H
                if (cpxexp)
                    plap += coef[ik][iz]*cur[(iz+s1[ik]+nz)%nz];
                else
                    plap += coef[ik][iz]*(cur[(iz+s1[ik]+nz)%nz]+cur[(iz-s1[ik]+nz)%nz]);
#else
                if (cpxexp)
                    plap = sf_cadd(plap,sf_cmul(coef[ik][iz],cur[(iz+s1[ik]+nz)%nz]));
                else
                    plap = sf_cadd(plap,sf_cmul(coef[ik][iz],sf_cadd(cur[(iz+s1[ik]+nz)%nz],cur[(iz-s1[ik]+nz)%nz])));
#endif
            }

            if (mode == 0) {
                nxt[iz] = plap;
            } else if (mode == 1) {
#ifdef SF_HAS_COMPLEX_H
                nxt[iz] = plap + old[iz]; 
#else
                nxt[iz] = sf_cadd(plap,old[iz]); 
#endif
            } else {
#ifdef SF_HAS_COMPLEX_H
                nxt[iz] = plap - old[iz]; 
#else
                nxt[iz] = sf_csub(plap,old[iz]); 
#endif 
            }

        }

        //for (iz=LEN2-1; iz < nz-(LEN2-1); iz++) {  
        for (iz=0; iz < nz; iz++) {  
            old[iz] = cur[iz];
            cur[iz] = nxt[iz];
        }

    }
    sf_warning(".");

    exit(0); 
} 
