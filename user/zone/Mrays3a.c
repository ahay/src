/* Ray tracing in general anisotropic media by a Runge-Kutta integrator in 3-D.

Takes: > rays.rsf

Rays can be plotted with sfplotrays.
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

#include "araytracegenani.h"

int main(int argc, char* argv[])
{
    bool escvar, wantphase, wantpolarization, wantsine;
    char *raytype, *velfile;
    int is, n[3], nm, order, nshot, ndim, two, na, nb;
    int nt, nt1, nr, ir, it, i, ia, ib;
    float t, dt, da=0., a0, amax, db=0., b0, bmax;
    float x[3], p[3], d[3], o[3], **traj, **s, *a, *b, *sine, **polarization;
    float *C11, *C12, *C13, *C14, *C15, *C16, *C22, *C23, *C24, *C25, *C26, *C33, *C34, *C35, *C36, *C44, *C45, *C46, *C55, *C56, *C66, *medthe, *medphi;
    araytrace rt;
    sf_file shots, rays, angles, mediumtheta, mediumphi, sinenu, phase, polar, tt, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, c55, c56, c66;

    sf_init (argc,argv);
    rays = sf_output("out");
    tt = sf_output("tt");
    c11 = sf_input("in");
    c12 = sf_input("c12");
    c13 = sf_input("c13");
    c22 = sf_input("c22");
    c23 = sf_input("c23");
    c33 = sf_input("c33");
    c44 = sf_input("c44");
    c55 = sf_input("c55");
    c66 = sf_input("c66");

    /* get 3-D grid parameters */
    if (!sf_histint(c11,"n1",n))     sf_error("No n1= in input");
    if (!sf_histint(c11,"n2",n+1))   sf_error("No n2= in input");
    if (!sf_histint(c11,"n3",n+2))   sf_error("No n3= in input");
    if (!sf_histfloat(c11,"d1",d))   sf_error("No d1= in input");
    if (!sf_histfloat(c11,"d2",d+1)) sf_error("No d2= in input");
    if (!sf_histfloat(c11,"d3",d+2)) sf_error("No d3= in input");
    if (!sf_histfloat(c11,"o1",o))   o[0]=0.;
    if (!sf_histfloat(c11,"o2",o+1)) o[1]=0.;
    if (!sf_histfloat(c11,"o3",o+2)) o[2]=0.;

    /* additional parameters */
    if(NULL == (raytype = sf_getstring("raytype"))) raytype='p';
    /* If p = qP, s = s1, and t = s2 */
    
    if(!sf_getint("order",&order)) order=4;
    /* Interpolation order */

    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    /* Number of time steps */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
    /* Sampling in time */

    if(!sf_getbool("escvar",&escvar)) escvar=false;
    /* If y - output escape values, n - trajectories */
    
     if(!sf_getbool("wantphase",&wantphase)) wantphase=false;
    /* If y - output phase direction */

     if(!sf_getbool("wantpolarization",&wantpolarization)) wantpolarization=false;
    /* If y - output polarization direction */
    
     if(!sf_getbool("wantsine",&wantsine)) wantsine=false;
    /* If y - output proximity to singularity */

    /* get shot locations */
    if (NULL != sf_getstring("shotfile")) {
        /* file with shot locations */
        shots = sf_input("shotfile");
        if (!sf_histint(shots,"n1",&ndim) || 3 != ndim) 
            sf_error("Must have n1=2 in shotfile");
        if (!sf_histint(shots,"n2",&nshot)) 
            sf_error("No n2= in shotfile");
  
        s = sf_floatalloc2 (ndim,nshot);
        sf_floatread(s[0],ndim*nshot,shots);
        sf_fileclose (shots);
    } else {
        nshot = 1;
        ndim = 3;

        s = sf_floatalloc2 (ndim,nshot);

        if (!sf_getfloat("zshot",&s[0][0])) s[0][0]=o[0];
        /* shot location in depth (if shotfile is not specified) */
        if (!sf_getfloat("yshot",&s[0][1])) s[0][1]=o[1] + 0.5*(n[1]-1)*d[1];
        /* shot location inline (if shotfile is not specified) */
        if (!sf_getfloat("xshot",&s[0][2])) s[0][2]=o[2] + 0.5*(n[2]-1)*d[2];
        /* shot location crossline (if shotfile is not specified) */

        sf_warning("Shooting from z=%g, y=%g, x=%g",s[0][0],s[0][1],s[0][2]);
    }


    if (NULL != sf_getstring("anglefile")) {
        /* file with initial angles */
        angles = sf_input("anglefile");

        if (!sf_histint(angles,"n1",&nr)) sf_error("No n1= in anglefile");
        if (!sf_histint(angles,"n2",&two) || 2 != two) sf_error("Need n2=2 in anglefile");

        a = sf_floatalloc(nr);
        b = sf_floatalloc(nr);
    } else {
        angles = NULL;

        if (!sf_getint("na",&na)) sf_error("Need na=");
        /* Number of azimuths (if anglefile is not specified) */
        if (!sf_getint("nb",&nb)) sf_error("Need nb=");
        /* Number of inclinations (if anglefile is not specified) */

        if (!sf_getfloat("a0",&a0)) a0 = 0.; 
        /* First azimuth angle in degrees (if anglefile is not specified) */
        if (!sf_getfloat("amax",&amax)) amax=360.;
        /* Maximum azimuth angle in degrees (if anglefile is not specified) */

        if (!sf_getfloat("b0",&b0)) b0 = 0.; 
        /* First inclination angle in degrees (if anglefile is not specified) */
        if (!sf_getfloat("bmax",&bmax)) bmax=180.;
        /* Maximum inclination angle in degrees (if anglefile is not specified) */
        
        /* convert degrees to radians */
        a0 *= SF_PI/180.;
        amax *= SF_PI/180.;
        b0 *= SF_PI/180.;
        bmax *= SF_PI/180.;

        /* figure out angle spacing */
        da = (na > 1)? (amax - a0)/(na-1) : 0.;
        db = (nb > 1)? (bmax - b0)/(nb-1) : 0.;

        nr = na*nb;

        a = sf_floatalloc(nr);
        b = sf_floatalloc(nr);

        for (ir=ib=0; ib < nb; ib++) {
            for (ia=0; ia < na; ia++, ir++) {
                b[ir] = b0 + ib*db;
                a[ir] = a0 + ia*da;
            }
        }
    }
    
    /* specify output dimensions */
    nt1 = nt+1;
    if (escvar) {
        sf_putint (rays, "n1", 4);
        sf_putfloat (rays, "o1", 0.0);
        sf_putfloat (rays, "d1", 1.0);
        sf_putstring (rays, "label1", "Escape variable");
        sf_putstring (rays, "unit1", "");
        if (NULL == angles) {
            sf_putint (rays, "n2", na);
            sf_putfloat (rays, "d2", da*180.0/SF_PI);
            sf_putfloat (rays, "o2", a0*180.0/SF_PI);
            sf_putstring (rays, "label2", "Azimuth");
            sf_putstring (rays, "unit2", "Degrees");
            sf_putint (rays, "n3", nb);
            sf_putfloat (rays, "d3", db*180.0/SF_PI);
            sf_putfloat (rays, "o3", b0*180.0/SF_PI);
            sf_putstring (rays, "label3", "Inclination");
            sf_putstring (rays, "unit3", "Degrees");
            sf_putint (rays,"n4",nshot);
            sf_putfloat (rays, "d4", 1.0);
            sf_putfloat (rays, "o4", 0.0);
            sf_putstring (rays, "label4", "Shots");
            sf_putstring (rays, "unit4", "");
        } else {
            sf_putint (rays,"n2",nr);
            sf_putfloat (rays, "d2", 1.0);
            sf_putfloat (rays, "o2", 0.0);
            sf_putstring (rays, "label2", "Angles");
            sf_putstring (rays, "unit2", "");
            sf_putint (rays,"n3",nshot);
            sf_putfloat (rays, "d3", 1.0);
            sf_putfloat (rays, "o3", 0.0);
            sf_putstring (rays, "label3", "Shots");
            sf_putstring (rays, "unit3", "");
        }
    } else {
        sf_putint (rays,"n1",ndim);
        sf_putfloat (rays, "o1", 0.0);
        sf_putfloat (rays, "d1", 1.0);
        sf_putstring (rays, "label1", "Dimension");
        sf_putstring (rays, "unit1", "");
        
        sf_putint (rays,"n2",nt1);
        sf_putfloat (rays, "o2", 0.0);
        sf_putfloat (rays, "d2", dt);
        sf_putstring (rays, "label2", "Time");
        sf_putstring (rays, "unit2", "s");
        
        sf_putint (rays,"n3",nr);
        sf_putstring (rays, "label3", "Angles"); /*Including both azimuths and inclinations*/
        sf_putstring (rays, "unit3", "Radians");
        
        sf_putint (rays,"n4",nshot);
        sf_putstring (rays, "label4", "Shots");
        sf_putstring (rays, "unit4", "");
        
        sf_putint(tt,"n1",nr*nshot);
        sf_putfloat (tt, "o1", 0.0);
        sf_putfloat (tt, "d1", 1.0);
        sf_putstring (tt, "label1", "Traveltime");
        sf_putstring (tt, "unit1", "s");
        sf_putint(tt,"n2",1);
        sf_putint(tt,"n3",1);
        
        if (raytype[0]!='p'){
            sinenu = sf_output("sine");
            sf_putint(sinenu,"n1",nt1);
            sf_putfloat (sinenu, "o1", 0.0);
            sf_putfloat (sinenu, "d1", 1.0);
            sf_putstring (sinenu, "unit1", "");
            sf_putint(sinenu,"n2",nr);
            sf_putfloat (sinenu, "d2", 1.0);
            sf_putfloat (sinenu, "o2", 0.0);
            sf_putstring (sinenu, "unit2", "");
            sf_putint(sinenu,"n3",1);
        }
        
        if (wantphase) {
            phase = sf_output("phase");
            sf_putint (phase,"n1",ndim);
            sf_putfloat (phase, "o1", 0.0);
            sf_putfloat (phase, "d1", 1.0);
            sf_putstring (phase, "label1", "Dimension");
            sf_putstring (phase, "unit1", "");
            
            sf_putint (phase,"n2",nt1);
            sf_putfloat (phase, "o2", 0.0);
            sf_putfloat (phase, "d2", dt);
            sf_putstring (phase, "label2", "Time");
            sf_putstring (phase, "unit2", "s");
            
            sf_putint (phase,"n3",nr);
            sf_putstring (phase, "label3", "Angles"); /*Including both azimuths and inclinations*/
            sf_putstring (phase, "unit3", "Radians");
            
            sf_putint (phase,"n4",nshot);
            sf_putstring (phase, "label4", "Shots");
            sf_putstring (phase, "unit4", "");
        }
        if (wantpolarization) {
            polar = sf_output("polarization");
            sf_putint (polar,"n1",ndim);
            sf_putfloat (polar, "o1", 0.0);
            sf_putfloat (polar, "d1", 1.0);
            sf_putstring (polar, "label1", "Dimension");
            sf_putstring (polar, "unit1", "");
            
            sf_putint (polar,"n2",nt1);
            sf_putfloat (polar, "o2", 0.0);
            sf_putfloat (polar, "d2", dt);
            sf_putstring (polar, "label2", "Time");
            sf_putstring (polar, "unit2", "s");
            
            sf_putint (polar,"n3",nr);
            sf_putstring (polar, "label3", "Angles"); /*Including both azimuths and inclinations*/
            sf_putstring (polar, "unit3", "Radians");
            
            sf_putint (polar,"n4",nshot);
            sf_putstring (polar, "label4", "Shots");
            sf_putstring (polar, "unit4", "");
        }
    }
            
    /* get all cij (default = orthorhombic) ---------------------------------------------------------*/
    nm = n[0]*n[1]*n[2];
    C11 = sf_floatalloc(nm); C12 = sf_floatalloc(nm); C13 = sf_floatalloc(nm);
    C22 = sf_floatalloc(nm); C23 = sf_floatalloc(nm); C33 = sf_floatalloc(nm);
    C44 = sf_floatalloc(nm); C55 = sf_floatalloc(nm); C66 = sf_floatalloc(nm);

    sf_floatread(C11,nm,c11); sf_floatread(C12,nm,c12); sf_floatread(C13,nm,c13); 
    sf_floatread(C22,nm,c22); sf_floatread(C23,nm,c23); sf_floatread(C33,nm,c33); 
    sf_floatread(C44,nm,c44); sf_floatread(C55,nm,c55); sf_floatread(C66,nm,c66);

    if (NULL != (velfile = sf_getstring("mediumtheta"))) {
        mediumtheta = sf_input("mediumtheta");
        medthe = sf_floatalloc(nm);
        sf_floatread(medthe,nm,mediumtheta);

        free(velfile);
        sf_fileclose(mediumtheta);
    } else {
        medthe = NULL;
    }
    
    if (NULL != (velfile = sf_getstring("mediumphi"))) {
        mediumphi = sf_input("mediumphi");
        medphi = sf_floatalloc(nm);
        sf_floatread(medphi,nm,mediumphi);

        free(velfile);
        sf_fileclose(mediumphi);
    } else {
        medphi = NULL;
    }

    if (NULL != (velfile = sf_getstring("c14"))) {
        c14 = sf_input(velfile);
        C14 = sf_floatalloc(nm);
        sf_floatread(C14,nm,c14);

        free(velfile);
        sf_fileclose(c14);
    } else {
        C14 = NULL;
    }
    if (NULL != (velfile = sf_getstring("c15"))) {
        c15 = sf_input(velfile);
        C15 = sf_floatalloc(nm);
        sf_floatread(C15,nm,c15);

        free(velfile);
        sf_fileclose(c15);
    } else {
        C15 = NULL;
    }
    if (NULL != (velfile = sf_getstring("c16"))) {
        c16 = sf_input(velfile);
        C16 = sf_floatalloc(nm);
        sf_floatread(C16,nm,c16);

        free(velfile);
        sf_fileclose(c16);
    } else {
        C16 = NULL;
    }
    if (NULL != (velfile = sf_getstring("c24"))) {
        c24 = sf_input(velfile);
        C24 = sf_floatalloc(nm);
        sf_floatread(C24,nm,c24);

        free(velfile);
        sf_fileclose(c24);
    } else {
        C24 = NULL;
    }
    if (NULL != (velfile = sf_getstring("c25"))) {
        c25 = sf_input(velfile);
        C25 = sf_floatalloc(nm);
        sf_floatread(C25,nm,c25);

        free(velfile);
        sf_fileclose(c25);
    } else {
        C25 = NULL;
    }
    if (NULL != (velfile = sf_getstring("c26"))) {
        c26 = sf_input(velfile);
        C26 = sf_floatalloc(nm);
        sf_floatread(C26,nm,c26);

        free(velfile);
        sf_fileclose(c26);
    } else {
        C26 = NULL;
    }
    if (NULL != (velfile = sf_getstring("c34"))) {
        c34 = sf_input(velfile);
        C34 = sf_floatalloc(nm);
        sf_floatread(C34,nm,c34);

        free(velfile);
        sf_fileclose(c34);
    } else {
        C34 = NULL;
    }
    if (NULL != (velfile = sf_getstring("c35"))) {
        c35 = sf_input(velfile);
        C35 = sf_floatalloc(nm);
        sf_floatread(C35,nm,c35);

        free(velfile);
        sf_fileclose(c35);
    } else {
        C35 = NULL;
    }
    if (NULL != (velfile = sf_getstring("c36"))) {
        c36 = sf_input(velfile);
        C36 = sf_floatalloc(nm);
        sf_floatread(C36,nm,c36);

        free(velfile);
        sf_fileclose(c36);
    } else {
        C36 = NULL;
    }
    if (NULL != (velfile = sf_getstring("c45"))) {
        c45 = sf_input(velfile);
        C45 = sf_floatalloc(nm);
        sf_floatread(C45,nm,c45);

        free(velfile);
        sf_fileclose(c45);
    } else {
        C45 = NULL;
    }
    if (NULL != (velfile = sf_getstring("c46"))) {
        c46 = sf_input(velfile);
        C46 = sf_floatalloc(nm);
        sf_floatread(C46,nm,c46);

        free(velfile);
        sf_fileclose(c46);
    } else {
        C46 = NULL;
    }
    if (NULL != (velfile = sf_getstring("c56"))) {
        c56 = sf_input(velfile);
        C56 = sf_floatalloc(nm);
        sf_floatread(C56,nm,c56);

        free(velfile);
        sf_fileclose(c56);
    } else {
        C56 = NULL;
    }

    /* initialize ray tracing object ----------------------------------------------------------------- */
    rt = araytracegenani_init (3, nt, dt, n, o, d, 
                         C11, C12, C13, C14, C15,
                         C16, C22, C23, C24, C25,
                         C26, C33, C34, C35, C36,
                         C44, C45, C46, C55, C56,
                         C66, medthe, medphi, order,
                         wantsine,wantpolarization);
    
    free (C11); free (C12); free (C13); 
    free (C22); free (C23); free (C33); 
    free (C44); free (C55); free (C66);
    
    if (NULL != medthe) free (medthe); if (NULL != medphi) free(medphi);
    if (NULL != C14) free(C14); if (NULL != C15) free(C15);
    if (NULL != C16) free(C16); if (NULL != C24) free(C24); 
    if (NULL != C25) free(C25); if (NULL != C34) free(C34); 
    if (NULL != C35) free(C35); if (NULL != C36) free(C36); 
    if (NULL != C36) free(C36); if (NULL != C45) free(C45); 
    if (NULL != C46) free(C46); if (NULL != C56) free(C56);
    
    traj = sf_floatalloc2 (2*ndim,nt1);
    polarization = sf_floatalloc2 (ndim,nt1);
    sine = sf_floatalloc (nt1);

    for( is = 0; is < nshot; is++) { /* loop over shots */
        /* initialize angles */
        if (NULL != angles) {
            sf_floatread(a,nr,angles);
            sf_floatread(b,nr,angles);
        } 
        
        for (ir = 0; ir < nr; ir++) { /* loop over rays */
            /* initialize position */
            x[0] = s[is][0]; 
            x[1] = s[is][1];
            x[2] = s[is][2];

            /* initialize direction */
            p[2] = cosf(a[ir])*sinf(b[ir]); // x
            p[1] = sinf(a[ir])*sinf(b[ir]); // y
            p[0] = cosf(b[ir]); // z
            
            /* Which raytype to trace */
            switch(raytype[0]) {
                case 'p': // qP
                    it = trace_p_aray (rt, x, p, traj, polarization);
                    break;
                case 's': // qS1
                    it = trace_s1_aray (rt, x, p, traj, sine, polarization);
                    break;
                case 't': // qS2
                    it = trace_s2_aray (rt, x, p, traj, sine, polarization);
                    break;
                case 'c': // qS1 follow most coupled via dot product
                    it = trace_s1c_aray (rt, x, p, traj, sine, polarization);
                    break;
                case 'k': // qS2 follow most couple via dot product
                    it = trace_s2c_aray (rt, x, p, traj, sine, polarization);
                    break;
            }
            
            if (it < 0) it = -it; /* keep side-exiting rays */
            

            if (escvar) {
                /* Write escape variables only */
                sf_floatwrite (traj[it],ndim,rays); /* z, y, x */
                if (wantphase) sf_floatwrite (traj[it]+ndim,ndim,phase); /* pz, py, px */
                if (wantpolarization) sf_floatwrite (polarization[it],ndim,polar); /* polarz, polary polarx */
                t = it*dt; /* t */
                sf_floatwrite (&t, 1, rays);
            } else {
                for (i=0; i < nt1; i++) {
                    if (0==it || it > i) {
                        sf_floatwrite (traj[i],ndim,rays);
                        if (wantphase) sf_floatwrite (traj[i]+ndim,ndim,phase);
                        if (wantpolarization) sf_floatwrite (polarization[i],ndim,polar); 
                    } else {
                        sf_floatwrite (traj[it],ndim,rays);
                        if (wantphase) sf_floatwrite (traj[it]+ndim,ndim,phase);
                        if (wantpolarization) sf_floatwrite (polarization[i],ndim,polar); 
                    }
                }
                if (raytype[0]!='p' && wantsine) sf_floatwrite (sine,nt1,sinenu);
                t = it*dt; /* t */
                sf_floatwrite (&t, 1, tt);
            }
        }
    }

    exit (0);
}

/*         $Id$         */
