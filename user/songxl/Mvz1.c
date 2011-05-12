/* constant-velocity prestack exploditing reflector. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "cosft3.h"
#include "abcpass.h"


int main(int argc, char* argv[])
{
    bool mig;
    int it, nt, ix, nx, iz, nz, ih, nh, it1, it2, its, snap;
    int nzb, nbt, nbb;
    float ct, cb;
    float dt, dx, dz, dh, kx, kz, kh,  h, x, c, dkx, dkz, dkh, v0, pi=SF_PI;
    float ***prev, ***curr, **img, **dat, ***tmp, ***tmpq, *v;
    float ***a, ***b3, ***d3, epsilon;
    sf_file data, image, snaps, vel; 
    int nvz;
    float dvz, vz;

    sf_init(argc,argv);

    if (!sf_getbool("mig",&mig)) mig=false;
    /* if n, modeling; if y, migration */

    if (!sf_getint("snap",&snap)) snap=0;
    /* interval for snapshots */
    
    snaps = (snap > 0)? sf_output("snaps"): NULL;
    /* (optional) snapshot file */
    if (!sf_getfloat("epsilon",&epsilon)) epsilon=0.00001;
    vel = sf_input("vel");   /* velocity */

    if (mig) { /* migration */
	data = sf_input("in");
	image = sf_output("out");

	if (!sf_histint(data,"n1",&nh));
	if (!sf_histfloat(data,"d1",&dh));

	if (!sf_histint(data,"n2",&nx));
	if (!sf_histfloat(data,"d2",&dx));

	if (!sf_histint(data,"n3",&nt));
	if (!sf_histfloat(data,"d3",&dt));

	if (!sf_getint("nz",&nz)) sf_error("Need nz=");
	/* time samples (if migration) */
	if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
	/* time sampling (if migration) */

	sf_putint(image,"n1",nz);
	sf_putfloat(image,"d1",dz);
	sf_putstring(image,"label1","Depth");

	sf_putint(image,"n3",1); /* stack for now */
    } else { /* modeling */
	image = sf_input("in");
	data = sf_output("out");

	if (!sf_histint(image,"n1",&nz));
	if (!sf_histfloat(image,"d1",&dz));

	if (!sf_histint(image,"n2",&nx));
	if (!sf_histfloat(image,"d2",&dx));

	if (!sf_getint("nt",&nt)) sf_error("Need nt=");
	/* time samples (if modeling) */
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
	/* time sampling (if modeling) */

	if (!sf_getint("nh",&nh)) sf_error("Need nh=");
        /* offset samples (if modeling) */
	if (!sf_getfloat("dh",&dh)) sf_error("Need dh=");
	/* offset sampling (if modeling) */

	sf_putint(data,"n1",nh);
	sf_putfloat(data,"d1",dh);
	sf_putstring(data,"label1","Half-Offset");

	sf_putint(data,"n3",nt);
	sf_putfloat(data,"d3",dt);
	sf_putstring(data,"label3","Time");
	sf_putstring(data,"unit3","s");
    }

    dat = sf_floatalloc2(nh,nx);
    if (!sf_histint(vel,"n1",&nvz));
    if (nvz!=nz) sf_error("nvz and nz do not match!");
    if (!sf_histfloat(vel,"d1",&dvz));
    if (dvz!=dz) sf_error("dvz and dz do not match!");

    if (!sf_getint("nbt",&nbt)) nbt=44;
    if (!sf_getint("nbb",&nbb)) nbb=44;

    if (!sf_getfloat("ct",&ct)) ct = 0.001; /*decaying parameter*/
    if (!sf_getfloat("cb",&cb)) cb = 0.001; /*decaying parameter*/

    nzb = nz + nbt + nbb;

    bdz_init(nx,nz,nh,nbt,nbb,0,0,ct,cb,0,0);
    img = sf_floatalloc2(nzb,nx);


    prev = sf_floatalloc3(nh,nx,nzb);
    curr = sf_floatalloc3(nh,nx,nzb);
    tmp  = sf_floatalloc3(nh,nx,nzb);
    tmpq  = sf_floatalloc3(nh,nx,nzb);


    v = sf_floatalloc(nzb);
    sf_floatread(v+nbt,nz,vel); 
    for (iz=0; iz<nbt; iz++) v[iz]=v[nbt];
    for (iz=nz+nbt; iz<nzb; iz++) v[iz]=v[nz+nbt-1]; 
    v0 =0.0;
    for (iz=0; iz < nzb; iz++) {
            v0 += v[iz]*v[iz];
    }
    v0 = sqrtf(v0/(nzb));
    /* velocity */


/*
    dx = cosft_dk(nx,dx);
    dz = cosft_dk(nz,dz);
    dh = cosft_dk(nh,dh);
*/

    dkx = cosft_dk(nx,dx);
    dkz = cosft_dk(nzb,dz);
    dkh = cosft_dk(nh,dh);

    a = sf_floatalloc3(nh,nx,nzb);
    b3 = sf_floatalloc3(nh,nx,nzb);
    d3 = sf_floatalloc3(nh,nx,nzb);
    for (iz=0; iz < nzb; iz++) {
        vz = v[iz];
        for (ih=0; ih < nh; ih++) {b3[iz][0][ih]=0.;d3[iz][0][ih]=0.;a[iz][0][ih]=0.;}
        for (ix=1; ix < nx; ix++) {
	    kx = ix*dkx*2.0*pi;
            b3[iz][ix][0]=0.;d3[iz][ix][0]=0.;a[iz][ix][0]=0.;
            //if(ix < 0.001*epsilon) kx = epsilon;
            for (ih=1; ih < nh; ih++) {
		kh = ih*dkh*2.0*pi;
              //  if(ih < 0.001*epsilon) kh = epsilon;
                b3[iz][ix][ih] = dt*cosf(sqrtf(kx*kh)*dz)/(sinf(sqrtf(kx*kh)*dz)*sinf(sqrtf(kx*kh)*dz)*sinf((kx+kh)*dt*v0/4.0)*sinf((kx+kh)*dt*v0/4.0)*2.0*(kx+kh)*dz*dz)
*(2.0*v0*cosf((kx+kh)*dt*v0/4.0)*sinf((kx+kh)*dt*vz/4.0)*sinf((kx+kh)*dt*vz/4.0)/sinf((kx+kh)*dt*v0/4.0)-vz*sinf((kx+kh)*dt*vz/2.0));
                d3[iz][ix][ih] = dt/(sinf(sqrtf(kx*kh)*dz)*sinf(sqrtf(kx*kh)*dz)*sinf((kx+kh)*dt*v0/4.0)*sinf((kx+kh)*dt*v0/4.0)*8.0*(kx+kh)*dz*dz)
*(-2.0*v0*cosf((kx+kh)*dt*v0/4.0)*sinf((kx+kh)*dt*vz/4.0)*sinf((kx+kh)*dt*vz/4.0)/sinf((kx+kh)*dt*v0/4.0)+vz*sinf((kx+kh)*dt*vz/2.0));
                a[iz][ix][ih] = sinf(((kh+kx)*dt*vz)/4.0)*sinf(((kh+kx)*dt*vz)/4.0)/(sinf(((kh+kx)*dt*v0)/4.0)*sinf(((kh+kx)*dt*v0)/4.0))-2.0*b3[iz][ix][ih]*cosf(sqrtf(kx*kh)*dz)-2.0*d3[iz][ix][ih]*cosf(2.0*sqrtf(kx*kh)*dz); 
            } 
        } 
    }
   

   if (NULL != snaps) {
	sf_putint(snaps,"n1",nh);

	sf_putfloat(snaps,"d1",dkh);
	sf_putstring(snaps,"label1","Half-Offset");

	sf_putint(snaps,"n2",nx);
	sf_putfloat(snaps,"d2",dkx);
	sf_putstring(snaps,"label2","Midpoint");

	sf_putint(snaps,"n3",nz);
	sf_putfloat(snaps,"d3",dz);
	sf_putstring(snaps,"label3","Depth");

	sf_putint(snaps,"n4",nt/snap);
	sf_putfloat(snaps,"d4",dt*snap);
	sf_putfloat(snaps,"o4",0.);
	sf_putstring(snaps,"label4","Time");
    }

    for (iz=0; iz < nzb; iz++) {
	for (ix=0; ix < nx; ix++) {
	    for (ih=0; ih < nh; ih++) {
		prev[iz][ix][ih] = 0.;
		curr[iz][ix][ih] = 0.;
	    }
	}
    }

    if (mig) { /* migration */
	/* initialize image */
	for (iz=0; iz < nzb; iz++) {
	    for (ix=0; ix < nx; ix++) {
		img[ix][iz] = 0.;
	    }
	}

	/* step backward in time */
	it1 = nt-1;
	it2 = -1;
	its = -1;

    } else { /* modeling */
	//sf_floatread(img[nbt],nz*nx,image);
	for (ix=0; ix < nx; ix++) {
            sf_floatread(img[ix]+nbt,nz,image);
	}
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nbt; iz++) {
		img[ix][iz] = img[ix][nbt];
	    }
	    for (iz=0; iz < nbb; iz++) {
		img[ix][nz+nbt+iz] = img[ix][nz+nbt-1];
	    }
	}

	/* Initialize model */
	//for (iz=1; iz < nz; iz++) {
	for (iz=0; iz < nzb; iz++) {
	    for (ix=0; ix < nx; ix++) {
		curr[iz][ix][0] = img[ix][iz];
	    }
	}
	//cosft3(false,nh,nx,nzb,curr);
	cosfthx(false,nh,nx,nzb,curr);
	
	/* step forward in time */
	it1 = 0;
	it2 = nt;
	its = +1;
    }


    /* time stepping */
    for (it=it1; it != it2; it += its) {
	sf_warning("it=%d;",it);

	if (mig) { /* migration <- read data */
	    sf_floatread(dat[0],nh*nx,data);
	    cosft2(false,nh,nx,dat);	    
	} else {
	    for (ix=0; ix < nx; ix++) {
		for (ih=0; ih < nh; ih++) {
		    dat[ix][ih] = 0.;
		}
	    }
	}

	if (NULL != snaps && 0 == it%snap) 
	    sf_floatwrite(curr[nbt][0],nh*nx*nz,snaps);

	for (ix=0; ix < nx; ix++) {
            for (ih=0; ih < nh; ih++) {
		if (mig) {
                    curr[nbt][ix][ih] += dat[ix][ih];
	//c += (iz==nz-1 || iz==0)? dat[ix][ih]*0.5: dat[ix][ih];
		    } else {
			//dat[ix][ih] += (iz==nz-1 || iz==0)? c*0.5: c;
			dat[ix][ih] = curr[nbt][ix][ih];
		    }
            }
        }
	for (iz=0; iz < nzb; iz++) {
	    for (ix=0; ix < nx; ix++) {
		for (ih=0; ih < nh; ih++) {
                    tmp[iz][ix][ih] = curr[iz][ix][ih];
                }
            }
        }
           
        //cosftz(false,nh,nx,nzb,curr);
        cosftz(false,nh,nx,nzb,tmp);
	//for (iz=1; iz < nz; iz++) {
	for (iz=0; iz < nzb; iz++) {
	    kz = iz*dkz*2.0*pi;
	    for (ix=0; ix < nx; ix++) {
		kx = ix*dkx*2.0*pi;
		//x = (kz*kz+kx*kx)/kz;
                   
	       if(kz < epsilon) 
		 x = ((kz+epsilon)*(kz+epsilon)+kx*kx)/(kz+epsilon);
               else 
		 x = (kz*kz+kx*kx)/(kz);

		for (ih=0; ih < nh; ih++) {
		    kh = ih*dkh*2.0*pi;
                    //if(kz < sqrtf(kx*kh)) { ih = nh; continue;}
                    if(kz < sqrtf(kx*kh)) { tmp[iz][ix][ih] = 0.0; continue;}
		    //h = (kz*kz+kh*kh)/kz;
	            if(kz < epsilon) 
		      h  = ((kz+epsilon)*(kz+epsilon)+kh*kh)/(kz+epsilon);
                    else 
		      h  = (kz*kz+kh*kh)/(kz);
		    
		    //c = curr[iz][ix][ih];
		    c = tmp[iz][ix][ih];

		    /*curr[iz][ix][ih] = 2*cosf(v*sqrtf(x*h))*c - prev[iz][ix][ih];*/
                   // if(kz < sqrtf(kx*kh)) { ih = nh; continue;}
	              //if(kz < sqrtf(kx*kh))
                       // { ih = nh; continue;}

                       // tmp[iz][ix][ih] = 1.0/((kz+epsilon)*(kz+epsilon))*c;
                       // tmp[iz][ix][ih] = 0;
                      //else{ 
		      tmp[iz][ix][ih] = 2*(cosf(v0*dt*sqrtf(x*h)/2.0)-1.0)*c;
                      //}
		    //prev[iz][ix][ih] = c;
		}
	    }
	}

	cosftz(true,nh,nx,nzb,tmp);
	for (iz=2; iz < nzb-2; iz++) {
	    for (ix=0; ix < nx; ix++) {
		for (ih=0; ih < nh; ih++) {
                    tmpq[iz][ix][ih] = a[iz][ix][ih]*tmp[iz][ix][ih] + b3[iz][ix][ih]*(tmp[iz+1][ix][ih]+tmp[iz-1][ix][ih]) + d3[iz][ix][ih]*(tmp[iz+2][ix][ih]+tmp[iz-2][ix][ih]);
                }
            }
        }
	for (ix=0; ix < nx; ix++) {
            for (ih=0; ih < nh; ih++) {
                /*
                tmpq[0][ix][ih] = a*tmp[0][ix][ih] + 1.0*b3*(tmp[1][ix][ih]) + 1.0*d3*(tmp[2][ix][ih]);
                tmpq[1][ix][ih] = a*tmp[1][ix][ih] + b3*(tmp[0][ix][ih]+tmp[2][ix][ih]) + 1.0*d3*(tmp[3][ix][ih]);
                tmpq[nzb-2][ix][ih] = a*tmp[nzb-2][ix][ih] + b3*(tmp[nzb-1][ix][ih]+tmp[nzb-3][ix][ih]) + 1.0*d3*(tmp[nzb-4][ix][ih]);
                tmpq[nzb-1][ix][ih] = a*tmp[nzb-1][ix][ih] + 1.0*b3*(tmp[nzb-2][ix][ih]) + 1.0*d3*(tmp[nzb-3][ix][ih]);
                
                tmpq[0][ix][ih] = a[0][ix][ih]*tmp[0][ix][ih] + 1.0*b3[0][ix][ih]*(tmp[1][ix][ih]) + 1.0*d3[0][ix][ih]*(tmp[2][ix][ih]);
                tmpq[1][ix][ih] = a[1][ix][ih]*tmp[1][ix][ih] + b3[1][ix][ih]*(tmp[0][ix][ih]+tmp[2][ix][ih]) + 1.0*d3[1][ix][ih]*(tmp[3][ix][ih]);
                tmpq[nzb-2][ix][ih] = a[nzb-2][ix][ih]*tmp[nzb-2][ix][ih] + b3[nzb-2][ix][ih]*(tmp[nzb-1][ix][ih]+tmp[nzb-3][ix][ih]) + 1.0*d3[nzb-2][ix][ih]*(tmp[nzb-4][ix][ih]);
                tmpq[nzb-1][ix][ih] = a[nzb-1][ix][ih]*tmp[nzb-1][ix][ih] + 1.0*b3[nzb-1][ix][ih]*(tmp[nzb-2][ix][ih]) + 1.0*d3[nzb-1][ix][ih]*(tmp[nzb-3][ix][ih]);
               */ 
               tmpq[0][ix][ih] = a[0][ix][ih]*tmp[0][ix][ih] + 2.0*b3[0][ix][ih]*(tmp[1][ix][ih]) + 2.0*d3[0][ix][ih]*(tmp[2][ix][ih]);
                tmpq[1][ix][ih] = a[1][ix][ih]*tmp[1][ix][ih] + b3[1][ix][ih]*(tmp[0][ix][ih]+tmp[2][ix][ih]) + 2.0*d3[1][ix][ih]*(tmp[3][ix][ih]);
                tmpq[nzb-2][ix][ih] = a[nzb-2][ix][ih]*tmp[nzb-2][ix][ih] + b3[nzb-2][ix][ih]*(tmp[nzb-1][ix][ih]+tmp[nzb-3][ix][ih]) + 2.0*d3[nzb-2][ix][ih]*(tmp[nzb-4][ix][ih]);
                tmpq[nzb-1][ix][ih] = a[nzb-1][ix][ih]*tmp[nzb-1][ix][ih] + 2.0*b3[nzb-1][ix][ih]*(tmp[nzb-2][ix][ih]) + 2.0*d3[nzb-1][ix][ih]*(tmp[nzb-3][ix][ih]);
            }
        }
       // cosftz(false,nh,nx,nz,tmpq);
	for (iz=0; iz < nzb; iz++) {
	    for (ix=0; ix < nx; ix++) {
		for (ih=0; ih < nh; ih++) {
                    tmpq[iz][ix][ih]  += (2.0*curr[iz][ix][ih]-prev[iz][ix][ih]);
                    prev[iz][ix][ih]  = curr[iz][ix][ih];
                    curr[iz][ix][ih]  = tmpq[iz][ix][ih];
                }
            }
        }
        //cosftz(true,nh,nx,nz,curr);
        bdz_decay(prev);
        bdz_decay(curr);
	if (!mig) { /* modeling -> write out data */
	    cosft2(true,nh,nx,dat);
	    sf_floatwrite(dat[0],nh*nx,data);
	}
	sf_warning("    tmp=          %f           ;",curr[5][50][50]);
    }
    sf_warning(".");

    if (mig) {
	//for (iz=1; iz < nz; iz++) {
	for (iz=0; iz < nzb; iz++) {
	    for (ix=0; ix < nx; ix++) {
		for (ih=0; ih < nh; ih++) {
		    c = curr[iz][ix][ih];
		    img[ix][iz] += (ih==nh-1 || ih==0)? c*0.5: c;
		}
	    }
	}

	cosft2(true,nzb,nx,img);
	//sf_floatwrite(img[0],nz*nx,image);
	for (ix=0; ix < nx; ix++) {
            sf_floatwrite(img[ix]+nbt,nz,image);
	}
    }
    bd_close();
    exit(0);
}
