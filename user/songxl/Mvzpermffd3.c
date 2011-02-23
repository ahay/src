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
#include "abcpass3.h"


int main(int argc, char* argv[])
{
    bool mig;
    int it, nt, ix, nx, iz, nz, ih, nh, it1, it2, its, snap;
    int nzb, nbt, nbb, nxb, nxl, nxr, nhb, nhr, nhl;
    float ct, cb, chl, chr, cxl, cxr;
    float dt, dx, dz, dh, kx, kz, kh,  h, x, c, dkx, dkz, dkh, v0, pi=SF_PI;
    float ***prev, ***curr, **img, **dat, ***tmp, ***tmpq, *v;
    float *a, *b1, *b2, *b3, *c1, *c2, *c3, *d3, epsilon;
    float x1, x2, x3, v02, vz2, t2;
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
	/* depth samples (if migration) */
	if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
	/* depth sampling (if migration) */

	sf_putint(image,"n1",nz);
	sf_putfloat(image,"d1",dz);
	sf_putstring(image,"label1","Depth");
	sf_putint(image,"n2",nx);
	sf_putfloat(image,"d2",dx);
	sf_putstring(image,"label2","Midpoint");

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
	sf_putstring(data,"unit1","m");

	sf_putint(data,"n2",nx);
	sf_putfloat(data,"d2",dx);
	sf_putstring(data,"label2","Midpoint");
	sf_putstring(data,"unit2","m");

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
    if (!sf_getint("nxl",&nxl)) nxl=44;
    if (!sf_getint("nxr",&nxr)) nxr=44;
    if (!sf_getint("nhl",&nhl)) nhl=44;
    if (!sf_getint("nhr",&nhr)) nhr=44;

    if (!sf_getfloat("ct",&ct)) ct = 0.001; /*decaying parameter*/
    if (!sf_getfloat("cb",&cb)) cb = 0.001; /*decaying parameter*/
    if (!sf_getfloat("chl",&chl)) chl = 0.001; /*decaying parameter*/
    if (!sf_getfloat("chr",&chr)) chr = 0.001; /*decaying parameter*/
    if (!sf_getfloat("cxl",&cxl)) cxl = 0.001; /*decaying parameter*/
    if (!sf_getfloat("cxr",&cxr)) cxr = 0.001; /*decaying parameter*/

    nzb = nz + nbt + nbb;
    nxb = nx + nxl + nxr;
    nhb = nh + nhl + nhr;
  
    bd3_init(nz,nx,nh,nbt,nbb,nxl,nxr,nhl,nhr,ct,cb,cxl,cxr,chl,chr);

    if(mig){
      img = sf_floatalloc2(nz,nx);
    }
    else{
      img = sf_floatalloc2(nzb,nxb);
    }


    prev = sf_floatalloc3(nhb,nxb,nzb);
    curr = sf_floatalloc3(nhb,nxb,nzb);
    tmp  = sf_floatalloc3(nhb,nxb,nzb);
    tmpq  = sf_floatalloc3(nhb,nxb,nzb);

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

    dkx = cosft_dk(nxb,dx);
    dkz = cosft_dk(nzb,dz);
    dkh = cosft_dk(nhb,dh);

    a  = sf_floatalloc(nzb);
    b1 = sf_floatalloc(nzb);
    b2 = sf_floatalloc(nzb);
    b3 = sf_floatalloc(nzb);
    c1 = sf_floatalloc(nzb);
    c2 = sf_floatalloc(nzb);
    c3 = sf_floatalloc(nzb);
    d3 = sf_floatalloc(nzb);
    x1 = dh*dh;
    x2 = dx*dx;
    x3 = dz*dz;
    v02 = v0*v0;
    t2 = dt*dt;
    for (iz=0; iz < nzb; iz++) {
        vz = v[iz];
        vz2= vz*vz;
        d3[iz] = vz2*(t2*v02-t2*vz2+4.0*x3)/(48.0*v02*x3*x3);
        c1[iz] = vz2*(t2*v02-t2*vz2+2.0*x3)/(24.0*v02*x2*x3)-x3/x2*d3[iz];
        c2[iz] = vz2*(t2*v02-t2*vz2+2.0*x3)/(24.0*v02*x1*x3)-x3/x1*d3[iz];
        c3[iz] = vz2*t2*(v02-vz2)/(48.0*v02*x1*x2);
        b1[iz] = -2.0*(c2[iz]+c3[iz]);
        b2[iz] = -2.0*(c1[iz]+c3[iz]);
        b3[iz] = -2.0*(c2[iz]+c1[iz]+2.0*d3[iz])-vz2/(v02*x3);
        a[iz]  = -2.0*(b1[iz]+b2[iz]+b3[iz]+d3[iz])-4.0*(c1[iz]+c2[iz]+c3[iz]);
       
    }

   if (NULL != snaps) {
	sf_putint(snaps,"n1",nh);

	sf_putfloat(snaps,"d1",dh);
	sf_putstring(snaps,"label1","Half-Offset");

	sf_putint(snaps,"n2",nx);
	sf_putfloat(snaps,"d2",dx);
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
	for (ix=0; ix < nxb; ix++) {
	    for (ih=0; ih < nhb; ih++) {
		prev[iz][ix][ih] = 0.;
		curr[iz][ix][ih] = 0.;
	    }
	}
    }

    if (mig) { /* migration */
	/* initialize image */
	for (iz=0; iz < nz; iz++) {
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
	for (ix=nxl; ix < nx+nxl; ix++) {
            sf_floatread(img[ix]+nbt,nz,image);
	}
	for (ix=0; ix < nxl; ix++) {
	    for (iz=nbt; iz < nbt+nz; iz++) {
		img[ix][iz] = img[nxl][iz];
            }
        }
	for (ix=nxl+nx; ix < nxb; ix++) {
	    for (iz=nbt; iz < nbt+nz; iz++) {
		img[ix][iz] = img[nxl+nx-1][iz];
            }
        }
	for (ix=0; ix < nxb; ix++) {
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
	    for (ix=0; ix < nxb; ix++) {
		curr[iz][ix][nhl] = img[ix][iz];
	    }
	}
	//cosft3(false,nh,nx,nzb,curr);
	//cosfthx(false,nhb,nxb,nzb,curr);
	
	/* step forward in time */
	it1 = 0;
	it2 = nt;
	its = +1;
	cosft3(false,nhb,nxb,nzb,curr);
	for (iz=0; iz < nzb; iz++) {
	    kz = iz*dkz*2.0*pi;
	    for (ix=0; ix < nxb; ix++) {
		kx = ix*dkx*2.0*pi;
		for (ih=0; ih < nhb; ih++) {
		    kh = ih*dkh*2.0*pi;
                    if(kz < sqrtf(kx*kh)) { curr[iz][ix][ih] = 0.; }
                }
            }
        }
	cosft3(true,nhb,nxb,nzb,curr);
    }


    /* time stepping */
    for (it=it1; it != it2; it += its) {
	sf_warning("it=%d;",it);

	if (mig) { /* migration <- read data */
	    sf_floatread(dat[0],nh*nx,data);
	} else {
	    for (ix=0; ix < nx; ix++) {
		for (ih=0; ih < nh; ih++) {
		    dat[ix][ih] = 0.;
		}
	    }
	}

	if (NULL != snaps && 0 == it%snap) 
            for (iz=nbt; iz < nbt+nz; iz++) {
                for (ix=nxl; ix < nxl+nx; ix++) {
	            sf_floatwrite(curr[iz][ix]+nhl,nh,snaps);
                }
            }

	for (ix=0; ix < nx; ix++) {
            for (ih=0; ih < nh; ih++) {
		if (mig) {
                    //curr[nbt][ix+nxl][ih+nhl] += dat[ix][ih];
                    curr[nbt][ix+nxl][ih+nhl] = dat[ix][ih];
		    } else {
			dat[ix][ih] = curr[nbt][ix+nxl][ih+nhl];
		    }
            }
        }
        if (mig){
	   cosft3(false,nhb,nxb,nzb,curr);
	   for (iz=0; iz < nzb; iz++) {
	       kz = iz*dkz*2.0*pi;
	       for (ix=0; ix < nxb; ix++) {
	   	   kx = ix*dkx*2.0*pi;
		   for (ih=0; ih < nhb; ih++) {
		       kh = ih*dkh*2.0*pi;
                       if(kz < sqrtf(kx*kh) || kz<epsilon*0.001) { 
                         curr[iz][ix][ih] = 0.; 
                        }
                       tmp[iz][ix][ih] = curr[iz][ix][ih];
                    }
                }
           }
	cosft3(true,nhb,nxb,nzb,curr);
        }
        if(!mig){
        //if(1){
	  for (iz=0; iz < nzb; iz++) {
	      for (ix=0; ix < nxb; ix++) {
	 	  for (ih=0; ih < nhb; ih++) {
                      tmp[iz][ix][ih] = curr[iz][ix][ih];
                  }
              }
          }
           
	  cosft3(false,nhb,nxb,nzb,tmp);
        }
	for (iz=0; iz < nzb; iz++) {
	    kz = iz*dkz*2.0*pi;
	    for (ix=0; ix < nxb; ix++) {
		kx = ix*dkx*2.0*pi;
                   
	       if(!kz) 
		 x = ((kz+epsilon)*(kz+epsilon)+kx*kx)/(kz+epsilon);
               else 
		 x = (kz*kz+kx*kx)/(kz);

		for (ih=0; ih < nhb; ih++) {
		    kh = ih*dkh*2.0*pi;
                    if(kz < sqrtf(kx*kh) || kz < epsilon*0.001) { tmp[iz][ix][ih] = 0. ; continue;}
        //            if(kz < sqrtf(kx*kh)) { ih = nh; continue;}
		    //h = (kz*kz+kh*kh)/kz;
	            if(!kz) 
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
	              if(!kz) 
                      tmp[iz][ix][ih] = 2*(cosf(v0*dt*sqrtf(x*h)/2.0)-1.0)/((kz+epsilon)*(kz+epsilon))*c;
                      else 
		      tmp[iz][ix][ih] = 2*(cosf(v0*dt*sqrtf(x*h)/2.0)-1.0)/(kz*kz)*c;
		}
	    }
	}

	cosft3(true,nhb,nxb,nzb,tmp);
	sf_warning("      tmp=          %f           ;",tmp[102][145][90]);
	for (iz=0; iz < nzb; iz++) {
	    for (ix=0; ix < nxb; ix++) {
		for (ih=0; ih < nhb; ih++) {
                    tmpq[iz][ix][ih] = 0.0;
                }
            }
        }
	for (iz=2; iz < nzb-2; iz++) {
	    for (ix=1; ix < nxb-1; ix++) {
		for (ih=1; ih < nhb-1; ih++) {
                    tmpq[iz][ix][ih] = a[iz]*tmp[iz][ix][ih] + b1[iz]*(tmp[iz][ix][ih-1]+tmp[iz][ix][ih+1])+b2[iz]*(tmp[iz][ix-1][ih]+tmp[iz][ix+1][ih])+b3[iz]*(tmp[iz+1][ix][ih]+tmp[iz-1][ix][ih]) + d3[iz]*(tmp[iz+2][ix][ih]+tmp[iz-2][ix][ih])+c1[iz]*(tmp[iz-1][ix-1][ih]+tmp[iz-1][ix+1][ih]+tmp[iz+1][ix-1][ih]+tmp[iz+1][ix+1][ih])+c2[iz]*(tmp[iz-1][ix][ih-1]+tmp[iz-1][ix][ih+1]+tmp[iz+1][ix][ih-1]+tmp[iz+1][ix][ih+1])+c3[iz]*(tmp[iz][ix-1][ih-1]+tmp[iz][ix-1][ih+1]+tmp[iz][ix+1][ih-1]+tmp[iz][ix+1][ih+1]);
                }
                    tmpq[iz][ix][0] = a[iz]*tmp[iz][ix][0] + 2.0*b1[iz]*(tmp[iz][ix][0+1])+b2[iz]*(tmp[iz][ix-1][0]+tmp[iz][ix+1][0])+b3[iz]*(tmp[iz+1][ix][0]+tmp[iz-1][ix][0]) + d3[iz]*(tmp[iz+2][ix][0]+tmp[iz-2][ix][0])+c1[iz]*(tmp[iz-1][ix-1][0]+tmp[iz-1][ix+1][0]+tmp[iz+1][ix-1][0]+tmp[iz+1][ix+1][0])+2.0*c2[iz]*(tmp[iz-1][ix][1]+tmp[iz+1][ix][1])+2.0*c3[iz]*(tmp[iz][ix-1][1]+tmp[iz][ix+1][1]);
                    tmpq[iz][ix][nhb-1] = a[iz]*tmp[iz][ix][nhb-1] + 2.0*b1[iz]*(tmp[iz][ix][nhb-2])+b2[iz]*(tmp[iz][ix-1][nhb-1]+tmp[iz][ix+1][nhb-1])+b3[iz]*(tmp[iz+1][ix][nhb-1]+tmp[iz-1][ix][nhb-1]) + d3[iz]*(tmp[iz+2][ix][nhb-1]+tmp[iz-2][ix][nhb-1])+c1[iz]*(tmp[iz-1][ix-1][nhb-1]+tmp[iz-1][ix+1][nhb-1]+tmp[iz+1][ix-1][nhb-1]+tmp[iz+1][ix+1][nhb-1])+2.0*c2[iz]*(tmp[iz-1][ix][nhb-2]+tmp[iz+1][ix][nhb-2])+2.0*c3[iz]*(tmp[iz][ix-1][nhb-2]+tmp[iz][ix+1][nhb-2]);
            }
        }
	for (iz=2; iz < nzb-2; iz++) {
            for (ih=1; ih < nhb-1; ih++) {
                tmpq[iz][0][ih] = a[iz]*tmp[iz][0][ih] + b1[iz]*(tmp[iz][0][ih-1]+tmp[iz][0][ih+1])+2.0*b2[iz]*(tmp[iz][1][ih])+b3[iz]*(tmp[iz+1][0][ih]+tmp[iz-1][0][ih]) + d3[iz]*(tmp[iz+2][0][ih]+tmp[iz-2][0][ih])+2.0*c1[iz]*(tmp[iz-1][1][ih]+tmp[iz+1][1][ih])+c2[iz]*(tmp[iz-1][0][ih-1]+tmp[iz-1][0][ih+1]+tmp[iz+1][0][ih-1]+tmp[iz+1][0][ih+1])+2.0*c3[iz]*(tmp[iz][1][ih-1]+tmp[iz][1][ih+1]);
                tmpq[iz][nxb-1][ih] = a[iz]*tmp[iz][nxb-1][ih] + b1[iz]*(tmp[iz][nxb-1][ih-1]+tmp[iz][nxb-1][ih+1])+2.0*b2[iz]*(tmp[iz][nxb-2][ih])+b3[iz]*(tmp[iz+1][nxb-1][ih]+tmp[iz-1][nxb-1][ih]) + d3[iz]*(tmp[iz+2][nxb-1][ih]+tmp[iz-2][nxb-1][ih])+2.0*c1[iz]*(tmp[iz-1][nxb-2][ih]+tmp[iz+1][nxb-2][ih])+c2[iz]*(tmp[iz-1][nxb-1][ih-1]+tmp[iz-1][nxb-1][ih+1]+tmp[iz+1][nxb-1][ih-1]+tmp[iz+1][nxb-1][ih+1])+2.0*c3[iz]*(tmp[iz][nxb-2][ih-1]+tmp[iz][nxb-2][ih+1]);
            }
        }
	for (ix=1; ix < nxb-1; ix++) {
            for (ih=1; ih < nhb-1; ih++) {
                tmpq[0][ix][ih] = a[0]*tmp[0][ix][ih] + b1[0]*(tmp[0][ix][ih-1]+tmp[0][ix][ih+1])+b2[0]*(tmp[0][ix-1][ih]+tmp[0][ix+1][ih])+2.0*b3[0]*(tmp[1][ix][ih]) + 2.0*d3[0]*(tmp[2][ix][ih])+2.0*c1[0]*(tmp[1][ix-1][ih]+tmp[1][ix+1][ih])+2.0*c2[0]*(tmp[1][ix][ih-1]+tmp[1][ix][ih+1])+c3[0]*(tmp[0][ix-1][ih-1]+tmp[0][ix-1][ih+1]+tmp[0][ix+1][ih-1]+tmp[0][ix+1][ih+1]);
                tmpq[nzb-1][ix][ih] = a[nzb-1]*tmp[nzb-1][ix][ih] + b1[nzb-1]*(tmp[nzb-1][ix][ih-1]+tmp[nzb-1][ix][ih+1])+b2[nzb-1]*(tmp[nzb-1][ix-1][ih]+tmp[nzb-1][ix+1][ih])+2.0*b3[nzb-1]*(tmp[nzb-2][ix][ih]) + 2.0*d3[nzb-1]*(tmp[nzb-3][ix][ih])+2.0*c1[nzb-1]*(tmp[nzb-2][ix-1][ih]+tmp[nzb-2][ix+1][ih])+2.0*c2[nzb-1]*(tmp[nzb-2][ix][ih-1]+tmp[nzb-2][ix][ih+1])+c3[nzb-1]*(tmp[nzb-1][ix-1][ih-1]+tmp[nzb-1][ix-1][ih+1]+tmp[nzb-1][ix+1][ih-1]+tmp[nzb-1][ix+1][ih+1]);
                    tmpq[1][ix][ih] = a[1]*tmp[1][ix][ih] + b1[1]*(tmp[1][ix][ih-1]+tmp[1][ix][ih+1])+b2[1]*(tmp[1][ix-1][ih]+tmp[1][ix+1][ih])+b3[1]*(tmp[2][ix][ih]+tmp[0][ix][ih]) + 2.0*d3[1]*(tmp[3][ix][ih])+c1[1]*(tmp[0][ix-1][ih]+tmp[0][ix+1][ih]+tmp[2][ix-1][ih]+tmp[2][ix+1][ih])+c2[1]*(tmp[0][ix][ih-1]+tmp[0][ix][ih+1]+tmp[2][ix][ih-1]+tmp[2][ix][ih+1])+c3[1]*(tmp[1][ix-1][ih-1]+tmp[1][ix-1][ih+1]+tmp[1][ix+1][ih-1]+tmp[1][ix+1][ih+1]);
                    tmpq[nzb-2][ix][ih] = a[nzb-2]*tmp[nzb-2][ix][ih] + b1[nzb-2]*(tmp[nzb-2][ix][ih-1]+tmp[nzb-2][ix][ih+1])+b2[nzb-2]*(tmp[nzb-2][ix-1][ih]+tmp[nzb-2][ix+1][ih])+b3[nzb-2]*(tmp[nzb-1][ix][ih]+tmp[nzb-3][ix][ih]) + 2.0*d3[nzb-2]*(tmp[nzb-4][ix][ih])+c1[nzb-2]*(tmp[nzb-3][ix-1][ih]+tmp[nzb-3][ix+1][ih]+tmp[nzb-1][ix-1][ih]+tmp[nzb-1][ix+1][ih])+c2[nzb-2]*(tmp[nzb-3][ix][ih-1]+tmp[nzb-3][ix][ih+1]+tmp[nzb-1][ix][ih-1]+tmp[nzb-1][ix][ih+1])+c3[nzb-2]*(tmp[nzb-2][ix-1][ih-1]+tmp[nzb-2][ix-1][ih+1]+tmp[nzb-2][ix+1][ih-1]+tmp[nzb-2][ix+1][ih+1]);
            }
        }
	for (iz=2; iz < nzb-2; iz++) {
            tmpq[iz][0][0] = a[iz]*tmp[iz][0][0] + 2.0*b1[iz]*(tmp[iz][0][1])+2.0*b2[iz]*(tmp[iz][1][0])+b3[iz]*(tmp[iz+1][0][0]+tmp[iz-1][0][0]) + d3[iz]*(tmp[iz+2][0][0]+tmp[iz-2][0][0])+2.0*c1[iz]*(tmp[iz-1][1][0]+tmp[iz+1][1][0])+2.0*c2[iz]*(tmp[iz-1][0][1]+tmp[iz+1][0][1])+4.0*c3[iz]*(tmp[iz][1][1]);
            tmpq[iz][nxb-1][nhb-1] = a[iz]*tmp[iz][nxb-1][nhb-1] + 2.0*b1[iz]*(tmp[iz][nxb-1][nhb-2])+2.0*b2[iz]*(tmp[iz][nxb-2][nhb-1])+b3[iz]*(tmp[iz+1][nxb-1][nhb-1]+tmp[iz-1][nxb-1][nhb-1]) + d3[iz]*(tmp[iz+2][nxb-1][nhb-1]+tmp[iz-2][nxb-1][nhb-1])+2.0*c1[iz]*(tmp[iz-1][nxb-2][nhb-1]+tmp[iz+1][nxb-2][nhb-1])+2.0*c2[iz]*(tmp[iz-1][nxb-1][nhb-2]+tmp[iz+1][nxb-1][nhb-2])+4.0*c3[iz]*(tmp[iz][nxb-2][nhb-2]);
            tmpq[iz][0][nhb-1] = a[iz]*tmp[iz][0][nhb-1] + 2.0*b1[iz]*(tmp[iz][0][nhb-2])+2.0*b2[iz]*(tmp[iz][1][nhb-1])+b3[iz]*(tmp[iz+1][0][nhb-1]+tmp[iz-1][0][nhb-1]) + d3[iz]*(tmp[iz+2][0][nhb-1]+tmp[iz-2][0][nhb-1])+2.0*c1[iz]*(tmp[iz-1][1][nhb-1]+tmp[iz+1][1][nhb-1])+2.0*c2[iz]*(tmp[iz-1][0][nhb-2]+tmp[iz+1][0][nhb-2])+4.0*c3[iz]*(tmp[iz][1][nhb-2]);
            tmpq[iz][nxb-1][0] = a[iz]*tmp[iz][nxb-1][0] + 2.0*b1[iz]*(tmp[iz][nxb-1][1])+2.0*b2[iz]*(tmp[iz][nxb-2][0])+b3[iz]*(tmp[iz+1][nxb-1][0]+tmp[iz-1][nxb-1][0]) + d3[iz]*(tmp[iz+2][nxb-1][0]+tmp[iz-2][nxb-1][0])+2.0*c1[iz]*(tmp[iz-1][nxb-2][0]+tmp[iz+1][nxb-2][0])+2.0*c2[iz]*(tmp[iz-1][nxb-1][1]+tmp[iz+1][nxb-1][1])+4.0*c3[iz]*(tmp[iz][nxb-2][1]);
        }
	for (ix=1; ix < nxb-1; ix++) {
            tmpq[0][ix][0] = a[0]*tmp[0][ix][0] + 2.0*b1[0]*(tmp[0][ix][1])+b2[0]*(tmp[0][ix-1][0]+tmp[0][ix+1][0])+2.0*b3[0]*(tmp[1][ix][0]) + 2.0*d3[0]*(tmp[2][ix][0])+2.0*c1[0]*(tmp[1][ix-1][0]+tmp[1][ix+1][0])+4.0*c2[0]*(tmp[1][ix][1])+2.0*c3[0]*(tmp[0][ix-1][1]+tmp[0][ix+1][1]);
            tmpq[0][ix][nhb-1] = a[0]*tmp[0][ix][nhb-1] + 2.0*b1[0]*(tmp[0][ix][nhb-2])+b2[0]*(tmp[0][ix-1][nhb-1]+tmp[0][ix+1][nhb-1])+2.0*b3[0]*(tmp[1][ix][nhb-1]) + 2.0*d3[0]*(tmp[2][ix][nhb-1])+2.0*c1[0]*(tmp[1][ix-1][nhb-1]+tmp[1][ix+1][nhb-1])+4.0*c2[0]*(tmp[1][ix][nhb-2])+2.0*c3[0]*(tmp[0][ix-1][nhb-2]+tmp[0][ix+1][nhb-2]);
            tmpq[nzb-1][ix][0] = a[nzb-1]*tmp[nzb-1][ix][0] + 2.0*b1[nzb-1]*(tmp[nzb-1][ix][1])+b2[nzb-1]*(tmp[nzb-1][ix-1][0]+tmp[nzb-1][ix+1][0])+2.0*b3[nzb-1]*(tmp[nzb-2][ix][0]) + 2.0*d3[nzb-1]*(tmp[nzb-3][ix][0])+2.0*c1[nzb-1]*(tmp[nzb-2][ix-1][0]+tmp[nzb-2][ix+1][0])+4.0*c2[nzb-1]*(tmp[nzb-2][ix][1])+2.0*c3[nzb-1]*(tmp[nzb-1][ix-1][1]+tmp[nzb-1][ix+1][1]);
            tmpq[nzb-1][ix][nhb-1] = a[nzb-1]*tmp[nzb-1][ix][nhb-1] + 2.0*b1[nzb-1]*(tmp[nzb-1][ix][nhb-2])+b2[nzb-1]*(tmp[nzb-1][ix-1][nhb-1]+tmp[nzb-1][ix+1][nhb-1])+2.0*b3[nzb-1]*(tmp[nzb-2][ix][nhb-1]) + 2.0*d3[nzb-1]*(tmp[nzb-3][ix][nhb-1])+2.0*c1[nzb-1]*(tmp[nzb-2][ix-1][nhb-1]+tmp[nzb-2][ix+1][nhb-1])+4.0*c2[nzb-1]*(tmp[nzb-2][ix][nhb-2])+2.0*c3[nzb-1]*(tmp[nzb-1][ix-1][nhb-2]+tmp[nzb-1][ix+1][nhb-2]);
            tmpq[1][ix][0] = a[1]*tmp[1][ix][0] + 2.0*b1[1]*(tmp[1][ix][1])+b2[1]*(tmp[1][ix-1][0]+tmp[1][ix+1][0])+b3[1]*(tmp[2][ix][0]+tmp[0][ix][0]) + 2.0*d3[1]*(tmp[3][ix][0])+c1[1]*(tmp[0][ix-1][0]+tmp[0][ix+1][0]+tmp[2][ix-1][0]+tmp[2][ix+1][0])+2.0*c2[1]*(tmp[0][ix][1]+tmp[2][ix][1])+2.0*c3[1]*(tmp[1][ix-1][1]+tmp[1][ix+1][1]);
            tmpq[1][ix][nhb-1] = a[1]*tmp[1][ix][nhb-1] + 2.0*b1[1]*(tmp[1][ix][nhb-2])+b2[1]*(tmp[1][ix-1][nhb-1]+tmp[1][ix+1][nhb-1])+b3[1]*(tmp[2][ix][nhb-1]+tmp[0][ix][nhb-1]) + 2.0*d3[1]*(tmp[3][ix][nhb-1])+c1[1]*(tmp[0][ix-1][nhb-1]+tmp[0][ix+1][nhb-1]+tmp[2][ix-1][nhb-1]+tmp[2][ix+1][nhb-1])+2.0*c2[1]*(tmp[0][ix][nhb-2]+tmp[2][ix][nhb-2])+2.0*c3[1]*(tmp[1][ix-1][nhb-2]+tmp[1][ix+1][nhb-2]);
            tmpq[nzb-2][ix][0] = a[nzb-2]*tmp[nzb-2][ix][0] + 2.0*b1[nzb-2]*(tmp[nzb-2][ix][1])+b2[nzb-2]*(tmp[nzb-2][ix-1][0]+tmp[nzb-2][ix+1][0])+b3[nzb-2]*(tmp[nzb-1][ix][0]+tmp[nzb-3][ix][0]) + 2.0*d3[nzb-2]*(tmp[nzb-4][ix][0])+c1[nzb-2]*(tmp[nzb-3][ix-1][0]+tmp[nzb-3][ix+1][0]+tmp[nzb-1][ix-1][0]+tmp[nzb-1][ix+1][0])+2.0*c2[nzb-2]*(tmp[nzb-3][ix][1]+tmp[nzb-1][ix][1])+2.0*c3[nzb-2]*(tmp[nzb-2][ix-1][1]+tmp[nzb-2][ix+1][1]);
            tmpq[nzb-2][ix][nhb-1] = a[nzb-2]*tmp[nzb-2][ix][nhb-1] + 2.0*b1[nzb-2]*(tmp[nzb-2][ix][nhb-2])+b2[nzb-2]*(tmp[nzb-2][ix-1][nhb-1]+tmp[nzb-2][ix+1][nhb-1])+b3[nzb-2]*(tmp[nzb-1][ix][nhb-1]+tmp[nzb-3][ix][nhb-1]) + 2.0*d3[nzb-2]*(tmp[nzb-4][ix][nhb-1])+c1[nzb-2]*(tmp[nzb-3][ix-1][nhb-1]+tmp[nzb-3][ix+1][nhb-1]+tmp[nzb-1][ix-1][nhb-1]+tmp[nzb-1][ix+1][nhb-1])+2.0*c2[nzb-2]*(tmp[nzb-3][ix][nhb-2]+tmp[nzb-1][ix][nhb-2])+2.0*c3[nzb-2]*(tmp[nzb-2][ix-1][nhb-2]+tmp[nzb-2][ix+1][nhb-2]);
        }
        for (ih=1; ih < nhb-1; ih++) {
            tmpq[0][0][ih] = a[0]*tmp[0][0][ih] + b1[0]*(tmp[0][0][ih-1]+tmp[0][0][ih+1])+2.0*b2[0]*(tmp[0][1][ih])+2.0*b3[0]*(tmp[1][0][ih]) + 2.0*d3[0]*(tmp[2][0][ih])+4.0*c1[0]*(tmp[1][1][ih])+2.0*c2[0]*(tmp[1][0][ih-1]+tmp[1][0][ih+1])+2.0*c3[0]*(tmp[0][1][ih-1]+tmp[0][1][ih+1]);
            tmpq[0][nxb-1][ih] = a[0]*tmp[0][nxb-1][ih] + b1[0]*(tmp[0][nxb-1][ih-1]+tmp[0][nxb-1][ih+1])+2.0*b2[0]*(tmp[0][nxb-2][ih])+2.0*b3[0]*(tmp[1][nxb-1][ih]) + 2.0*d3[0]*(tmp[2][nxb-1][ih])+4.0*c1[0]*(tmp[1][nxb-2][ih])+2.0*c2[0]*(tmp[1][nxb-1][ih-1]+tmp[1][nxb-1][ih+1])+2.0*c3[0]*(tmp[0][nxb-2][ih-1]+tmp[0][nxb-2][ih+1]);
            tmpq[nzb-1][0][ih] = a[nzb-1]*tmp[nzb-1][0][ih] + b1[nzb-1]*(tmp[nzb-1][0][ih-1]+tmp[nzb-1][0][ih+1])+2.0*b2[nzb-1]*(tmp[nzb-1][0+1][ih])+2.0*b3[nzb-1]*(tmp[nzb-2][0][ih]) + 2.0*d3[nzb-1]*(tmp[nzb-3][0][ih])+4.0*c1[nzb-1]*(tmp[nzb-2][1][ih])+2.0*c2[nzb-1]*(tmp[nzb-2][0][ih-1]+tmp[nzb-2][0][ih+1])+2.0*c3[nzb-1]*(tmp[nzb-1][1][ih-1]+tmp[nzb-1][1][ih+1]);
            tmpq[nzb-1][nxb-1][ih] = a[nzb-1]*tmp[nzb-1][nxb-1][ih] + b1[nzb-1]*(tmp[nzb-1][nxb-1][ih-1]+tmp[nzb-1][nxb-1][ih+1])+2.0*b2[nzb-1]*(tmp[nzb-1][nxb-2][ih])+2.0*b3[nzb-1]*(tmp[nzb-2][nxb-1][ih]) + 2.0*d3[nzb-1]*(tmp[nzb-3][nxb-1][ih])+4.0*c1[nzb-1]*(tmp[nzb-2][nxb-2][ih])+2.0*c2[nzb-1]*(tmp[nzb-2][nxb-1][ih-1]+tmp[nzb-2][nxb-1][ih+1])+2.0*c3[nzb-1]*(tmp[nzb-1][nxb-2][ih-1]+tmp[nzb-1][nxb-2][ih+1]);
            tmpq[1][0][ih] = a[1]*tmp[1][0][ih] + b1[1]*(tmp[1][0][ih-1]+tmp[1][0][ih+1])+2.0*b2[1]*(tmp[1][0+1][ih])+b3[1]*(tmp[2][0][ih]+tmp[0][0][ih]) + 2.0*d3[1]*(tmp[3][0][ih])+2.0*c1[1]*(tmp[0][1][ih]+tmp[2][1][ih])+c2[1]*(tmp[0][0][ih-1]+tmp[0][0][ih+1]+tmp[2][0][ih-1]+tmp[2][0][ih+1])+2.0*c3[1]*(tmp[1][1][ih-1]+tmp[1][1][ih+1]);
            tmpq[1][nxb-1][ih] = a[1]*tmp[1][nxb-1][ih] + b1[1]*(tmp[1][nxb-1][ih-1]+tmp[1][nxb-1][ih+1])+2.0*b2[1]*(tmp[1][nxb-2][ih])+b3[1]*(tmp[2][nxb-1][ih]+tmp[0][nxb-1][ih]) + 2.0*d3[1]*(tmp[3][nxb-1][ih])+2.0*c1[1]*(tmp[0][nxb-2][ih]+tmp[2][nxb-2][ih])+c2[1]*(tmp[0][nxb-1][ih-1]+tmp[0][nxb-1][ih+1]+tmp[2][nxb-1][ih-1]+tmp[2][nxb-1][ih+1])+2.0*c3[1]*(tmp[1][nxb-2][ih-1]+tmp[1][nxb-2][ih+1]);
            tmpq[nzb-2][0][ih] = a[nzb-2]*tmp[nzb-2][0][ih] + b1[nzb-2]*(tmp[nzb-2][0][ih-1]+tmp[nzb-2][0][ih+1])+2.0*b2[nzb-2]*(tmp[nzb-2][0+1][ih])+b3[nzb-2]*(tmp[nzb-1][0][ih]+tmp[nzb-3][0][ih]) + 2.0*d3[nzb-2]*(tmp[nzb-4][0][ih])+2.0*c1[nzb-2]*(tmp[nzb-3][1][ih]+tmp[nzb-1][1][ih])+c2[nzb-2]*(tmp[nzb-3][0][ih-1]+tmp[nzb-3][0][ih+1]+tmp[nzb-1][0][ih-1]+tmp[nzb-1][0][ih+1])+2.0*c3[nzb-2]*(tmp[nzb-2][1][ih-1]+tmp[nzb-2][1][ih+1]);
            tmpq[nzb-2][nxb-1][ih] = a[nzb-2]*tmp[nzb-2][nxb-1][ih] + b1[nzb-2]*(tmp[nzb-2][nxb-1][ih-1]+tmp[nzb-2][nxb-1][ih+1])+2.0*b2[nzb-2]*(tmp[nzb-2][nxb-2][ih])+b3[nzb-2]*(tmp[nzb-1][nxb-1][ih]+tmp[nzb-3][nxb-1][ih]) + 2.0*d3[nzb-2]*(tmp[nzb-4][nxb-1][ih])+2.0*c1[nzb-2]*(tmp[nzb-3][nxb-2][ih]+tmp[nzb-1][nxb-2][ih])+c2[nzb-2]*(tmp[nzb-3][nxb-1][ih-1]+tmp[nzb-3][nxb-1][ih+1]+tmp[nzb-1][nxb-1][ih-1]+tmp[nzb-1][nxb-1][ih+1])+2.0*c3[nzb-2]*(tmp[nzb-2][nxb-2][ih-1]+tmp[nzb-2][nxb-2][ih+1]);
        }
            tmpq[0][0][0] = a[0]*tmp[0][0][0] + 2.0*b1[0]*(tmp[0][0][1])+2.0*b2[0]*(tmp[0][1][0])+2.0*b3[0]*(tmp[1][0][0]) + 2.0*d3[0]*(tmp[2][0][0])+4.0*c1[0]*(tmp[1][1][0])+4.0*c2[0]*(tmp[1][0][1])+4.0*c3[0]*(tmp[0][1][1]);
            tmpq[0][0][nhb-1] = a[0]*tmp[0][0][nhb-1] + 2.0*b1[0]*(tmp[0][0][nhb-2])+2.0*b2[0]*(tmp[0][1][nhb-1])+2.0*b3[0]*(tmp[1][0][nhb-1]) + 2.0*d3[0]*(tmp[2][0][nhb-1])+4.0*c1[0]*(tmp[1][1][nhb-1])+4.0*c2[0]*(tmp[1][0][nhb-2])+4.0*c3[0]*(tmp[0][1][nhb-2]);
            tmpq[0][nxb-1][0] = a[0]*tmp[0][nxb-1][0] + 2.0*b1[0]*(tmp[0][nxb-1][1])+2.0*b2[0]*(tmp[0][nxb-2][0])+2.0*b3[0]*(tmp[1][nxb-1][0]) + 2.0*d3[0]*(tmp[2][nxb-1][0])+4.0*c1[0]*(tmp[1][nxb-2][0])+4.0*c2[0]*(tmp[1][nxb-1][1])+4.0*c3[0]*(tmp[0][nxb-2][1]);
            tmpq[0][nxb-1][nhb-1] = a[0]*tmp[0][nxb-1][nhb-1] + 2.0*b1[0]*(tmp[0][nxb-1][nhb-2])+2.0*b2[0]*(tmp[0][nxb-2][nhb-1])+2.0*b3[0]*(tmp[1][nxb-1][nhb-1]) + 2.0*d3[0]*(tmp[2][nxb-1][nhb-1])+4.0*c1[0]*(tmp[1][nxb-2][nhb-1])+4.0*c2[0]*(tmp[1][nxb-1][nhb-2])+4.0*c3[0]*(tmp[0][nxb-2][nhb-2]);
            tmpq[nzb-1][0][0] = a[nzb-1]*tmp[nzb-1][0][0] + 2.0*b1[nzb-1]*(tmp[nzb-1][0][1])+2.0*b2[nzb-1]*(tmp[nzb-1][1][0])+2.0*b3[nzb-1]*(tmp[nzb-2][0][0]) +2.0* d3[nzb-1]*(tmp[nzb-3][0][0])+4.0*c1[nzb-1]*(tmp[nzb-2][1][0])+4.0*c2[nzb-1]*(tmp[nzb-2][0][1])+2.0*c3[nzb-1]*(tmp[nzb-1][1][1]);
            tmpq[nzb-1][0][nhb-1] = a[nzb-1]*tmp[nzb-1][0][nhb-1] + 2.0*b1[nzb-1]*(tmp[nzb-1][0][nhb-2])+2.0*b2[nzb-1]*(tmp[nzb-1][1][nhb-1])+2.0*b3[nzb-1]*(tmp[nzb-2][0][nhb-1]) + 2.0*d3[nzb-1]*(tmp[nzb-3][0][nhb-1])+4.0*c1[nzb-1]*(tmp[nzb-2][1][nhb-1])+4.0*c2[nzb-1]*(tmp[nzb-2][0][nhb-2])+4.0*c3[nzb-1]*(tmp[nzb-1][1][nhb-2]);
            tmpq[nzb-1][nxb-1][0] = a[nzb-1]*tmp[nzb-1][nxb-1][0] + 2.0*b1[nzb-1]*(tmp[nzb-1][nxb-1][1])+2.0*b2[nzb-1]*(tmp[nzb-1][nxb-2][0])+2.0*b3[nzb-1]*(tmp[nzb-2][nxb-1][0]) + 2.0*d3[nzb-1]*(tmp[nzb-3][nxb-1][0])+4.0*c1[nzb-1]*(tmp[nzb-2][nxb-2][0])+4.0*c2[nzb-1]*(tmp[nzb-2][nxb-1][1])+4.0*c3[nzb-1]*(tmp[nzb-1][nxb-2][1]);
            tmpq[nzb-1][nxb-1][nhb-1] = a[nzb-1]*tmp[nzb-1][nxb-1][nhb-1] +2.0*b1[nzb-1]*(tmp[nzb-1][nxb-1][nhb-2])+2.0*b2[nzb-1]*(tmp[nzb-1][nxb-2][nhb-1])+2.0*b3[nzb-1]*(tmp[nzb-2][nxb-1][nhb-1]) + 2.0*d3[nzb-1]*(tmp[nzb-3][nxb-1][nhb-1])+4.0*c1[nzb-1]*(tmp[nzb-2][nxb-2][nhb-1])+4.0*c2[nzb-1]*(tmp[nzb-2][nxb-1][nhb-2])+4.0*c3[nzb-1]*(tmp[nzb-1][nxb-2][nhb-2]);
            tmpq[1][0][0] = a[1]*tmp[1][0][0] + 2.0*b1[1]*(tmp[1][0][0+1])+2.0*b2[1]*(tmp[1][0+1][0])+b3[1]*(tmp[2][0][0]+tmp[0][0][0]) + 2.0*d3[1]*(tmp[3][0][0])+2.0*c1[1]*(tmp[0][1][0]+tmp[2][0+1][0])+2.0*c2[1]*(tmp[0][0][1]+tmp[2][0][1])+4.0*c3[1]*(tmp[1][1][1]);
            tmpq[1][0][nhb-1] = a[1]*tmp[1][0][nhb-1] + 2.0*b1[1]*(tmp[1][0][nhb-2])+2.0*b2[1]*(tmp[1][1][nhb-1])+b3[1]*(tmp[2][0][nhb-1]+tmp[0][0][nhb-1]) + 2.0*d3[1]*(tmp[3][0][nhb-1])+2.0*c1[1]*(tmp[0][1][nhb-1]+tmp[2][1][nhb-1])+2.0*c2[1]*(tmp[0][0][nhb-2]+tmp[2][0][nhb-2])+4.0*c3[1]*(tmp[1][1][nhb-2]);
            tmpq[1][nxb-1][0] = a[1]*tmp[1][nxb-1][0] + 2.0*b1[1]*(tmp[1][nxb-1][1])+2.0*b2[1]*(tmp[1][nxb-2][0])+b3[1]*(tmp[2][nxb-1][0]+tmp[0][nxb-1][0]) + 2.0*d3[1]*(tmp[3][nxb-1][0])+2.0*c1[1]*(tmp[0][nxb-2][0]+tmp[2][nxb-2][0])+2.0*c2[1]*(tmp[0][nxb-1][1]+tmp[2][nxb-1][1])+4.0*c3[1]*(tmp[1][nxb-2][1]);
            tmpq[1][nxb-1][nhb-1] = a[1]*tmp[1][nxb-1][nhb-1] + 2.0*b1[1]*(tmp[1][nxb-1][nhb-2])+2.0*b2[1]*(tmp[1][nxb-2][nhb-1])+b3[1]*(tmp[2][nxb-1][nhb-1]+tmp[0][nxb-1][nhb-1]) + 2.0*d3[1]*(tmp[3][nxb-1][nhb-1])+2.0*c1[1]*(tmp[0][nxb-2][nhb-1]+tmp[2][nxb-2][nhb-1])+2.0*c2[1]*(tmp[0][nxb-1][nhb-2]+tmp[2][nxb-1][nhb-2])+4.0*c3[1]*(tmp[1][nxb-2][nhb-2]);
            tmpq[nzb-2][0][0] = a[nzb-2]*tmp[nzb-2][0][0] + 2.0*b1[nzb-2]*(tmp[nzb-2][0][1])+2.0*b2[nzb-2]*(tmp[nzb-2][1][0])+b3[nzb-2]*(tmp[nzb-1][0][0]+tmp[nzb-3][0][0]) + 2.0*d3[nzb-2]*(tmp[nzb-4][0][0])+2.0*c1[nzb-2]*(tmp[nzb-3][1][0]+tmp[nzb-1][1][0])+2.0*c2[nzb-2]*(tmp[nzb-3][0][1]+tmp[nzb-1][0][1])+4.0*c3[nzb-2]*(tmp[nzb-2][1][1]);
            tmpq[nzb-2][0][nhb-1] = a[nzb-2]*tmp[nzb-2][0][nhb-1] + 2.0*b1[nzb-2]*(tmp[nzb-2][0][nhb-2])+2.0*b2[nzb-2]*(tmp[nzb-2][1][nhb-1])+b3[nzb-2]*(tmp[nzb-1][0][nhb-1]+tmp[nzb-3][0][nhb-1]) + 2.0*d3[nzb-2]*(tmp[nzb-4][0][nhb-1])+2.0*c1[nzb-2]*(tmp[nzb-3][1][nhb-1]+tmp[nzb-1][1][nhb-1])+2.0*c2[nzb-2]*(tmp[nzb-3][0][nhb-2]+tmp[nzb-1][0][nhb-2])+4.0*c3[nzb-2]*(tmp[nzb-2][1][nhb-2]);
            tmpq[nzb-2][nxb-1][0] = a[nzb-2]*tmp[nzb-2][nxb-1][0] + 2.0*b1[nzb-2]*(tmp[nzb-2][nxb-1][1])+2.0*b2[nzb-2]*(tmp[nzb-2][nxb-2][0])+b3[nzb-2]*(tmp[nzb-1][nxb-1][0]+tmp[nzb-3][nxb-1][0]) + 2.0*d3[nzb-2]*(tmp[nzb-4][nxb-1][0])+2.0*c1[nzb-2]*(tmp[nzb-3][nxb-2][0]+tmp[nzb-1][nxb-2][0])+2.0*c2[nzb-2]*(tmp[nzb-3][nxb-1][1]+tmp[nzb-1][nxb-1][1])+4.0*c3[nzb-2]*(tmp[nzb-2][nxb-2][1]);
            tmpq[nzb-2][nxb-1][nhb-1] = a[nzb-2]*tmp[nzb-2][nxb-1][nhb-1] + 2.0*b1[nzb-2]*(tmp[nzb-2][nxb-1][nhb-2])+2.0*b2[nzb-2]*(tmp[nzb-2][nxb-2][nhb-1])+b3[nzb-2]*(tmp[nzb-1][nxb-1][nhb-1]+tmp[nzb-3][nxb-1][nhb-1]) + 2.0*d3[nzb-2]*(tmp[nzb-4][nxb-1][nhb-1])+2.0*c1[nzb-2]*(tmp[nzb-3][nxb-2][nhb-1]+tmp[nzb-1][nxb-2][nhb-1])+2.0*c2[nzb-2]*(tmp[nzb-3][nxb-1][nhb-2]+tmp[nzb-1][nxb-1][nhb-2])+4.0*c3[nzb-2]*(tmp[nzb-2][nxb-2][nhb-2]);
   
	for (iz=0; iz < nzb; iz++) {
	    for (ix=0; ix < nxb; ix++) {
		for (ih=0; ih < nhb; ih++) {
                    tmpq[iz][ix][ih]  += (2.0*curr[iz][ix][ih]-prev[iz][ix][ih]);
                    prev[iz][ix][ih]  = curr[iz][ix][ih];
                    curr[iz][ix][ih]  = tmpq[iz][ix][ih];
                }
            }
        }
        bd3_decay(prev);
        bd3_decay(curr);
	if (!mig) { /* modeling -> write out data */
	    //cosft2(true,nh,nx,dat);
	    sf_floatwrite(dat[0],nh*nx,data);
	}
       sf_warning("              tmpq=          %f           ;",curr[102][145][90]);
    }
    sf_warning(".");

    if (mig) {
	//for (iz=1; iz < nz; iz++) {
	for (iz=0; iz < nz; iz++) {
	    for (ix=0; ix < nx; ix++) {
		    img[ix][iz] = curr[iz+nbt][ix+nxl][nhl];
	    }
	}

        sf_floatwrite(img[0],nz*nx,image);
    }
    bd3_close();
    exit(0);
}
