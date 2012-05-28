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

#include "fftc3.h"
#include "abcpass3b.h"


int main(int argc, char* argv[])
{
    bool mig;
    int it, nt, ix, nx, iz, nz, ih, nh, it1, it2, its, snap, is, ir;
    int nx2, nz2, nh2, nk;
    float kz0, kx0, kh0;
    int nzb, nbt, nbb, nxb, nxl, nxr, nhb, nhr, nhl;
    float ct, cb, chl, chr, cxl, cxr;
    float dt, dx, dz, dh, kx, kz, kh,  h, x, z, dkx, dkz, dkh, pi=SF_PI;
    float ***prev, ***curr, **img, **dat, ***tmp, ***tmpq, **v;
    float ***a, ***b1, ***b2, ***b3, ***d1, ***d2, ***d3, epsilon;
    float x1, x2, x3, v0, v02, vz2, t2, vs, vr, vp2, tmpdt;
    kiss_fft_cpx ***uk;
    float pp;
    sf_file data, image, snaps, vel; 
    int nvz;
    float dvz; 
    float ww, gg;

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

    if (!sf_getint("nbt",&nbt)) nbt=0;
    if (!sf_getint("nbb",&nbb)) nbb=0;
    if (!sf_getint("nxl",&nxl)) nxl=0;
    if (!sf_getint("nxr",&nxr)) nxr=0;
    if (!sf_getint("nhl",&nhl)) nhl=0;
    if (!sf_getint("nhr",&nhr)) nhr=0;

    if (!sf_getfloat("ct",&ct)) ct = 0.0; /*decaying parameter*/
    if (!sf_getfloat("cb",&cb)) cb = 0.0; /*decaying parameter*/
    if (!sf_getfloat("chl",&chl)) chl = 0.0; /*decaying parameter*/
    if (!sf_getfloat("chr",&chr)) chr = 0.0; /*decaying parameter*/
    if (!sf_getfloat("cxl",&cxl)) cxl = 0.0; /*decaying parameter*/
    if (!sf_getfloat("cxr",&cxr)) cxr = 0.0; /*decaying parameter*/

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
    uk  =  (kiss_fft_cpx ***) sf_complexalloc3(nhb,nxb,nzb); 
    

    /* velocity */

    v = sf_floatalloc2(nzb,nxb);
    for (ix=nxl; ix < nx+nxl; ix++) {
        sf_floatread(v[ix]+nbt,nz,vel);
    }
    for (ix=0; ix < nxl; ix++) {
	for (iz=nbt; iz < nbt+nz; iz++) {
	    v[ix][iz] = v[nxl][iz];
        }
    }
    for (ix=nxl+nx; ix < nxb; ix++) {
	for (iz=nbt; iz < nbt+nz; iz++) {
	    v[ix][iz] = v[nxl+nx-1][iz];
        }
    }
    for (ix=0; ix < nxb; ix++) {
	for (iz=0; iz < nbt; iz++) {
            v[ix][iz] = v[ix][nbt];
	}
	for (iz=0; iz < nbb; iz++) {
	    v[ix][nz+nbt+iz] = v[ix][nz+nbt-1];
	}
    }
    v02 =0.0;
    for (ix=0; ix < nxb; ix++) {
        for (iz=0; iz < nzb; iz++) {
            v02 += v[ix][iz]*v[ix][iz];
        }
    }
    v02 /= (8.0*nxb*nzb);
    v0 = sqrtf(v02);
    sf_warning("v0=%f",v0);

   /* dkx = cosft_dk(nxb,dx);
    dkz = cosft_dk(nzb,dz);
    dkh = cosft_dk(nhb,dh);
   */
    dkx = 1./(nxb*dx);
    kx0 = -0.5/dx;
    dkh = 1./(nhb*dh);
    kh0 = -0.5/dh;
    dkz = 1./(nzb*dz);
    kz0 = -0.5/dz;

    a  = sf_floatalloc3(nhb,nxb,nzb);
    b1 = sf_floatalloc3(nhb,nxb,nzb);
    b2 = sf_floatalloc3(nhb,nxb,nzb);
    b3 = sf_floatalloc3(nhb,nxb,nzb);
    d1 = sf_floatalloc3(nhb,nxb,nzb);
    d2 = sf_floatalloc3(nhb,nxb,nzb);
    d3 = sf_floatalloc3(nhb,nxb,nzb);
    x1 = dh*dh;
    x2 = dx*dx;
    x3 = dz*dz;
    t2 = dt*dt;
  
    for (iz=0; iz < nzb; iz++) {
	for (ix=0; ix < nxb; ix++) {
            for (ih=0; ih < nhb; ih++) {
	        is = SF_MIN(SF_MAX(0,ix-ih),nxb-1);
	        ir = SF_MIN(SF_MAX(0,ix+ih),nxb-1);
	        vs = v[is][iz];  
	        vr = v[ir][iz]; 
	        vp2 = 0.5*(vs+vr);
	        vp2 *= vp2;
    	        vs *= vs;
	        vr *= vr;
                vz2= vs*vr/(8.0*vp2);
                ww=vz2/v02;
                gg=ww/72.0*((v02-vz2)*t2);
                pp=(t2*vz2)*(t2*(vz2*ww/90.0+v02/60.0-vz2/36.0));
                d1[iz][ix][ih] = pp/(x1*x1)+gg/x1;
                d2[iz][ix][ih] = pp/(x2*x2)+gg/x2;
                d3[iz][ix][ih] = pp/(x3*x3)+gg/x3;
                b1[iz][ix][ih] = -12.0*gg/x1-4.0*d1[iz][ix][ih];
                b2[iz][ix][ih] = -12.0*gg/x2-4.0*d2[iz][ix][ih];
                b3[iz][ix][ih] = -12.0*gg/x3-4.0*d3[iz][ix][ih];
                a[iz][ix][ih]  = ww-2.0*(b1[iz][ix][ih]+b2[iz][ix][ih]+b3[iz][ix][ih]+d1[iz][ix][ih]+d2[iz][ix][ih]+d3[iz][ix][ih]);
            }
        }
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
	for (ix=0; ix < nxb; ix++) {
	    for (iz=0; iz < nzb; iz++) {
		img[ix][iz] = 0;
            }
        }
	for (ix=nxl; ix < nx+nxl; ix++) {
            sf_floatread(img[ix]+nbt,nz,image);
	}
	/* Initialize model */
	for (iz=0; iz < nzb; iz++) {
	    for (ix=0; ix < nxb; ix++) {
		curr[iz][ix][nhl] = img[ix][iz];
	    }
	}
	
	/* step forward in time */
	it1 = 0;
	it2 = nt;
	its = +1;
    }

    
    nk = fft3_init(true,1,nhb,nxb,nzb,&nh2,&nx2,&nz2);
    if((nh2!=nhb) || (nx2!=nxb) || (nz2!=nzb)) sf_warning("not fft number");

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


	/* Initialize model */

	for (ix=0; ix < nx; ix++) {
            for (ih=0; ih < nh; ih++) {
		if (mig) {
                    curr[nbt][ix+nxl][ih+nhl] += dat[ix][ih];
		    } else {
			dat[ix][ih] = curr[nbt][ix+nxl][ih+nhl];
		    }
            }
        }
	if (NULL != snaps && 0 == it%snap) 
            for (iz=nbt; iz < nbt+nz; iz++) {
                for (ix=nxl; ix < nxl+nx; ix++) {
	            sf_floatwrite(prev[iz][ix]+nhl,nh,snaps);
	            //sf_floatwrite(curr[iz][ix]+nhl,nh,snaps);
                }
            }
/*
        for (iz=0; iz < nzb; iz++) {
            for (ix=0; ix < nxb; ix++) {
       	  for (ih=0; ih < nhb; ih++) {
                    tmp[iz][ix][ih] = curr[iz][ix][ih];
                }
            }
        }
 */        
      //  cosft3(false,nhb,nxb,nzb,tmp);
	fft3(curr,(sf_complex ***)uk);

	for (iz=0; iz < nzb; iz++) {
	    kz = (kz0+iz*dkz)*2.0*pi;
             z = kz*kz;
	    for (ix=0; ix < nxb; ix++) {
		kx = (kx0+ix*dkx)*2.0*pi;
                 x = kx*kx;
		for (ih=0; ih < nhb; ih++) {
		    kh = (kh0+ih*dkh)*2.0*pi;
                     h = kh*kh;
		     //c = tmp[iz][ix][ih];
                    tmpdt = 2*(cosf(v0*dt*sqrtf(x+h+z+sqrtf(4.0*x*h+(x+h+z)*(x+h+z))))-1); 
//#ifdef SF_HAS_COMPLEX_H
//                    uk[iz][ix][ih] = uk[iz][ix][ih]*tmpdt;
//#else
                    uk[iz][ix][ih] = sf_crmul(uk[iz][ix][ih],tmpdt);
//#endif
//		    tmp[iz][ix][ih] = 2*(cosf(v0*dt*sqrtf(x+h+z+sqrtf(4.0*x*h+(x+h+z)*(x+h+z))))-1)*c;
		}
	    }
	}
        
	ifft3(tmp,(sf_complex ***) uk);
//	cosft3(true,nhb,nxb,nzb,tmp);
	for (iz=0; iz < nzb; iz++) {
	    for (ix=0; ix < nxb; ix++) {
		for (ih=0; ih < nhb; ih++) {
                    tmpq[iz][ix][ih] = 0.0;
                }
            }
        }
	for (iz=2; iz < nzb-2; iz++) {
	    for (ix=2; ix < nxb-2; ix++) {
		for (ih=2; ih < nhb-2; ih++) {
                    tmpq[iz][ix][ih] = a[iz][ix][ih]*tmp[iz][ix][ih] + b1[iz][ix][ih]*(tmp[iz][ix][ih-1]+tmp[iz][ix][ih+1])+b2[iz][ix][ih]*(tmp[iz][ix-1][ih]+tmp[iz][ix+1][ih])+b3[iz][ix][ih]*(tmp[iz+1][ix][ih]+tmp[iz-1][ix][ih])+ d1[iz][ix][ih]*(tmp[iz][ix][ih-2]+tmp[iz][ix][ih+2])+d2[iz][ix][ih]*(tmp[iz][ix-2][ih]+tmp[iz][ix+2][ih])+d3[iz][ix][ih]*(tmp[iz+2][ix][ih]+tmp[iz-2][ix][ih]);
                }
            }
        }
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
	    sf_floatwrite(dat[0],nh*nx,data);
	}
    }
    sf_warning(".");

    if (mig) {
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
