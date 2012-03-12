/* RTM */

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

#include <rsf.h>
#include <math.h>
#include "rtmsub.h"

#define max(a,b) (a>b) ? (a) : (b)
#define min(a,b) (a<b) ? (a) : (b)

int main(int argc, char* argv[])
{
   
    int ncoefz,ncoefx,ncoefy,nz,nx,ny,izs,ixs,iys,nreport,nsnap,ntw,nze;
    int nye,nxe,nabsorbz,nabsorby,nabsorbx,itype,nttotal,j,iz1,iz2,ix1;
    int ix2,iy1,iy2,norig,ngrid,nn,nm,nn1,nn2,i,k,jsnap,jreport;
    int ss1,ss2,kt,ii,ij,ik,l,m,n,nt,jk,first=1,last=0;
    float ***pnow=NULL,***lp=NULL,***v=NULL,***v2=NULL,*w=NULL, ***psave=NULL;
    float ***plast=NULL,***pnext=NULL,***slices1=NULL,***slices2=NULL,***slices3=NULL;
    float **zx=NULL,**zy=NULL,**xy=NULL,**slicexy,**slicezx,**slicezy;
    float *dampz,*dampx,*dampy,*slice,*seis,**temp,decayx,decayy,decayz,fpeak,dx,dy;
    float walpha[3],wbeta[3],zs,xs,ys,dz,dt,dt2;
    sf_file V1=NULL,Szx=NULL;
    const int iz=1,ix=2,iy=3;
    
    sf_init(argc,argv);
   

    if(!sf_getint("nz",&nz)) sf_error("need nz=");
    if(!sf_getint("nx",&nx)) sf_error("need nx=");
    if(!sf_getint("ny",&ny)) sf_error("need ny=");
    if(!sf_getint("ncoefz",&ncoefz)) sf_error("need nz=");
    if(!sf_getint("ncoefx",&ncoefx)) sf_error("need nz=");
    if(!sf_getint("ncoefy",&ncoefy)) sf_error("need nz=");
    if(!sf_getint("itype",&itype)) sf_error("need nz=");
    if(!sf_getint("ntw",&ntw)) sf_error("need nz=");
    if(!sf_getint("nttotal",&nttotal)) sf_error("need nz=");
    if(!sf_getint("nreport",&nreport)) sf_error("need nz=");
    if(!sf_getint("nsnap",&nsnap)) sf_error("need nz=");
    

    
    if(!sf_getfloat("zs",&zs)) sf_error("need nz=");
    if(!sf_getfloat("xs",&xs)) sf_error("need nx=");
    if(!sf_getfloat("ys",&ys)) sf_error("need ny=");
    if(!sf_getfloat("dz",&dz)) sf_error("need nz=");
    if(!sf_getfloat("dx",&dx)) sf_error("need nx=");
    if(!sf_getfloat("dy",&dy)) sf_error("need ny=");
    if(!sf_getfloat("decayz",&decayz)) sf_error("need nz=");
    if(!sf_getfloat("decayx",&decayx)) sf_error("need nx=");
    if(!sf_getfloat("decayy",&decayy)) sf_error("need ny=");
    if(!sf_getfloat("dt",&dt)) sf_error("need nz=");
    if(!sf_getfloat("fpeak",&fpeak)) sf_error("need nz=");
    if(!sf_getint("nabsorbz",&nabsorbz)) sf_error("need nz=");
    if(!sf_getint("nabsorbx",&nabsorbx)) sf_error("need nx=");
    if(!sf_getint("nabsorby",&nabsorby)) sf_error("need ny=");

    
    
    izs=max(1,(int)(zs/dz)+1);
    izs=min(izs,nz);
    ixs=max(1,(int)(xs/dx)+1);
    ixs=min(ixs,nx);
    iys=max(1,(int)(ys/dy)+1);
    iys=min(iys,ny);
    nreport=max(nreport,1);
    nsnap=max(nsnap,1);

    for(j=0;j<3;j++){
	walpha[j]=0.56;
	wbeta[j]=6;
    }
    if(ny<=1) {
	ny=1;
	dy=1;
	nye=1;
	ys=0;
	
	decayy=0;
	ncoefy=1;
    }
    
    w=sf_floatalloc(ntw);
    
    wavelet(w,ntw,dt,fpeak,itype);
    dampz=sf_floatalloc(nabsorbz);
    dampx=sf_floatalloc(nabsorbx);
    dampy=sf_floatalloc(nabsorby);

    getabsorb(dampz,nabsorbz,decayz);
    getabsorb(dampx,nabsorbx,decayx);
    getabsorb(dampy,nabsorby,decayy);

    nze=nz+2*nabsorbz;
    nxe=nx+2*nabsorbx;
    nye=ny+2*nabsorby;

    

    iz1=nabsorbz+1;
    iz2=iz1+nz-1;
    ix1=nabsorbx+1;
    ix2=ix1+nx-1;
    iy1=nabsorby+1;
    iy2=iy1+ny-1;

    
    if(ny<=1) {
	iy1=1;
	iy2=1;
	nye=1;
	nabsorby=1;
    }

    norig=nz*nx*ny;
    ngrid=nze*nxe*nye;
    Szx=sf_output("out");
    
    V1=sf_input("in");
    v=sf_floatalloc3(ny,nx,nz);
    sf_floatread(v[0][0],ny*nx*nz,V1);
    v2=sf_floatalloc3(nye,nxe,nze);
   
    nn=nx*ny*nz;
    nm=nx*ny;
    nn1=nx*nz;
    nn2=ny*nz;
    
    for(i=0;i<nye;i++) {
	for(j=0;j<nxe;j++) {
	    for(k=0;k<nze;k++){
	    
		v2[k][j][i]=0;
	    }
	}
    }

       
    
    sf_floatwrite(v2[0][0],nye*nxe*nze,Szx);
    exit(0);

    
	
   

}
	    
