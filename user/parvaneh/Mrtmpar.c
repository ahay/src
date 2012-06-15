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
    int ix2,iy1,iy2,nn,nm,nn1,nn2,i,k,jsnap,jreport;
    int ss1,ss2,kt,ii,ij,ik,l,m,n,nt,jk,first=1,last=0;
    float ***pnow=NULL,***lp=NULL,***v=NULL,***v2=NULL,*w=NULL, ***psave=NULL;
    float ***plast=NULL,***pnext=NULL,***slices1=NULL,***slices2=NULL,***slices3=NULL;
    float **slicexy,**slicezx,**slicezy;
    float *dampz,*dampx,*dampy,decayx,decayy,decayz,fpeak,dx,dy;
    float walpha[3],wbeta[3],zs,xs,ys,dz,dt,dt2;
    sf_file V1=NULL,Szx=NULL,Szy=NULL,Sxy=NULL;
    
    
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
    
    wavelet(w,dt,fpeak,itype,ntw);
    
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
    
    dampz=sf_floatalloc(nabsorbz);
    dampx=sf_floatalloc(nabsorbx);
    dampy=sf_floatalloc(nabsorby);

    getabsorb(dampz,nabsorbz,decayz);
    getabsorb(dampx,nabsorbx,decayx);
    getabsorb(dampy,nabsorby,decayy);

   
       
    V1=sf_input("in");
    v=sf_floatalloc3(ny,nx,nz);
    sf_floatread(v[0][0],ny*nx*nz,V1);
    v2=sf_floatalloc3(nye,nxe,nze);
    nn=nx*ny*nz;
    nm=nx*ny;
    nn1=nx*nz;
    nn2=ny*nz;
    Szx=sf_output("out");
   
    
    
    for(i=0;i<nye;i++) {
	for(j=0;j<nxe;j++) {
	    for(k=0;k<nze;k++){
	    
		v2[k][j][i]=0;
	    }
	}
    }

       
    extend(v,v2,nz,nx,ny,iz1,iz2,ix1,ix2,iy1,iy2,nze,nxe,nye);
    free(**v);
        
    dt2=dt*dt;
    for(i=0;i<nye;i++) {
	for(j=0;j<nxe;j++) {
	    for(k=0;k<nze;k++){
	    
		v2[k][j][i]=(v2[k][j][i]*v2[k][j][i])*dt2;
	    }
	}
    }
    
    psave=sf_floatalloc3(nye,nxe,nze);
    pnow=sf_floatalloc3(nye,nxe,nze);
    plast=sf_floatalloc3(nye,nxe,nze);
    pnext=sf_floatalloc3(nye,nxe,nze);

    lp=sf_floatalloc3(nye,nxe,nze);
    slicexy=sf_floatalloc2(ny,nx);
    slicezx=sf_floatalloc2(nx,nz);
    slicezy=sf_floatalloc2(ny,nz);

    jsnap=0;
    jreport=0;

    ss1=(int)(nttotal/nsnap);
    ss2=(int)(nttotal/nreport);
    sf_putint(Szx,"n3",ss2);
    slices1=sf_floatalloc3(ss2,ny,nx);
    slices2=sf_floatalloc3(ss2,nx,nz);
    slices3=sf_floatalloc3(ss2,ny,nz);

    
    for(ii=0;ii<nye;ii++) {
	for(ij=0;ij<nxe;ij++) {
	    for(ik=0;ik<nze;ik++) {
		psave[ik][ij][ii]=0;
		pnow[ik][ij][ii]=0;
		plast[ik][ij][ii]=0;
		pnext[ik][ij][ii]=0;
		lp[ik][ij][ii]=0;
	    }
	}
    }



    for(kt=0;kt<nttotal;kt++) {
	for(ii=0;ii<nye;ii++) {
	    for(ij=0;ij<nxe;ij++) {
		for(ik=0;ik<nze;ik++) {
		    psave[ik][ij][ii]=pnow[ik][ij][ii];
		    pnow[ik][ij][ii]=plast[ik][ij][ii];
		    plast[ik][ij][ii]=psave[ik][ij][ii];
		    pnext[ik][ij][ii]=plast[ik][ij][ii];
		}
	    }
	}

	if(kt<=ntw){
	    sourcetogrid(pnow,nze,nxe,nye,dz,dx,dy,iz1,iz2,ix1,ix2,iy1,iy2,zs,xs,ys,w[kt]);
	}
	firlaplace(pnow,lp,walpha,wbeta,ncoefz,ncoefy,ncoefx,dz,dx,dy,nze,nxe,nye,first,last);
	absorb(lp,nze,nxe,nye,dampz,nabsorbz,dampx,nabsorbx,dampy,nabsorby);

	for(l=0;l<nye;l++) {
	    for(m=0;m<nxe;m++) {
		for(n=0;n<nze;n++) {
		  lp[n][m][l]=v2[n][m][l]*lp[n][m][l];
		  pnext[n][m][l]=(2*pnow[n][m][l])-plast[n][m][l]+lp[n][m][l];
		}
	    }
	}
	
	if(kt%nreport==0) {
	    xycollect(pnext,slicexy,nze,nxe,nye,nx,ny,ix1,iy1,izs);
	    for(ii=0;ii<ny;ii++) {
		for(ij=0;ij<nx;ij++)
		slices1[ij][ii][jreport]=slicexy[ij][ii];
	    jreport=jreport+1;
	    }
	}
	
	if(kt%nsnap==0) {
	    zxcollect(pnext,slicezx,nze,nxe,nye,nz,nx,iz1,ix1,iys);
	    for(i=0;i<nx;i++) {
		for(j=0;j<nz;j++)
		    slices2[j][i][jsnap]=slicezx[j][i];
	    }

	    zycollect(pnext,slicezy,nze,nxe,nye,nz,ny,iz1,iy1,ixs);
	    for(ik=0;ik<ny;ik++) {
		for(jk=0;jk<nz;jk++)
		    slices3[jk][ik][jsnap]=slicezy[ik][jk];
	    }
	    jsnap=jsnap+1;
	}
    }
    
    
    sf_floatwrite(slices2[0][0],ss2*nx*nz,Szx);
   
    
   

    nt=jreport;
   

    /* sei=sf_floatalloc2(nt,nx*ny);
       temp=sf_floatalloc2(nx*ny,nt);

       trans(seis,temp,nx*ny,nt);
    sf_floatwrite(seis[0],nt*nm,Slices);*/

    free(**pnow);
    free(**psave);
    free(**plast);
    free(**pnext);
    
    exit(0);

    
	
   

}
	    
