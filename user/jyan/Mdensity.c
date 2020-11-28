/* Compute density */
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
#include <math.h>

static int inRegion(float val, float min, float max, int n)
{
    int ind=0;
    int i;
    float xs, xe;
    
    for (i=0;i<n; i++)    {
	xs=min+i*(max-min)/n;
	xe=min+(i+1)*(max-min)/n;
	if (val>xs && val<=xe) 	{
	    ind=i; 
	    break;
	}
	
    }
    return ind;
    
}

int main(int argc, char* argv[])
{
    bool verb;

    /* I/O files */
    sf_file Fin=NULL,FinY=NULL,FinZ=NULL,FinW=NULL; 
    sf_file Fout=NULL; 
    
  
    sf_axis az,ax,ay;      /* cube axes */
    sf_axis denax,denay,denaz,denaw;
    float ***inX,***inY,***inZ,***inW;              
    float ****out;
    
/*    float *tmpx,*tmpy,*tmpz,*tmpw;*/
   
    int ix,iy,iz,nx,nz,ny,nn,nnn;
    float dx,dz,dy;
    float dw;
 
    int n1,n2,n3,n4;

    float minx=1,maxx=0;
    float miny=1,maxy=0;
    float minz=90,maxz=-90;
    float minw=180,maxw=-180;
    
    int ii,jj,kk,ll;
    int indx,indy,indz,indw;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getint("n1",&n1)) n1=20; /* verbosity flag */
    if(! sf_getint("n2",&n2)) n2=20; /* verbosity flag */
    if(! sf_getint("n3",&n3)) n3=20; /* verbosity flag */
    if(! sf_getint("n4",&n4)) n4=20; /* verbosity flag */
    /*------------------------------------------------------------*/
    /* input files and axes */    
    Fin  = sf_input ("in");
    FinY = sf_input ("inY");
    FinZ = sf_input ("inZ");
    FinW = sf_input ("inW");

    az = sf_iaxa(Fin,1);
    ax = sf_iaxa(Fin,2);
    ay = sf_iaxa(Fin,3);
    nz=sf_n(az);
    nx=sf_n(ax);
    ny=sf_n(ay);

    dz=sf_d(az);
    dx=sf_d(ax);
    dy=sf_d(ay);

    nn=nx*ny*nz;

    inX=sf_floatalloc3(nz,nx,ny);   
    sf_floatread(inX[0][0],nn,Fin); 

    inY=sf_floatalloc3(nz,nx,ny);   
    sf_floatread(inY[0][0],nn,FinY);

    inZ=sf_floatalloc3(nz,nx,ny);   
    sf_floatread(inZ[0][0],nn,FinZ);

    inW=sf_floatalloc3(nz,nx,ny);   
    sf_floatread(inW[0][0],nn,FinW);

    if(verb) { 
	sf_raxa(az); 
	sf_raxa(ax); 
	sf_raxa(ay); 
    }
    
   
    /* sort values */
    
/*     smax = sf_quantile(ns-1,ns,ss2); */
/*     smin = sf_quantile(   0,ns,ss2); */
/*     nr = SF_MIN(nr,1+(smax-smin)/ds); */
   
    /*------------------------------------------------------------*/

    for (iz=0;iz<nz;iz++){
	for (ix=0;ix<nx;ix++){
	    for (iy=0;iy<ny;iy++){
		if (inX[iy][ix][iz] < minx ) minx=inX[iy][ix][iz];
		if (inY[iy][ix][iz] < miny ) miny=inY[iy][ix][iz];
		if (inZ[iy][ix][iz] < minz ) minz=inZ[iy][ix][iz];
		if (inW[iy][ix][iz] < minw ) minw=inW[iy][ix][iz];

		if (inX[iy][ix][iz] > maxx ) maxx=inX[iy][ix][iz];
		if (inY[iy][ix][iz] > maxy ) maxy=inY[iy][ix][iz];
		if (inZ[iy][ix][iz] > maxz ) maxz=inZ[iy][ix][iz];
		if (inW[iy][ix][iz] > maxw ) maxw=inW[iy][ix][iz];


	    }
	}
    }

    sf_warning("minx=%f maxx=%f",minx,maxx);
    sf_warning("miny=%f maxy=%f",miny,maxy);
    sf_warning("minz=%f maxz=%f",minz,maxz);
    sf_warning("minw=%f maxw=%f",minw,maxw);


    /* a four-dimensional problem    */

/*    int n=20;*/
    dx=(maxx-minx)/n1;
    dy=(maxy-miny)/n2;
    dz=(maxz-minz)/n3;
    dw=(maxw-minw)/n4;
    /*    /\*------------------------------------------------------------*\/ */
    /*     /\* output files*\/ */
    Fout = sf_output ("out");   
    out=sf_floatalloc4(n1,n2,n3,n4);
    nnn=n1*n2*n3*n4;
    denax=sf_maxa( n1 , minx , dx);
    denay=sf_maxa( n2 , miny , dy);
    denaz=sf_maxa( n3 , minz , dz);
    denaw=sf_maxa( n4 , minw , dw);

    sf_oaxa(Fout,denax,1);
    sf_oaxa(Fout,denay,2);
    sf_oaxa(Fout,denaz,3);
    sf_oaxa(Fout,denaw,4);
    
    for(ii=0;ii<n1;ii++){
	for(jj=0;jj<n2;jj++){
	    for(kk=0;kk<n3;kk++){
		for(ll=0;ll<n4;ll++){
		    out[ll][kk][jj][ii]=0;
		}
	    }
	}
    }

    for (iz=0;iz<nz;iz++){
	for (ix=0;ix<nx;ix++){
	    for (iy=0;iy<ny;iy++){
		indx=inRegion(inX[iy][ix][iz],minx,maxx,n1);
		indy=inRegion(inY[iy][ix][iz],miny,maxy,n2);
		indz=inRegion(inZ[iy][ix][iz],minz,maxz,n3);
		indw=inRegion(inW[iy][ix][iz],minw,maxw,n4);	
		out[indw][indz][indy][indx]+=1;
	    }
	}
    }

    for(ii=0;ii<n1;ii++){
	for(jj=0;jj<n2;jj++){
	    for(kk=0;kk<n3;kk++){
		for(ll=0;ll<n4;ll++){
		    out[ll][kk][jj][ii]/=nn;		    
		}
	    }
	}
    }



    sf_floatwrite(out[0][0][0],nnn,Fout);
    exit (0);
}
