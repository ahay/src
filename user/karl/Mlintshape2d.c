#include <string.h>
#include <rsf.h>
#include <math.h>

/* global variables */
float* x;
float* y;
float* z;
int max_npnts;


void free2d(void** p){
  free(p[0]);
  free(p);
}

void bilinear(int adjoint, int add,
               float** a, float xorg, float yorg,
               float deltax, float deltay,
               int nx, int ny,
               float xp, float yp, 
               float* out){
    int ix, iy;
    float u, v;
    float rx,ry;
    
    if(!add){
       if(!adjoint)*out=0;
       else {
          for(iy=0; iy<ny; iy++){
             for(ix=0; ix<nx; ix++){
                a[iy][ix]=0.0;
             }
          }
       }
    }
    rx=(xp-xorg)/deltax;
    if(rx<0)rx=0;
    if(rx>nx-1)rx=nx-1;
    ry=(yp-yorg)/deltay;
    if(ry<0)ry=0;
    if(ry>ny-1)ry=ny-1;
    ix=(int)rx;
    iy=(int)ry;
    if(0){
      fprintf(stderr,"x=%f, y=%f, rx=%f, ry=%f, xorg=%f, deltax=%f, yorg=%f, deltay=%f\n",xp,yp,rx,ry, xorg,deltax,yorg,deltay);
    }
    u=1-(rx-ix);
    v=1-(ry-iy);
    if(0)fprintf(stderr,"before clip u=%f,ix=%d,v=%f,iy=%d\n",u,ix,v,iy);
    if(ix<0){
       ix=0;
       u=1.0;
    } else if(ix>=nx-1){
       ix=nx-1;
       u=1.0;
    }
    if(iy<0){
       iy=0;
       v=1.0;
    } else if(iy>=ny-1){
       iy=ny-1;
       v=1.0;
    }
    if(u>1.0)u=1.0;
    if(v>1.0)v=1.0;
    if(ix>nx-2){
        if(iy>ny-2){
            if(!adjoint)*out+=a[iy][ix];
            else a[iy][ix]+=*out;
        } else {
            if(!adjoint){
                *out+=u*v*a[iy][ix]+u*(1-v)*a[iy+1][ix  ];
            } else{
               a[iy  ][ix  ]+= u   *   v * *out;
               a[iy+1][ix  ]+= u   *(1-v)* *out;
            }
        }
    }else{
        if(iy>ny-2){
            if(!adjoint)*out+=u*v*a[iy][ix]+(1-u)*v*a[iy][ix+1];
            else{
               a[iy  ][ix  ]+= u   *   v * *out;
               a[iy  ][ix+1]+=(1-u)*   v * *out;
            }
        } else {
            if(!adjoint){
               *out+=u* v   *a[iy  ][ix  ]+(1-u)*   v *a[iy  ][ix+1]+
                     u*(1-v)*a[iy+1][ix  ]+(1-u)*(1-v)*a[iy+1][ix+1];
               if(0)fprintf(stderr," rx=%f,ix=%d,u=%f\n ry=%f,iy=%d,v=%f\n",
                               rx,ix,u,ry,iy,v);
            } else {
              a[iy  ][ix  ]+= u   *   v * *out;
              a[iy  ][ix+1]+=(1-u)*   v * *out;         
              a[iy+1][ix  ]+= u   *(1-v)* *out;
              a[iy+1][ix+1]+=(1-u)*(1-v)* *out;         
            }        
        }
    }
}

void bilinearop(int adjoint, int add,
                float xorg, float yorg,
                float deltax, float deltay,
                int nx, int ny, int npnt,
                float** a, float* z){
   int ipnt, ix, iy;
   
   if(!add){
      if(!adjoint){
         for(ipnt=0; ipnt<npnt; ipnt++)z[ipnt]=0;
      } else {
         for(iy=0; iy<ny; iy++){
            for(ix=0; ix<nx; ix++){
               a[iy][ix]=0.0;
            }                      
         }
      }
   }
   
   for(ipnt=0; ipnt<npnt; ipnt++){
      bilinear(adjoint, 1,
               a,xorg,yorg,
               deltax,deltay,
               nx,ny,
               x[ipnt],y[ipnt], &z[ipnt]);
   }
}


void boxsmooth(float* in, int n, int len,
               float* out){
    int lenh=len/2;
    int j;
    double sum=0.0;
 
    if(lenh>n-1)lenh=n-1;   
    for(j=0; j<lenh; j++){
        sum+=in[j];
    }
    for(j=0; j<n; j++){
       if(j+lenh<n)sum+=in[j+lenh];
       if(j-lenh-1>=0)sum-=in[j-lenh-1];
       out[j]=sum;
    }
}
void boxsmooth2d(int adj, 
                 int nx, int ny, 
                 int lenx, int leny, 
                 float** in, float** out){
   int i,j;
   float* xnorm;
   float* ynorm;
   float* xtemparray; 
   float* xtemparray1; 
   float* ytemparray;
   float* ytemparray1;
   fprintf(stderr,"in boxsmooth2d lenx=%d,leny=%d\n",lenx,leny);
   /* compute the xnorm and ynorm */
   xnorm=(float*)malloc(nx*sizeof(float));
   ynorm=(float*)malloc(ny*sizeof(float));
   xtemparray =(float*)malloc(nx*sizeof(float));
   xtemparray1=(float*)malloc(nx*sizeof(float));
   ytemparray =(float*)malloc(ny*sizeof(float));
   ytemparray1=(float*)malloc(ny*sizeof(float));
   for(i=0; i<nx; i++)xtemparray[i]=1.0;
   boxsmooth(xtemparray , nx, lenx  , xtemparray1); 
   boxsmooth(xtemparray1, nx, (int)(lenx/1.5+.5), xnorm);  
   for(i=0; i<nx; i++)xnorm[i]=1.0/xnorm[i];
   
   for(i=0; i<ny; i++)ytemparray[i]=1.0;
   boxsmooth(ytemparray , ny, leny  , ytemparray1); 
   boxsmooth(ytemparray1, ny, (int)(leny/1.5+.5), ynorm); 
   for(i=0; i<ny; i++)ynorm[i]=1.0/ynorm[i];
   
   if(!adj){
      fprintf(stderr,"!adj\n");
      /* smoothx, normx, smoothy, normy */
      for(j=0; j<ny; j++){
         boxsmooth(in[j]      ,nx,lenx  ,xtemparray);
         boxsmooth(xtemparray ,nx,(int)(lenx/1.5+.5),out[j]);
         for(i=0; i<nx; i++)out[j][i]*=xnorm[i];
      }
 
      for(i=0; i<nx; i++){
         for(j=0; j<ny; j++)ytemparray[j]=out[j][i];
         boxsmooth(ytemparray ,ny,leny  ,ytemparray1);
         boxsmooth(ytemparray1,ny,(int)(leny/1.5+.5),ytemparray );
         for(j=0; j<ny; j++)out[j][i]=ytemparray[j]*ynorm[j];
      }
   } else {
      fprintf(stderr,"adj\n");
      /* normy, smoothy, normx, smoothx */
      for(i=0; i<nx; i++){
         for(j=0; j<ny; j++)ytemparray[j]=out[j][i]*ynorm[j];
         boxsmooth(ytemparray ,ny,(int)(leny/1.5+.5),ytemparray1);
         boxsmooth(ytemparray1,ny,leny  ,ytemparray );
         for(j=0; j<ny; j++)in[j][i]=ytemparray[j];
      }
      for(j=0; j<ny; j++){
         for(i=0; i<nx; i++)xtemparray[i]=in[j][i]*xnorm[i];
         boxsmooth(xtemparray ,nx,(int)(lenx/1.5+.5),xtemparray1);
         boxsmooth(xtemparray1,nx,lenx  ,in[j]);               }
      }
   free((void*)xnorm);
   free((void*)ynorm);
   free((void*)xtemparray);
   free((void*)xtemparray1);
   free((void*)ytemparray);
   free((void*)ytemparray1);
}

double dot2d(float** a2d, float** b2d,int nx, int ny){
       int ix,iy;
       double dotprod;
       dotprod=0.0;
       for(iy=0; iy<ny; iy++)
          for(ix=0; ix<nx; ix++)
             dotprod+=a2d[iy][ix]*b2d[iy][ix];
       return dotprod;      
}
double dot1d(float* a1d,float* b1d,int npnt){
       int ipnt;
       double dotprod;
       dotprod=0.0;
       for(ipnt=0; ipnt<npnt; ipnt++)dotprod+=a1d[ipnt]*b1d[ipnt];
       return dotprod;
}

void cgsolve(float* z, int npnt,
             float** a, float xorg, float yorg,
             float deltax, float deltay,
             int nx, int ny, int smoothlen,
             int niter){
     float** delta_x;
     float** aatemp;
     float** z1;
     float* r;
     int ipnt, iter;
     int ix, iy;
     double beta, alpha, gamma, deltax_dot_z1,deltar_dot_z2;
     float* z2; /* no need to allocate since no left conditioner */
     float** s;
     float* deltar;
     
     r=(float*)malloc(npnt*sizeof(float));
     delta_x=(float**)sf_floatalloc2(nx,ny);
     aatemp=(float**)sf_floatalloc2(nx,ny);
     z1=(float**)sf_floatalloc2(nx,ny);
     s=(float**)sf_floatalloc2(nx,ny);
     deltar=(float*)malloc(npnt*sizeof(float));

     bilinearop(0,0,
                xorg,yorg,
                deltax,deltay,
                nx,ny,npnt,
                a,r);
     for(ipnt=0; ipnt<npnt; ipnt++){
        r[ipnt]=z[ipnt]-r[ipnt];
     }
     fprintf(stderr,"before loop r_dot_r=%e\n",dot1d(r,r,npnt));
 
     beta=0;
     for(iter=0; iter<niter; iter++){
         fprintf(stderr,"iter=%d\n",iter); 
         bilinearop(1,0,
                    xorg,yorg,
                    deltax,deltay,
                    nx,ny,npnt,
                    delta_x,r);
         boxsmooth2d(0,
                     nx,ny,
                     smoothlen,smoothlen,
                     delta_x,aatemp);
         boxsmooth2d(1,
                     nx,ny,
                     smoothlen,smoothlen,
                     z1,aatemp);
         deltax_dot_z1=dot2d(delta_x,z1,nx,ny);
         if(iter>0){
            beta=deltax_dot_z1/gamma;
         }
         gamma=deltax_dot_z1;
         for(iy=0; iy<ny; iy++)
            for(ix=0; ix<nx; ix++)s[iy][ix]=z1[iy][ix]+beta*s[iy][ix];
         for(ipnt=0; ipnt<npnt; ipnt++)deltar[ipnt]*=beta;
         bilinearop(0,1,
                    xorg,yorg,
                    deltax,deltay,
                    nx,ny,npnt,
                    z1,deltar);
         z2=deltar;
         deltar_dot_z2=dot1d(deltar,z2,npnt);
         if(-1e-45<deltar_dot_z2 && deltar_dot_z2<1e-45 ) break;
         alpha=gamma/deltar_dot_z2;
	 fprintf(stderr,"gamma=%e,alpha=%e,deltar_dot_z2=%e\n",
		         gamma   ,alpha    ,deltar_dot_z2);
         for(iy=0; iy<ny; iy++)
            for(ix=0; ix<nx; ix++)a[iy][ix]+=alpha*s[iy][ix];
         for(ipnt=0; ipnt<npnt; ipnt++)r[ipnt]-=alpha*z2[ipnt];
         fprintf(stderr,"after iteration r_dot_r=%e\n",dot1d(r,r,npnt));
     }
     fprintf(stderr,"need to free arrays\n");
     free(r);
     free2d((void**)delta_x);
     free2d((void**)aatemp);
     free2d((void**)z1);
     free2d((void**)s);
     free(deltar);

}

int main(int argc, char *argv[]){
  int verbose;
  sf_file in=NULL, out=NULL;
  int n1_in, n2_in;
  int i2_in;
  int i;
  float** a;
  int nx,ny;
  float xmin,xmax,ymin,ymax;
  float xmin_data,xmax_data,ymin_data,ymax_data,zmin_data,zmax_data;
  float dx,dy;
  int max_nx_ny;
  int smoothlen;
  int numsmooth;

  sf_init (argc,argv);
  
  /*****************************/
  /* initialize verbose switch */
  /*****************************/
  /* verbose flag controls ammount of print */
  /*( verbose=1 0 terse, 1 informative, 2 chatty, 3 debug ) */
  /* fprintf(stderr,"read verbose switch.  getint reads command line.\n"); */
  if(!sf_getint("verbose",&verbose))verbose=1;
  /* \n
     flag to control amount of print
     0 terse, 1 informative, 2 chatty, 3 debug
  */
  fprintf(stderr,"verbose=%d\n",verbose);
  
  /******************************************/
  /* input and output data are stdin/stdout */
  /******************************************/
  
  if(verbose>0)fprintf(stderr,"read infile name\n");  
  in = sf_input ("in");
  
  if(verbose>0)
    fprintf(stderr,"read outfile name (probably default to stdout\n");
  out = sf_output ("out");
  
  if (!sf_histint(in,"n1",&n1_in))
    sf_error("input data not define n1_traces");
  if (!sf_histint(in,"n2",&n2_in))
    sf_error("input data not define n1_traces");
  
  if(!sf_getfloat("xmin",&xmin))sf_error("xmin is required parameter\n");
  if(!sf_getfloat("xmax",&xmax))sf_error("xmax is required parameter\n");
  if(!sf_getfloat("ymin",&ymin))sf_error("ymin is required parameter\n");
  if(!sf_getfloat("ymax",&ymax))sf_error("ymax is required parameter\n");
  if(!sf_getfloat("dx"  ,&dx  ))sf_error("dx is required parameter\n");
  if(!sf_getfloat("dy"  ,&dy  ))sf_error("dy is required parameter\n");
  if(!sf_getint  ("nx"  ,&nx  ))sf_error("nx is required parameter\n");
  if(!sf_getint  ("ny"  ,&ny  ))sf_error("ny is required parameter\n");
  
  sf_putfloat(out,"o1",xmin);
  sf_putint  (out,"n1",nx);
  sf_putfloat(out,"d1",dx);
  
  sf_putfloat(out,"o2",ymin);
  sf_putint  (out,"n2",ny);
  sf_putfloat(out,"d2",dy);
  
  sf_fileflush(out,in);
  
  /* now read the x,y,z data */
  
  x=(float*)malloc(n2_in*sizeof(float));
  y=(float*)malloc(n2_in*sizeof(float));
  z=(float*)malloc(n2_in*sizeof(float));
  
  fprintf(stderr,"n2_in=%d\n",n2_in);
  for(i2_in=0; i2_in<n2_in; i2_in++){
    sf_floatread(&(x[i2_in]),1,in);
    sf_floatread(&(y[i2_in]),1,in);
    sf_floatread(&(z[i2_in]),1,in);
  }
  
  



  xmin_data=x[0]; xmax_data=x[0];
  ymin_data=y[0]; ymax_data=y[0];
  zmin_data=z[0]; zmax_data=z[0];
  for(i=0; i<n2_in; i++){   
    if(xmin_data>x[i])xmin_data=x[i]; if(xmax_data<x[i])xmax_data=x[i];
    if(ymin_data>y[i])ymin_data=y[i]; if(ymax_data<y[i])ymax_data=y[i];
    if(zmin_data>z[i])zmin_data=z[i]; if(zmax_data<z[i])zmax_data=z[i];
  }
  fprintf(stderr,"xmin_data=%f, xmax_data=%f\n",xmin_data,xmax_data);
  fprintf(stderr,"ymin_data=%f, ymax_data=%f\n",ymin_data,ymax_data);
  fprintf(stderr,"zmin_data=%f, zmax_data=%f\n",zmin_data,zmax_data);
  
  a=(float**)sf_floatalloc2(nx,ny);
  if(ny>nx)max_nx_ny=ny;
  else max_nx_ny=nx;
  
  smoothlen=2.0*max_nx_ny+1.0;
  numsmooth=0;
  for( ; smoothlen>0.0; smoothlen=(int)(smoothlen/1.5)){
    cgsolve(z, n2_in,
	    a,xmin,ymin,
	    dx, dy, 
	    nx, ny, smoothlen,
	    5);
    numsmooth++;
    if(0 && numsmooth==1)break;
  }
  sf_floatwrite(&(a[0][0]),nx*ny,out);
  
  exit(0);
}
