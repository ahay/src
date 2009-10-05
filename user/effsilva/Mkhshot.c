/* Kirchhoff shot migration */

#include <rsf.h>
#include "alias1.h"
#include "aarots.h"
#include "slice.h"
#include "migzrots.h"

/*
improvements do be done:
 i.   use the correct weight function,
 ii.  review anti-aliasing filtering.

Golden, 11/28/2007. Eduardo Filpo.
*/

int main (int argc, char ** argv)
{
   /* RSF variables */
   bool verb;
   float *dat=NULL, **datf=NULL, **image=NULL, **tr=NULL, **ts=NULL, **trs=NULL;
   float **px=NULL, **pz=NULL;
   sf_file Fin=NULL, Fout=NULL, Ftt=NULL;
   sf_axis a1,a2,a1t,a2t,a3t,ax,az;
   int     nf,ntaper,i1,i2,i2t,itt;
   float   fmax,df,theta,dtheta,tg,tgtap,tmin,xs,xr,eps;
   aalias  aa;
   fslice  tabtt=NULL;
   int     n1,n2,n1t,n2t,n3t,nx,nz;
   float   o1,d1,o2,d2,o1t,d1t,o2t,d2t,o3t,d3t,ox,oz,dx,dz;

   sf_init(argc,argv);

   Fin  = sf_input  ("in" );
   Fout = sf_output ("out" );
   Ftt  = sf_input  ("ttfile" );

   /* axes */
   a1  = sf_iaxa(Fin,1);   /* data */
   a2  = sf_iaxa(Fin,2);   /* data */
   a1t = sf_iaxa(Ftt,1);   /* traveltime */
   a2t = sf_iaxa(Ftt,2);   /* traveltime */
   a3t = sf_iaxa(Ftt,3);   /* traveltime */

   o1  = sf_o(a1);  n1   = sf_n(a1);  d1  = sf_d(a1);
   o2  = sf_o(a2);  n2   = sf_n(a2);  d2  = sf_d(a2);
   o1t = sf_o(a1t); n1t  = sf_n(a1t); d1t = sf_d(a1t);
   o2t = sf_o(a2t); n2t  = sf_n(a2t); d2t = sf_d(a2t);
   o3t = sf_o(a3t); n3t  = sf_n(a3t); d3t = sf_d(a3t);

   /* migration parameters */
   if(! sf_getbool("verb",&verb))     verb = false; /* verbosity flag */
   if(! sf_getfloat("theta",&theta))  theta = 30.;  /* maximum dip */
   if(! sf_getfloat("dtheta",&dtheta))  dtheta = theta/3;  /* taper zone */
   if(dtheta>theta) dtheta=theta;
   if(! sf_getfloat("df",&df))        df = 5.;    /* anti-aliasing sampling */
   if(!sf_getfloat("fmax",&fmax)) fmax=.5/d1;
   if(fmax>(.5/d1)) fmax=.5/d1;
   if (!sf_getint("ntaper",&ntaper)) ntaper=11;
   if(!sf_getfloat("tmin",&tmin)) tmin=3*d1;
   if(!sf_getfloat("xs",&xs)) sf_error("missing xs parameter\n");

   /* image parameters */
   if (!sf_getint("nx",&nx))   nx=n2t;
   if(!sf_getfloat("ox",&ox))  ox=o2t;
   if(!sf_getfloat("dx",&dx))  dx=d2t;
   if (!sf_getint("nz",&nz))   nz=n1t;
   if(!sf_getfloat("oz",&oz))  oz=o1t;
   if(!sf_getfloat("dz",&dz))  dz=d1t;

   /* checking dimensions */
   if((dx!=d2t)||(dz!=d1t)) 
     sf_error("sampling interval have to be the same in:\n"
	      " image and traveltime file\n");
   if(ox<o2t) ox=o2t; 
   if(oz<o1t) oz=o1t;
   if((ox+(nx-1)*dx)>(o2t+(n2t-1)*d2t)) nx=floor(((o2t+(n2t-1)*d2t)-ox)/dx)+1;
   if((oz+(nz-1)*dz)>(o1t+(n1t-1)*d1t)) nz=floor(((o1t+(n1t-1)*d1t)-oz)/dz)+1;

   /* output axis */
   ax = sf_maxa(nx,ox,dx); if(verb) sf_raxa(ax);
   az = sf_maxa(nz,oz,dz); if(verb) sf_raxa(az);
   sf_oaxa(Fout,az,1);
   sf_oaxa(Fout,ax,2);

   /* anti-aliasing */
   /* df = fmax; */
   nf = initAalias(-1,verb,fmax,df,n1,d1,&aa);
   /* fprintf(stderr,"forcing nf=%d df=%f\n",nf,df); */

   /* aperture angle */
   tg    = tan(SF_PI*theta/180);
   tgtap = tan(SF_PI*(theta-dtheta)/180);
   if(verb) sf_warning("tgmax=%f tgtap=%f",tg,tgtap);

   /* allocating */
   dat   = sf_floatalloc(n1);
   image = sf_floatalloc2(nz,nx);
   datf  = sf_floatalloc2 (n1,nf);
   ts    = sf_floatalloc2(n1t,n2t);
   tr    = sf_floatalloc2(n1t,n2t);
   trs   = sf_floatalloc2(n1t,n2t);
   px    = sf_floatalloc2(n1t,n2t);
   pz    = sf_floatalloc2(n1t,n2t);

   if(verb) sf_warning("initializing traveltime loading");
   /* initializing traveltime maps */
   tabtt = fslice_init(n1t*n2t,n3t,sizeof(float));
   fslice_load(Ftt,tabtt,SF_FLOAT);
   if(verb) sf_warning("traveltime loading has finished");

   /* reading the source-slice from traveltime table */
   eps = .01*d2;
   itt = floor((xs+eps-o3t)/d3t);
   fslice_get(tabtt,itt,ts[0]);
   if(verb) sf_warning("traveltime table from source was read");

   for(i2=0;i2<n2;i2++) {

     sf_floatread(dat,n1,Fin);
     xr = o2+i2*d2;

     if((xr>=o3t)&&(xr<(o3t+n3t*d3t))) {

       for (i1=n1-ntaper;i1<n1;i1++) dat[i1]=(n1-i1-1)*dat[i1]/ntaper;

       loadBank(aa,dat,datf);

       /* reading the receiver-slice of traveltime table */
       /* itt = floor((o2-o2t)/d2)+i2; */
       itt = floor((xr+eps-o3t)/d3t);
       fslice_get(tabtt,itt,tr[0]);

       for(i2t=0;i2t<n2t;i2t++)
         for(i1=0;i1<n1t;i1++)
            trs[i2t][i1]=ts[i2t][i1]+tr[i2t][i1];

       derive_1(n1t,n2t,d1t,trs,pz);

       derive_2(n1t,n2t,d2t,trs,px);

       spreadSR(nf,fmax,df,tg,tgtap,
                n1,o1,d1, nx,ox,dx, nz,oz,dz, o1t,o2t, 
                px,pz,tr,ts,datf,image);

       if(verb) fprintf(stderr,"+");

     }else{ if(verb) fprintf(stderr,".");}

   }

   sf_floatwrite(image[0],n1*n2,Fout);
   fprintf( stderr," \n finished processing \n");

 
   exit(0);
}
