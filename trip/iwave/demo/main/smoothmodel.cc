/* Author: Mario Bencomo
 * Purpose: Build 2D smooth models.
 * Modeled after standardmodel.cc
 */

#include <par.h>
#include <parser.h>
#include <math.h>

/* define number of model types */
#define NMODEL 4

#define cmax_def 4
#define cmin_def 2
#define rho_def  2.3

#define o1_def 0
#define o2_def 0
#define o3_def 0
#define d1_def 5
#define d2_def 5
#define d3_def 5
#define n1_def 101
#define n2_def 101
#define n3_def 1

#define cycles_def 2

#ifndef PI
#define PI 3.14159265359
#endif

/*********************** self documentation **********************/
const char *sdoc[] = {
    "                                                                ",
    " SMOOTHMODEL - build 2D smooth models",
    "                                                                ",
    " This module builds some standard 2D models.                    ",
    " Units:		                                             ",
    " 	densities    [g/cm^3]                                        ",
    "   velocities   [m/ms] = [km/s]                                 ",
    " 	bulk moduli  [GPa]  = [g/cm^3][m/ms]^2                       ",
    "   buoyancies   [cm^3/g]                                        ",
    "                                                                ",
    " The output format is native binary floats.                     ",
    "                                                                ",
    " Example run:   smoothmodel model=4 choose=2 > vp.bin           ",
    "          or:   smoothmodel model=4 choose=2 hfile=vp.rsf       ",
    "                                                                ",
    " Optional parameters:                                           ",
    "   model=1                    = choice of standard model        ",
    "                              =0; homogeneous                   ",
    "                              =1; linear velocity in depth, with",
    "                                  const density                 ",
    "                              =2; negative gaussian lense in    ",
    "                                  velocity, with const density  ",
    "                              =3; Sinusodial velocity in depth, ",
    "                                  with const density            ",
    "                                                                ",
    "   choose=1                   =0; output vp                     ",
    "                              =1; output bulk modulus           ",
    "                              =2; output density                ",
    "                              =3; output buoyancy               ",
    "                                                                ",
    "   hfile=\"\"                 =<hfile>: output in rsf/sep format",
    "                                  consisting of header file     ",
    "                                  hfile containing grid info,   ",
    "                                  and binary data file. If null ",
    "                                  string, output to stdout. If  ",
    "                                  non-null hfile given, path to ",
    "                                  data file is of form          ",
    "                                                                ",
    "                                  DATAPATH<hfile>@              ",
    "                                                                ",
    "                                  with DATAPATH read from env   ",
    "                                  (must end with / if set).     ",
    "                                                                ",
    "                                  Typical example: hfile=vp.rsf,",
    "                                  DATAPATH=/var/tmp/, so data   ",
    "                                  written to /var/tmp/vp.rsf@.  ",
    " ",
    "                                                                ",
    " Fast Dimension (depth)                                         ",
    "   o1=0.                      starting coordinate on grid       ",
    "   d1=5.0                     increment on grid in meters       ",
    "   e1=500.0                   end coordinate on grid            ",
    "   n1=(int)(1.5+(e1-o1)/d1)   number of gridpoints              ",
    "                                                                ",
    " Middle Dimension (x)                                           ",
    "   o2=0.                      starting coordinate on grid       ",
    "   d2=5.0                     increment on grid in meters       ",
    "   e2=500.0                   end coordinate on grid            ",
    "   n2=(int)(1.5+(e2-o2)/d2)   number of gridpoints              ",
    "                                                                ",
    " Slow Dimension (y)                                             ",
    "   o3=0.                      starting coordinate on grid       ",
    "   d3=5.0                     increment on grid in meters       ",
    "   e3=500.0                   end coordinate on grid            ",
    "   n3=(int)(1.5+(e3-o3)/d3)   number of gridpoints              ",
    "                                                                ",
    NULL};


//----------------------------------------------------------------------------//
// model 0: Homogeneous model
static inline float homogeneous(int choose, float cmin, float cmax, float rho)
{
    float vel = (cmax+cmin)/2;
    
    if (choose==1) return vel*vel*rho; 
    if (choose==2) return rho;
    if (choose==3) return 1./rho;
    
    return vel;
}
//----------------------------------------------------------------------------//
// model 1: linear depth velocity model with constant density
static inline float lineardepth(float x1, int choose, 
				float cmin, float cmax, 
				float rho, float * modelpars)
{
    float top = modelpars[0];
    float bot = modelpars[1];
    
    float vel = (bot-x1)/(bot-top)*cmin + (x1-top)/(bot-top)*cmax;
    
    if (choose==1) return vel*vel*rho;
    if (choose==2) return rho;
    if (choose==3) return 1./rho;
        
    return vel;
}
//----------------------------------------------------------------------------//
// model 2: negative gaussian lense in velocity, with constant density
static inline float gaussianlense2D(float x1, float x2, int choose, 
				    float cmin, float cmax, float rho, float * modelpars)
{
    float cx1 = modelpars[0];
    float cx2 = modelpars[1];
    float sigmax1 = modelpars[2];
    float sigmax2 = modelpars[3];
    
    float gaussx1 = exp( -(x1-cx1)*(x1-cx1)/(2.*sigmax1*sigmax1) );
    float gaussx2 = exp( -(x2-cx2)*(x2-cx2)/(2.*sigmax2*sigmax2) );
    
    float vel = cmax - (cmax-cmin)*gaussx1*gaussx2;
    
    if (choose==1) return vel*vel*rho;
    if (choose==2) return rho;
    if (choose==3) return 1./rho;
    
    return vel;
}
//----------------------------------------------------------------------------//
// model 3: Sinusodial-depth velocity with constant density
static inline float Sindepth(float x1, int choose, 
			     float cmin, float cmax, 
			     float rho, float * modelpars)
{
    float top = modelpars[0];
    float bot = modelpars[1];
    float cycles = modelpars[2];
    
    float wnum = cycles*2*PI;
    float vel = (cmax-cmin)*.5*sin(wnum*(x1-bot)/(bot-top)) + (cmax+cmin)*.5;
    
    if (choose==1) return vel*vel*rho;
    if (choose==2) return rho;
    if (choose==3) return 1./rho;
    
    return vel;
}




/******************************************************************************/
float get_zvalue(float x1, float x2, float x3,
                 int model, int choose, 
		 float cmin, float cmax, float rho,
		 float * modelpars)
{
    float v;
	
    switch (model){
    case 0: v = homogeneous(choose,cmin,cmax,rho); break;
    case 1: v = lineardepth(x1,choose,cmin,cmax,rho,modelpars); break;
    case 2: v = gaussianlense2D(x1,x2,choose,cmin,cmax,rho,modelpars); break;
    case 3: v = Sindepth(x1,choose,cmin,cmax,rho,modelpars); break;
    default: v = homogeneous(choose,cmin,cmax,rho);
    }
    return v;
}


/******************** end self doc **************/

int writemodel(
               int choose,     /* =1 output density; =2 output vp =3 output vm */
               int model,      /* choice of model                              */
               int n1,         /* number of samples in depth (fast dimension)  */
               int n2,         /* number of samples in x (middle dimension)    */
               int n3,         /* number of samples in y (slow dimenssion)     */
               float o1,       /* start depth                                  */
               float o2,       /* start x                                      */
               float o3,       /* start y                                      */
               float e1,       /* end depth                                    */
               float e2,       /* end x                                        */
               float e3,       /* end y                                        */
               float d1,       /* depth increment                              */
               float d2,       /* x increment                                  */
               float d3,       /* y increment                                  */
               float cmin,     /* min velocity                                 */
	       float cmax,     /* max velocity                                 */
	       float rho,
	       FILE * fp,      /* output stream                                */
               float * omv     /* other model values                           */
)
{
    int j1, j2, j3;
    float x1, x2, x3;
    float *v   = NULL;
           
    v   = (float *)emalloc(n1*sizeof(float));
    
    for(j3 = 0; j3 < n3; j3++){
        x3 = o3 + j3 * d3;
        for(j2 = 0; j2 < n2; j2++){
            x2 = o2 + j2 * d2;
            for(j1 = 0; j1 < n1; j1++){
                x1 = o1 + j1 * d1;

                v[j1] = get_zvalue(x1, x2, x3, model, choose, cmin, cmax, rho, omv);
            }
            if (fwrite(v, sizeof(float), n1, fp) != n1){
                fprintf(stderr, "write error\n");
                free(v);
                return 1;
            }
        }
    }
    free(v);
    return 0;
}

/******************************************************************************/
int main(int argc, char **argv) {
    
    float o1,o2,o3,d1,d2,d3,e1,e2,e3;
    int n1,n2,n3,model,choose;
    float * omv; //other model values
    float cmin, cmax, rho;

    /* WWS */
    char * fname;
    char * dname;
    FILE * fp;
    /* end WWS */
    
    //    int i,j;
    
    char * cwdpath;
    char * pathptr;
    
    
    /******************
     * get parameters
     ******************/
    
    PARARRAY * par = ps_new();
    
    if (ps_createargs(par,argc-1,argv+1)) {
        fprintf(stderr,
                "ERROR. could not process command line\n");
        exit(1);
    }
    
    /***********************
     * end get parameters
     ***********************/
    
    xargc=argc; xargv=argv;
    requestdoc(0);
    
    /*   if (ps_getparint("model",&model)) { */
    if (!(ps_ffint(*par,"model",&model))) {
        if (model<0 || model >= NMODEL) {
            fprintf(stderr,"Error: standardmodel.x\n");
            fprintf(stderr,"model index %d not defined: must lie in range [1,%d]\n",model,NMODEL);
            exit(1);
        }
    }
    else {
        fprintf(stdout,"Warning: standardmodel.x\n");
        fprintf(stdout,"no model index given, so using model=1\n");
        model=1;
    }
    /*   if(ps_getparint("choose",&choose)) { */
    if (!(ps_ffint(*par,"choose",&choose))) {
        if (choose<0 || choose > 3) {
            fprintf(stderr,"Error: smoothmodel.x\n");
            fprintf(stderr,"choose index must be 0 (veloxity),  1 (bulkmod), or 2 (density), or 3 (buoyancy) \n");
            exit(1);
        }
    }
    else {
        fprintf(stdout,"Warning: smoothmodel.x\n");
        fprintf(stdout,"no choose index given, so using choose=1 (density)\n");
        choose=1;
    }
    
    if(ps_fffloat(*par,"o1",&o1)) {o1=o1_def; fprintf(stdout, "WARNING: o1 is set to the default value!\n"); }
    if(ps_fffloat(*par,"o2",&o2)) {o2=o2_def; fprintf(stdout, "WARNING: o2 is set to the default value!\n"); }
    if(ps_fffloat(*par,"o3",&o3)) {o3=o3_def; fprintf(stdout, "WARNING: o3 is set to the default value!\n"); }
    
    if(ps_fffloat(*par,"d1",&d1)) {d1=d1_def; fprintf(stdout, "WARNING: d1 is set to the default value!\n"); }
    if(ps_fffloat(*par,"d2",&d2)) {d2=d2_def; fprintf(stdout, "WARNING: d2 is set to the default value!\n"); }
    if(ps_fffloat(*par,"d3",&d3)) {d3=d3_def; fprintf(stdout, "WARNING: d3 is set to the default value!\n"); }
    
    
    if(ps_ffint(*par,"n1",&n1)) {n1=n1_def; fprintf(stdout, "WARNING: n1 is set to the default value!\n"); }
    if(ps_ffint(*par,"n2",&n2)) {n2=n2_def; fprintf(stdout, "WARNING: n2 is set to the default value!\n"); }
    if(ps_ffint(*par,"n3",&n3)) {n3=n3_def; fprintf(stdout, "WARNING: n3 is set to the default value!\n"); }
    
    e1 = o1 + d1 * n1;
    e2 = o2 + d2 * n2;
    e3 = o3 + d3 * n3;
    
    if(ps_fffloat(*par,"CMAX",&cmax)) {cmax=cmax_def; fprintf(stdout,"WARNING: cmax is set to default value!\n");}
    if(ps_fffloat(*par,"CMIN",&cmin)) {cmin=cmin_def; fprintf(stdout,"WARNING: cmin is set to default value!\n");}
    if(ps_fffloat(*par,"mdl_rho",&rho)) {rho=rho_def; fprintf(stdout,"WARNING: rho is set to default value!\n");}

    if (model==1) {
      omv = (float*)emalloc(2*sizeof(float));
      if(ps_fffloat(*par,"top",omv  )) {omv[0]=o1; fprintf(stdout,"WARNING: top is set to o1!\n");}
      if(ps_fffloat(*par,"bot",omv+1)) {omv[1]=e1; fprintf(stdout,"WARNING: bot is set to e1!\n");}
    }
    if (model==2) {
        omv = (float*)emalloc(4*sizeof(float));
        if(ps_fffloat(*par,"cx1",omv  )) {
	  omv[0]=(e1-o1)*.5; 
	  fprintf(stdout, "WARNING: cx1 is set to the default value!\n"); 
	}
        if(ps_fffloat(*par,"cx2",omv+1)) {
	  omv[1]=(e2-o2)*.5; 
	  fprintf(stdout, "WARNING: cx2 is set to the default value!\n"); 
	}
        if(ps_fffloat(*par,"sigmax1",omv+2)) {
	  omv[2]=(e1-o1)*.1; 
	  fprintf(stdout, "WARNING: sigmax1 is set to the default value!\n"); 
	}
        if(ps_fffloat(*par,"sigmax2",omv+3)) {
	  omv[3]=(e2-o2)*.1; 
	  fprintf(stdout, "WARNING: sigmax2 is set to the default value!\n"); 
	}
    }
    if (model==3) {
        omv = (float*) emalloc(3*sizeof(float));
	if(ps_fffloat(*par,"top",omv  )) {omv[0]=o1; fprintf(stdout,"WARNING: top is set to o1!\n");}
	if(ps_fffloat(*par,"bot",omv+1)) {omv[1]=e1; fprintf(stdout,"WARNING: bot is set to e1!\n");}
        if(ps_fffloat(*par,"cycles",omv+2)) {
	  omv[2]=cycles_def; 
	  fprintf(stdout, "WARNING: cycles is set to the default value!\n"); 
	}
    }
    
    fprintf(stdout," o1=%f e1=%f d1=%f n1=%d\n",o1,e1,d1,n1);
    fprintf(stdout," o2=%f e2=%f d2=%f n2=%d\n",o2,e2,d2,n2);
    fprintf(stdout," o3=%f e3=%f d3=%f n3=%d\n",o3,e3,d3,n3);
  
    /* WWS */
    if (!(ps_ffcstring(*par,"hfile",&fname))) {
        
        /* DATAPATH deprecated - WWS 08.01.12*/
        /* DATAPATH revived - WWS 28.05.12 */
        /*    if (getenv("DATAPATH")) {
              dname=malloc(strlen(getenv("DATAPATH"))+strlen(fname)+2);
         //      strcpy(dname,getenv("DATAPATH"));
         // DATAPATH made preferentially a parameter WWS 08.28.12
         // else env variable else */
        cwdpath = NULL;
        /* if you can't get it from the parfile */
        if (ps_flcstring(*par,"datapath",&cwdpath)) {
            /* try to get it from the environment */
            cwdpath = (char *)malloc(128*sizeof(char));
            memset(cwdpath,'\0',128);
            
            pathptr = getenv("DATAPATH");
            if (pathptr) strcpy(cwdpath,pathptr);
            /* otherwise set to cwd */
            else strcpy(cwdpath,".");
            //strcpy(cwdpath,".");
        }
        dname=(char *)malloc(sizeof(char)*(strlen(cwdpath)+strlen(fname)+2));
        strcpy(dname,cwdpath);
        if (cwdpath[strlen(cwdpath)-1] != '/') strcat(dname,"/");
        strcat(dname,fname);
        strcat(dname,"@");
        
        fprintf(stdout,"writing header file %s\n",fname);
        if (!(fp=fopen(fname,"w"))) {
            fprintf(stderr,"Error: standardmodel\n");
            fprintf(stderr,"failed to open new header file %s\n",fname);
            exit(1);
        }
        fprintf(fp,"n1=%d d1=%e o1=%e\n",
                n1,d1,o1);
        fprintf(fp,"n2=%d d2=%e o2=%e\n",
                n2,d2,o2);
        fprintf(fp,"n3=%d d3=%e o3=%e\n",
                n3,d3,o3);
        switch (choose){
            case 0: fprintf(fp,"data_type=velocity\n"); break;
            case 1: fprintf(fp,"data_type=bulkmod\n"); break;
            case 2: fprintf(fp,"data_type=density\n"); break;
            case 3: fprintf(fp,"data_type=buoyancy\n"); break;
            default: fprintf(stderr,"choose failed\n"); exit(1);
        }
        
        fprintf(fp,"data_format=native_float\n");
        fprintf(fp,"in=%s\n",dname);
        fclose(fp);
        
        fprintf(stdout,"writing data to %s\n",dname);
        if (!(fp=fopen(dname,"w"))) {
            fprintf(stderr,"Error: standardmodel\n");
            fprintf(stderr,"failed to open new data file %s\n",dname);
            exit(1);
        }
    }
    
    else fp=stdout;
    
    ps_delete(&par);
    
    writemodel(choose,model,n1,n2,n3,o1,o2,o3,e1,e2,e3,d1,d2,d3,cmin,cmax,rho,fp,omv);
        
    if (model>0) free(omv);


    exit(0);
}
