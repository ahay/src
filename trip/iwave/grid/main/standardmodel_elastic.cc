/* build 3D velocity cube for SEAM standard models 

 * Original by John E. Anderson, after Matlab implementation by
 * Joakim Blanch.
 * RSF output by William W. Symes, Jan 08.
 * switch to SEAM getpar functions WWS Jan 08
 * units fixed to m, ms, and kg -  WWS Jan 08 
 * get_zvalue_1d(), writemodel_1d() by Dong Sun, March 09
 * camembert() by Dong Sun, Sep 10
 * elastic model generator by Xin Wang, Dec 10
 */

#include <par.h>
#include <parser.h>

/* define data path - suitable for use in demo dirs */
#define DATAPATH "./"

/* define number of model types */
#define NMODEL 11

#define vp_top			2.3
#define vs_top                  0.93897
#define density_top		2100
#define lambda_top	        7406
#define mu_top                  1851.5

#define vp_bottom		3.0
#define vs_bottom               1.2247
#define density_bottom	        2300
#define lambda_bottom		13800
#define mu_bottom               3450

#define vp_wedge		2.6
#define vs_wedge                1.3
#define density_wedge	        2200
#define lambda_wedge		7436
#define mu_wedge                3718

#define vp_wedge2		2.1
#define vs_wedge2               0.85734
#define density_wedge2	        2050.0
#define lambda_wedge2		6027 
#define mu_wedge2               1506.8

#define vp_h2o			1.5
#define vs_h2o                  0.0
#define density_h2o		1000.0
#define lambda_h2o		2250
#define mu_h2o                  0.0

float Vdefault[6] = {2100.0, 2.3, 0.93897, 7406, 1851.5, 1.0/2100};
float Vbgdefault[6] ={4000.0, 2.5, 1.25, 12500.0, 6250, 1.0/4000}; 
float default_relpert = 0.01;
float default_radius = 250.0;

/*********************** self documentation **********************/
const char *sdoc[] = {
  "                                                                ",
  " STANDARDMODEL_elastic - build 3D velocity cube for SEAM standard models",
  "                                                                ",
  " This module builds some standard velocity and density cubes    ",
  " useful for comparing analytic solutions to finite difference   ",
  " variable-density elastic simulators.                           ",
  "                                                                ", 
  " Units: m, kg, ms                                               ",
  " densities: kg/m^3.                                             ",
  " velocities: m/ms = km/s.                                       ",
  " lambda, mu: MPa                                                ",
  " buoyancies: m^3/kg                                             ",
  "                                                                ",
  " The output format is native binary floats.                     ",
  "                                                                ",
  " Example run:   standardmodel model=4 choose=2 > vp.bin         ",
  "          or:   standardmodel model=4 choose=2 hfile=vp.rsf     ",
  "                                                                ",
  " Optional parameters:                                           ",
  "   model=1                    = choice of standard model        ",
  "                              =1; layer over half-space         ",
  "                              =2; dipping layer over half-space ",
  "                                  (15 deg dip in x, 0 deg in y) ",
  "                              =3; step structure (100 m high)   ",
  "                              =4; thin layer one sample thick   ",
  "                              =5; wedge                         ",
  "                              =6; wedge2                        ",
  "                              =7; ABS/PML test (homogeneous)    ",
  "                              =8; dome                          ",
  "                              =9; SEAM dipping layer model      ",
  "                                  (25 deg dip in x, 15 deg in y)",
  "                              =10; layered model                ",
  "                              =11; Camembert model from Gauthier,",
  "                                  Virieux, and Tarantola, Geophysics",
  "                                  vol. 51, 1986                 ",
  "                                                                ",
  "   choose=1                   =1; output density                ",
  "                              =2; output vp                     ",
  "                              =3; output vs                     ",
  "                              =4; output lambda                 ",
  "                              =5; output mu                     ",
  "                              =6; output buoyancy               ",
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
  "   f1=0.                      starting coordinate on grid       ",
  "   d1=5.0                     increment on grid in meters       ",
  "   e1=1800.0                  end coordinate on grid            ",
  "   n1=(int)(1.5+(e1-f1)/d1)   number of gridpoints              ",
  "                                                                ",
  " Middle Dimension (x)                                           ",
  "   f2=0.                      starting coordinate on grid       ",
  "   d2=5.0                     increment on grid in meters       ",
  "   e2=6600.0                  end coordinate on grid            ",
  "   n2=(int)(1.5+(e2-f2)/d2)   number of gridpoints              ",
  "                                                                ",
  " Slow Dimension (y)                                             ",
  "   f3=0.                      starting coordinate on grid       ",
  "   d3=5.0                     increment on grid in meters       ",
  "   e3=6600.0                  end coordinate on grid            ",
  "   n3=(int)(1.5+(e3-f3)/d3)   number of gridpoints              ",
  "                                                                ",
  " Layered medium parameters (model=10):                          ",
  "   each parameter prefaced by m<choose>_, to indicate parameter ",
  "   type: eg. m3_numl indicates number of layers for bulk modulus",
  "     numl: number of layers                                     ",
  "     val<n>: value of parameter in layer n, n=1,...,numl        ",
  "     rf<n>: depth in m of bottom of layer n, n=1,...,numl-1     ",
  "   example: two density layers, bottom of first at 1000 m:      ",
  "     m1_numl=2 m1_val1=1000 m1_val2=2000 m1_rf1=1000            ",
  "                                                                ",
  NULL};

/*
  ####################################################################
  # draft implementation version 1 
  # based on SEAM Numerical Verification Procedure Proposal by J. Blanch
  ####################################################################
  # credits: J.E. Anderson  2007

  # This builds as an SU getpar module with make similar to that
  # for modules in:  $CWPSRC/src/par/main

  # compile and link:  
  cmtc standardmodel

  ####################################################################
  # Below are some SU style commands that might be useful for testing
  # this code and looking at example output.
  ####################################################################
  # Set choice for model in range of 1 to 7
  MODEL=5

  ####################################################################
  # example run to create density grid for model on 5 x 5 x 5 m grid
  ####################################################################
  ./standardmodel model=${MODEL} choose=1 d1=5. d2=5. d3=5. >model_density.bin

  ####################################################################
  # display just first panel of density grid for model 
  # (others are same for models 1-7)
  ####################################################################
  ximage n1=361 n2=1321 legend=1 wclip=1.0 bclip=3.0 cmap=hue \
  f1=0 d1=5 label1="depth (meters)" n1tic=5 \
  f2=0 d2=5 label2="x (meters)" n2tic=5 \
  wbox=1321 hbox=361 \
  title="model ${MODEL} density" <model_density.bin &

  ####################################################################
  # example run to create vp grid for model 1 on 5 x 5 x 5 m grid
  ####################################################################
  ./standardmodel model=${MODEL} choose=2 d1=5. d2=5. d3=5. >model_vp.bin

  ####################################################################
  # display just first panel of vp grid for model 
  # (others are same for models 1-7)
  ####################################################################
  ximage n1=361 n2=1321 legend=1 wclip=1500.0 bclip=4000.0 cmap=hue \
  f1=0 d1=5 label1="depth (meters)" n1tic=5 \
  f2=0 d2=5 label2="x (meters)" n2tic=5 \
  wbox=1321 hbox=361 \
  title="model ${MODEL} vp" <model_vp.bin &

  ####################################################################
  # example run to create dome model on 5 x 5 x 5 m grid
  ####################################################################
  ./standardmodel model=8 choose=2 d1=5. d2=5. d3=5. >model_vp.bin

  ####################################################################
  # display dome model with SU xmovie tool
  ####################################################################
  xmovie n1=361 n2=1321 n3=1321 \
  f1=0 d1=5 label1="depth (meters)" nTic1=5 \
  f2=0 d2=5 label2="x (meters)" nTic2=5 \
  wclip=1400. bclip=4000.0 width=1321 height=361 \
  fframe=0 dframe=5 title="y=%g meters" <model_vp.bin &


*/
/******************************************************************************/
/* model 1:      layer over half space */
static inline float half_space(float x1c, float x1, float valuea, float valueb)
{
  if (x1 < x1c) return valuea;
  else return valueb;
}
/******************************************************************************/
/* model 2:            15 degree dip over half space*/
static inline float dip(float x1c, float x2c, float x1, float x2,
			float valuea, float valueb)
{
	
  if ( x1 < x1c + (x2 - x2c) * tan(PI*15.0/180.) ) return valuea;
  else return valueb;

}
/******************************************************************************/
/* model 3:                100 m step */
static inline float step(float x1c, float x2c, float x1, float x2, 
			 float valuea, float valueb)
{
  if ( (x2 <  x2c) && (x1 < x1c + 50.0)  ) return valuea;
  if ( (x2 >= x2c) && (x1 < x1c - 100.0) ) return valuea;
  else return valueb;
}
/******************************************************************************/
/* model 4:         thin layer of water 1 sample thick */
static inline float thin_layer(int n1, int j1, float valuea, float valueb)
{
  if (j1 == n1/2) return valuea;
  else return valueb;
}
/******************************************************************************/
/* model 5:         wedge */
static inline float wedge(float x1c, float x2c, float x1, float x2, 
			  float valuea, float valueb, float valuec)
{
	
  if ( (x1 < x1c) && (x1 < x1c + (x2-x2c) * tan(PI*15.0/180.)) ) 
    return valuea;
  else if ( (x1 < x1c) && (x1 >= x1c + (x2-x2c) * tan(PI*15.0/180.)) ) 
    return valuec;
  else return valueb;
}
/******************************************************************************/
/* model 7:			homogeneous model for ABS test */
static inline float homogeneous(float valuea)
{	
  return valuea;
}
/******************************************************************************/
/* model 8:			homogeneous model for ABS test */
static inline float dome(float x1c, float x2c, float x3c, float x1, float x2, float x3, int choose)
{
  float value1, value2, value3, value4, value5, value6, value7;
  float rr, curve;
	
  rr  = (x2 - x2c) * (x2 - x2c) + (x3 - x3c) * (x3 - x3c) + (x1 - x1c) * (x1 - x1c);
  rr *= 0.000001;
  rr  = -rr * rr;
	
  curve = 1800.0-1200.0 * exp(rr);
  
  /* vp */  
  value1=1500.0*1.e-3;
  value2=3000.0*1.e-3;
  value3=1800.0*1.e-3;
  value4=2100.0*1.e-3;
  value5=3500.0*1.e-3;
  value6=2600.0*1.e-3;
  value7=4000.0*1.e-3;

  /* density */
  if (choose == 1){
    value1=1.0*1.e+3;
    value2=2.1*1.e+3;
    value3=1.8*1.e+3;
    value4=2.0*1.e+3;
    value5=2.2*1.e+3;
    value6=2.05*1.e+3;
    value7=2.3*1.e+3;
  }
  /* vs */
  else if (choose == 3) {
    value1=0.0;
    value2=1500.0*1.e-3;
    value3=900.0*1.e-3;
    value4=1050.0*1.e-3;
    value5=1750.0*1.e-3;
    value6=1300.0*1.e-3;
    value7=2000.0*1.e-3;
  }
  /* lambda */
  else if (choose == 4){
    value1 = 2250.0;
    value2 = 9450.0;
    value3 = 2916.0;
    value4 = 4410.0;
    value5 = 13475.0;
    value6 = 6929.0;
    value7 = 18400.0;
  }
  /* mu */
  else if (choose == 5) {
    value1 = 0.0;
    value2 = 4725.0;
    value3 = 1458.0;
    value4 = 2205.0;
    value5 = 6737.5;
    value6 = 3464.5;
    value7 = 9200.0;
  }
  
  if (x1 < 350.0) return value1;
  else if ( (x1 >=  350.0) && (x1 < 1000.0) && (x1 <= curve) ) return value2;
  else if ( (x1 >= 1000.0) && (x1 < 1400.0) && (x1 <= curve) ) return value5;
  else if ( (x1 >= 1400.0) && (x1 <= 1800.0) && (x1 <= curve)) return value7;
  else if ( (x1 >   curve) && (x1 <  800.0)				   ) return value3;
  else if ( (x1 >   curve) && (x1 >= 800.0) && (x1 < 1200.0) ) return value4;
  else if ( (x1 >   curve) && (x1 >= 1200.0)				   ) return value6;
  /*	
    if (x1 < curve) return 1.0;
    else return 2.0;
  */	
  return 0.0;
}
/******************************************************************************/
/* model 9:			dipping layer model */
static inline float seam_dip(int choose, float x1, float x2, float x3)
{
  float valuea, valueb;
  if (choose==1) {
    valuea=2100;
    valueb=2300;
  }
  else {
    valuea=2.3;
    valueb=3.0;
  }
  
  /*  if  ( ((x2 - 7500.0) *  sin(PI*25.0/180.0) + (x3 - 7500.0) * cos(PI*25.0/180.)) * tan(PI*15.0/180.)  - x1 + 7500.0 > 0) {
      return valuea; */
  if  ( ((x2 - 7500.0) *  sin(PI*25.0/180.0) + (x3 - 7500.0) * cos(PI*25.0/180.)) * tan(PI*15.0/180.)  - x1 + 7000.0 > 0) {
    return valuea;
  }
  else {
    return valueb;
  }
}

/******************************************************************************/
/* model 11:			Camembert Model */
static inline float camembert(float x1, float x2, float x3, float o1, float o2, float o3, float * values)
{
  float dis1, dis2, dis3, radius, bgval, rpert, valpert;
  dis1 = x1 - o1; 
  dis2 = x2 - o2;
  dis3 = x3 - o3;
  radius = values[0];
  bgval = values[1];
  rpert = values[2];
  valpert = bgval * rpert;
  if (dis1*dis1 + dis2*dis2 + dis3*dis3 - radius*radius > 0)
    return bgval;
  else 
    return bgval + valpert;
}

/******************************************************************************/
float get_zvalue(float x1c, float x2c, float x3c, float x1, float x2, float x3, 
		 int j1, int j2, int j3, int n1, int n2, int n3, 
		 int model, int choose, float valuea, float valueb, float valuec)
{	
  float v;
	
  switch (model){
  case 1: v = half_space(x1c, x1, valuea, valueb); break;
  case 2: v = dip(x1c, x2c, x1, x2, valuea, valueb); break;
  case 3: v = step(x1c, x2c, x1, x2, valuea, valueb); break;
  case 4: v = thin_layer(n1, j1, valuea, valueb); break;
  case 5: v = wedge(x1c, x2c, x1, x2, valuea, valueb, valuec); break;
  case 6: v = wedge(x1c, x2c, x1, x2, valuea, valueb, valuec); break;
  case 7: v = homogeneous(valuea); break;
  case 8: v = dome(x1c, x2c, x3c, x1, x2, x3, choose); break;
  case 9: v = seam_dip(choose,x1, x2, x3); break;
  default: v = half_space(x1c, x1, valuea, valueb); 
  }
  return v;
}
/******************************************************************************/
/* get values from layered model with 'numl' layers */
float get_zvalue_1d(float * x1rs,       /* array storing relector depths */
		    float * values,     /* array storing layered model values */
		    float x1,           /* current depth */ 
		    int numl)           /* number of layers = number of reflectors + 1 */
{	
  float v;
  int i;
  v=values[numl-1];
  for(i=0;i<numl-1;i++){
    if (x1< x1rs[i]) {
      v=values[i]; 
      break;
    }
  }
  return v;
}
/******************** end self doc **************/

int writemodel(
	       int choose,     /* =1 output density; =2 output vp =3 output vs  =4 output lambda =5 output mu*/
	       int model,      /* choice of model                     */
	       /*  =1; layer over half-space          */
	       /*  =2; dipping layer over half-space  */
	       /*  =3; step structure                 */
	       /*  =4; thin layer                     */
	       /*  =5; wedge                          */
	       /*  =6; wedge2                         */
	       /*  =7; ABS/PML test                   */
	       /*  =8; dome                           */
	       /*  =10; layered model                 */
	       /*  =11; camembert model               */
	       int n1,         /* number of samples in depth (fast dimension) */
	       int n2,         /* number of samples in x (middle dimension)   */
	       int n3,         /* number of samples in y (slow dimenssion)    */
	       float f1,       /* start depth */
	       float f2,       /* start x */
	       float f3,       /* start y */
	       float e1,       /* end depth */
	       float e2,       /* end x */
	       float e3,       /* end y */
	       float d1,       /* depth increment */
	       float d2,       /* x increment */
	       float d3,       /* y increment */
	       FILE * fp,       /* output stream */
	       int numl,        /* number of layers*/
	       float *x1rs,     /* reflector depth array */
	       float *values    /* model value array */
	       )
{
  int j1, j2, j3;
  float valuea, valueb, valuec;
  float x1, x2, x3, x1c, x2c, x3c;
  float *v   = NULL;
             
  valuea = vp_top;
  valueb = vp_bottom;
  valuec = vp_wedge;
  if (model == 6 ) valuec = vp_wedge2;

  if (choose == 1){
    valuea = density_top;
    valueb = density_bottom;
    valuec = density_wedge;
    if (model == 6 ) valuec = density_wedge2;
  }
  else if (choose == 3){
    valuea = vs_top;
    valueb = vs_bottom;
    valuec = vs_wedge;
    if (model == 6 ) valuec = vs_wedge2;
  }
  else if (choose == 4){
    valuea = lambda_top;
    valueb = lambda_bottom;
    valuec = lambda_wedge;
    if (model == 6 ) valuec = lambda_wedge2;
  }
  else if (choose == 5){
    valuea = mu_top;
    valueb = mu_bottom;
    valuec = mu_wedge;
    if (model == 6 ) valuec = mu_wedge2;
  }
  
  v   = (float *)emalloc(n1*sizeof(float));

  x1c = 0.5 * (f1 + e1);
  x2c = 0.5 * (f2 + e2);
  x3c = 0.5 * (f3 + e3);

  for(j3 = 0; j3 < n3; j3++){
    x3 = f3 + j3 * d3;
    for(j2 = 0; j2 < n2; j2++){
      x2 = f2 + j2 * d2;
      for(j1 = 0; j1 < n1; j1++){
	x1 = f1 + j1 * d1;
	if (model==10)
	  v[j1] = get_zvalue_1d(x1rs,values,x1,numl);            
	else if (model==11) 
	  v[j1] = camembert(x1,x2,x3,x1c,x2c,x3c,values);
	else
	  v[j1] = get_zvalue(x1c, x2c, x3c, x1, x2, x3, j1, j2, j3, 
			     n1, n2, n3, model, choose, valuea, valueb, valuec);
      }
      if (fwrite(v, sizeof(float), n1, fp) != (unsigned)n1){
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

  float f1,f2,f3,d1,d2,d3,e1,e2,e3;
  int n1,n2,n3,model,choose;
  
  /* WWS */
  char * fname;
  char * dname;
  FILE * fp;
  /* end WWS */
  
  /* DS */
  int numl[4];                               /* layer number array */
  float * x1rs[4] = {NULL,NULL,NULL,NULL};   /* pointers to reflector depth arrays */
  float * values[4] = {NULL,NULL,NULL,NULL}; /* pointers to layered value arrays in model 10*/
  
  unsigned int i,j;
  char numlname[10];
  char valname[10];
  char x1rsname[10];

  /* end DS */   
  
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
    if (model<1 || model > NMODEL) {
      fprintf(stderr,"Error: standardmodel.x\n");
      fprintf(stderr,"model index %d not defined: must lie in range [1,%d]\n",model,NMODEL);
      exit(1);
    }
  }
  else {
    fprintf(stderr,"Warning: standardmodel.x\n");
    fprintf(stderr,"no model index given, so using model=1\n");
    model=1;
  }
  /*   if(ps_getparint("choose",&choose)) { */
  if (!(ps_ffint(*par,"choose",&choose))) {
    if (choose<1 || choose > 5) {
      fprintf(stderr,"Error: standardmodel.x\n");
      fprintf(stderr,"choose index must be 1 (density),  2 (vp), or 3 (vs), or 4 (lambda), or 5 (mu)\n");
      exit(1);
    }
  }
  else {
    fprintf(stderr,"Warning: standardmodel.x\n");
    fprintf(stderr,"no choose index given, so using choose=1 (density)\n");
    choose=1;
  }
  
  /*
    if(!ps_getparfloat("f1",&f1)) {f1=0.; fprintf(stderr, "WARNING: f1 is set to the default value!\n"); }
    if(!ps_getparfloat("f2",&f2)) {f2=0.; fprintf(stderr, "WARNING: f2 is set to the default value!\n"); }
    if(!ps_getparfloat("f3",&f3)) {f3=0.; fprintf(stderr, "WARNING: f3 is set to the default value!\n"); }
    
    if(!ps_getparfloat("d1",&d1)) {d1=5.0; fprintf(stderr, "WARNING: d1 is set to the default value!\n"); }
    if(!ps_getparfloat("d2",&d2)) {d2=5.0; fprintf(stderr, "WARNING: d2 is set to the default value!\n"); }
    if(!ps_getparfloat("d3",&d3)) {d3=5.0; fprintf(stderr, "WARNING: d3 is set to the default value!\n"); }


    if(!ps_getparint("n1",&n1)) {n1 = 361; fprintf(stderr, "WARNING: n1 is set to the default value!\n"); }
    if(!ps_getparint("n2",&n2)) {n2 = 1321; fprintf(stderr, "WARNING: n2 is set to the default value!\n"); }
    if(!ps_getparint("n3",&n3)) {n3 = 1321; fprintf(stderr, "WARNING: n3 is set to the default value!\n"); }
  */

  if(ps_fffloat(*par,"f1",&f1)) {f1=0.; fprintf(stderr, "WARNING: f1 is set to the default value!\n"); }
  if(ps_fffloat(*par,"f2",&f2)) {f2=0.; fprintf(stderr, "WARNING: f2 is set to the default value!\n"); }
  if(ps_fffloat(*par,"f3",&f3)) {f3=0.; fprintf(stderr, "WARNING: f3 is set to the default value!\n"); }

  if(ps_fffloat(*par,"d1",&d1)) {d1=5.0; fprintf(stderr, "WARNING: d1 is set to the default value!\n"); }
  if(ps_fffloat(*par,"d2",&d2)) {d2=5.0; fprintf(stderr, "WARNING: d2 is set to the default value!\n"); }
  if(ps_fffloat(*par,"d3",&d3)) {d3=5.0; fprintf(stderr, "WARNING: d3 is set to the default value!\n"); }


  if(ps_ffint(*par,"n1",&n1)) {n1 = 361; fprintf(stderr, "WARNING: n1 is set to the default value!\n"); }
  if(ps_ffint(*par,"n2",&n2)) {n2 = 1321; fprintf(stderr, "WARNING: n2 is set to the default value!\n"); }
  if(ps_ffint(*par,"n3",&n3)) {n3 = 1321; fprintf(stderr, "WARNING: n3 is set to the default value!\n"); }

  e1 = f1 + d1 * (n1 - 1);
  e2 = f2 + d2 * (n2 - 1);
  e3 = f3 + d3 * (n3 - 1);
   
  fprintf(stderr," f1=%f e1=%f d1=%f n1=%d\n",f1,e1,d1,n1);
  fprintf(stderr," f2=%f e2=%f d2=%f n2=%d\n",f2,e2,d2,n2);
  fprintf(stderr," f3=%f e3=%f d3=%f n3=%d\n",f3,e3,d3,n3);


  /* DS */
  if (model == 10){
    for(i = 0; i< 4 ; i++){
      snprintf(numlname,10,"m%d_numl",i+1);
      /* set up layer numbers for velocity, density, bulkmod and buoyancy */
      if (ps_ffint(*par,numlname,numl+i)) numl[i] = 1;
       
      fprintf(stderr,"%s, numl[%d]=%d \n",numlname, i, numl[i]);

      /* allocate merroy for reflector depths and model values */
      if(numl[i] > 1)
	x1rs[i] = (float *) malloc((numl[i]-1)*sizeof(float));	   
      values[i] = (float *) malloc(numl[i]*sizeof(float));
       
      for (j =0; j< numl[i]; j++){
	if (j < numl[i]-1){
	  snprintf(x1rsname,10,"m%d_rf%d",i+1,j+1);
	  if(ps_fffloat(*par,x1rsname,x1rs[i]+j)) 
	    { x1rs[i][j] = e1; 
	      fprintf(stderr, "WARNING: x1rs[%d][%d] is set to the default value %f!\n",i,j,e1);
	    } 
	  else 
	    fprintf(stderr,"%s=%e\n",x1rsname,*(x1rs[i]+j));

	}
	snprintf(valname,10,"m%d_val%d",i+1,j+1);
	if(ps_fffloat(*par,valname,values[i]+j)) 
	  { values[i][j] = Vdefault[i]; 
	    fprintf(stderr, "WARNING: values[%d][%d] is set to the default value %f!\n",i,j,Vdefault[i]);
	  }
	else 
	  fprintf(stderr,"%s=%e\n",valname,*(values[i]+j));

      }  
    } 
  }

  if (model == 11){
    values[choose-1] = (float *) malloc(3*sizeof(float));
       
    if(ps_fffloat(*par,"radius",values[choose-1])) { 
      values[choose -1][0] = default_radius; 
      fprintf(stderr, "WARNING: camembert_radius is set to the default value %f!\n",default_radius);
    }

    if(ps_fffloat(*par,"bgvalue",values[choose-1]+1)) { 
      values[choose -1][1] = Vbgdefault[choose-1]; 
      fprintf(stderr, "WARNING: camembert_bgvalue is set to the default value %f!\n",Vbgdefault[choose-1]);
    }
     
    if(ps_fffloat(*par,"relpert",values[choose-1]+2)) { 
      values[choose-1][2] = default_relpert; 
      fprintf(stderr, "WARNING: camembert_relpert is set to the default value %f!\n",default_relpert);
    } 
  }

  /* end DS*/

  /* WWS */
  if (!(ps_ffcstring(*par,"hfile",&fname))) {

    dname=(char *)malloc(sizeof(char)*(strlen(DATAPATH)+strlen(fname)+2));
    strcpy(dname,DATAPATH);
    strcat(dname,fname);
    strcat(dname,"@");

    fprintf(stderr,"writing header file %s\n",fname);
    if (!(fp=fopen(fname,"w"))) {
      fprintf(stderr,"Error: standardmodel\n");
      fprintf(stderr,"failed to open new header file %s\n",fname);
      exit(1);
    }
    fprintf(fp,"n1=%d d1=%e o1=%e\n",
	    n1,d1,f1);
    fprintf(fp,"n2=%d d2=%e o2=%e\n",
	    n2,d2,f2);
    fprintf(fp,"n3=%d d3=%e o3=%e\n",
	    n3,d3,f3);
    switch (choose){
    case 1: fprintf(fp,"data_type=density\n"); break;
    case 2: fprintf(fp,"data_type=p-velocity\n"); break;
    case 3: fprintf(fp,"data_type=s-velocity\n"); break;
    case 4: fprintf(fp,"data_type=lambda\n"); break;
    case 5: fprintf(fp,"data_type=mu\n"); break;
    default: fprintf(stderr,"choose failed\n"); exit(1);
    }

    fprintf(fp,"data_format=native_float\n");
    fprintf(fp,"in=%s\n",dname);
    fclose(fp);

    fprintf(stderr,"writing data to %s\n",dname);
    if (!(fp=fopen(dname,"w"))) {
      fprintf(stderr,"Error: standardmodel\n");
      fprintf(stderr,"failed to open new data file %s\n",fname);
      exit(1);
    }
  }
  else fp=stdout;
  /* end WWS */
     
  /* original: 
     writemodel(choose,model,n1,n2,n3,f1,f2,f3,d1,d2,d3);
  */
  /* WWS */
  writemodel(choose,model,n1,n2,n3,f1,f2,f3,e1,e2,e3,d1,d2,d3,fp,numl[choose-1],x1rs[choose-1],values[choose-1]);
  /* end WWS */
  for (i=1; i<4; i++){
    free(x1rs[i]);
    free(values[i]);
  }
  exit(0);
}
