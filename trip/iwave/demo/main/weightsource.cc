#include <par.h>
#include <parser.h>
#include <math.h>

/* define number of model types */
#define NMODEL 4

#define vp_max 4
#define vp_min 2
#define rho    2.3

#define o1_def 0
#define o2_def 0
#define o3_def 0
#define d1_def 5
#define d2_def 5
#define d3_def 5
#define n1_def 101
#define n2_def 101
#define n3_def 1

#define relsigmax1_def .2
#define relsigmax2_def .2

#define cycles_def 2

#ifndef PI
#define PI 3.14159265359
#endif

#define EMPTY -999

/*********************** self documentation **********************/
const char *sdoc[] = {
    "                                                                ",
    " WEIGHTSOURCE - weights traces of a source with model dependent ",
    "                values, either by den, vel, bulkmod, or buoy.   ",
    "                                                                ",
    " Units:		                                             ",
    " 	densities    [g/cm^3]                                        ",
    "   velocities   [m/ms] = [km/s]                                 ",
    " 	bulk moduli  [GPa] = [g/cm^3][m/ms]^2                        ",
    "   buoyancies   [cm^3/g]                                        ",
    "                                                                ",
    " Example run:   weightsource file.par in.bin out.bin            ",
    "                                                                ",
    " Required parameters:                                           ",
    "            file.par = file containing all parameters needed to ",
    "                       generate mdl and source                  ",
    "              in.bin = binary file containing input source data ",
    "             out.bin = binary output file                       ",
    NULL};

//--------------------------------------------------------------------------//
// model 0: Homogeneous model
static inline float homogeneous(int choose)
{
    //float vel = (vp_max+vp_min)/2;
    float vel = (vp_max+vp_min)/2;
    
    if (choose==1) return vel*vel*rho; 
    if (choose==2) return rho;
    if (choose==3) return 1./rho;
    
    return vel;
}
//--------------------------------------------------------------------------//
// model 1: linear depth velocity model with constant density
static inline float lineardepth(float x1, int choose, float * modelpars)
{
    float o1 = modelpars[0];
    float e1 = modelpars[1];
    
    float vel = (e1-x1)/(e1-o1)*vp_min + (x1-o1)/(e1-o1)*vp_max;
    
    if (choose==1) return vel*vel*rho;
    if (choose==2) return rho;
    if (choose==3) return 1./rho;
        
    return vel;
}
//--------------------------------------------------------------------------//
// model 2: negative gaussian lense in velocity, with constant density
static inline float gaussianlense2D(float x1, float x2, 
				    int choose, float * modelpars)
{
    float o1 = modelpars[0];
    float e1 = modelpars[1];
    float o2 = modelpars[2];
    float e2 = modelpars[3];
    float cx1 = modelpars[4];
    float cx2 = modelpars[5];
    float relsigmax1 = modelpars[6];
    float relsigmax2 = modelpars[7];
    
    float sigmax1 = relsigmax1*(e1-o1);
    float sigmax2 = relsigmax2*(e2-o2);
    
    float gaussx1 = exp( -(x1-cx1)*(x1-cx1)/(2.*sigmax1*sigmax1) );
    float gaussx2 = exp( -(x2-cx2)*(x2-cx2)/(2.*sigmax2*sigmax2) );
    
    float vel = vp_max - (vp_max-vp_min)*gaussx1*gaussx2;
    
    if (choose==1) return vel*vel*rho;
    if (choose==2) return rho;
    if (choose==3) return 1./rho;
    
    return vel;
}
//--------------------------------------------------------------------------//
// model 3: Sinusodial-depth velocity with constant density
static inline float Sindepth(float x1, int choose, float * modelpars)
{
   float o1 = modelpars[0];
    float e1 = modelpars[1];
    float cycles = modelpars[2];
    
    float wnum = cycles*2*PI;
    float vel = (vp_max-vp_min)*.5*sin(wnum*(x1-e1)/(e1-o1)) 
              + (vp_max+vp_min)*.5;
    
    if (choose==1) return vel*vel*rho;
    if (choose==2) return rho;
    if (choose==3) return 1./rho;
    
    return vel;
}

//--------------------------------------------------------------------------//
float get_zvalue(float x1, float x2,
                 int model, int choose, float * modelpars)
{
   float v;
	
    switch (model) {
        case 0: v = homogeneous(choose);                     break;
        case 1: v = lineardepth(x1,choose,modelpars);        break;
        case 2: v = gaussianlense2D(x1,x2,choose,modelpars); break;
        case 3: v = Sindepth(x1,choose,modelpars);           break;
        default: v = homogeneous(choose);
    }
    return v;
}

////////////////////////////////////////////////////////
// MAIN ////////////////////////////////////////////////
////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  PARARRAY *par;    
  char     *par_file;    //input parameter file name
  char     *src_bin;     //input binary file name
  char     *out_bin;
  FILE     *fp;
  int       err;

  int       horder;      //half order of FD method
  float     CMIN,CMAX;   //min./max. vel      [m/ms]
  float     fpeak;       //peak frequency     [Hz]
  
  int       nts;         //# of time steps, for src
  float     xcent;       //center of source   [m]
  float     zcent;       //center of source   [m]
  float     tcent;       //center of source   [s]
  float     lx,lz;       //diam. of source    [m]
  int       nx,nz;       //# of spatial steps
  float     src_dx;      //spatial discretization for generating source
  float     src_dz;       
  float     src_dt;      

  int       model;
  float     o1,o2;
  float     e1,e2;
  int       n1,n2;      
  float    *omv;
  float     cx1,cx2;
  float     relsigmax1,relsigmax2;
  float     cycles;
  
  float     GPC;         //Gridcells per cycle
  float     wavelength;  //min. wavelength    [m/cycle]
  float     mdl_dx;      //spatial discretization for generating model
  float     mdl_dz;       
  

  if (argc!=4)
  {
    fprintf(stderr,"Proper use:\nweightsource file.par in.bin out.bin\n");
    exit(1);
  }

  par_file = argv[1];
  src_bin  = argv[2];
  out_bin  = argv[3];
  
  //Extract parameters -------------------------------//
  par = ps_new();
  if (ps_createfile(par,par_file))
  {
    fprintf(stderr,"ERROR. Could not process par_file.\n");
    exit(1);
  }

  //Parsing ------------------------------------------//
  //find horder
  if ((ps_ffint(*par,"horder",&horder))) 
  {
    fprintf(stderr,"Error: could not find horder.\n");
    exit(1);
  }
  //find CMIN and CMAX
  if ((ps_fffloat(*par,"CMIN",&CMIN))) 
  {
    fprintf(stderr,"Error: could not find CMIN.\n");
    exit(1);
  }
  if ((ps_fffloat(*par,"CMAX",&CMAX))) 
  {
    fprintf(stderr,"Error: could not find CMAX.\n");
    exit(1);
  }
  //find fpeak
  if ((ps_fffloat(*par,"fpeak",&fpeak))) 
  {
    fprintf(stderr,"Error: could not find fpeak.\n");
    exit(1);
  }
  //find GPC or mdl_dx,mdl_dz
  if ((ps_fffloat(*par,"GPC",&GPC))) {
    if ((ps_fffloat(*par,"mdl_dx",&mdl_dx))&&
	(ps_fffloat(*par,"mdl_dz",&mdl_dz))) {
      fprintf(stderr,"Error: could not find either "
                     "GPC or mdl_dx,mdl_dz.\n");
      exit(1);
    } else GPC = (float)EMPTY;
  }
  else {
    mdl_dx = (float)EMPTY;
    mdl_dz = (float)EMPTY;
  }
  //find nts and src_dt
  if ((ps_ffint(*par,"nts",&nts))) {
    fprintf(stderr,"Error: could not find nts.\n");
    exit(1);
  }
  if ((ps_fffloat(*par,"src_dt",&src_dt))) {
    fprintf(stderr,"Error: could not find src_dt.\n");
    exit(1);
  }
  //find xcent, ycent, tcent
  if ((ps_fffloat(*par,"xcent",&xcent))) {
    fprintf(stderr,"Error: could not find xcent.\n");
    exit(1);
  }
  if ((ps_fffloat(*par,"zcent",&zcent))) {
    fprintf(stderr,"Error: could not find zcent.\n");
    exit(1);
  }
  if ((ps_fffloat(*par,"tcent",&tcent))) {
    fprintf(stderr,"Error: could not find tcent.\n");
    exit(1);
  }
  //find lx, lz
  if ((ps_fffloat(*par,"lx",&lx))) {
    fprintf(stderr,"Error: could not find lx.\n");
    exit(1);
  }
  if ((ps_fffloat(*par,"lz",&lz))) {
    fprintf(stderr,"Error: could not find lz.\n");
    exit(1);
  }
  //find src_dx, src_dz
  if ((ps_fffloat(*par,"src_dx",&src_dx))) {
    fprintf(stderr,"Error: could not find src_dx.\n");
    exit(1);
  }
  if ((ps_fffloat(*par,"src_dz",&src_dz))) {
    fprintf(stderr,"Error: could not find src_dz.\n");
    exit(1);
  }
  //find model 
  if ((ps_ffint(*par,"model",&model))) {
    fprintf(stderr,"Error: could not find model.\n");
    exit(1);
  }
  //find o1,o2
  if ((ps_fffloat(*par,"o1",&o1))) {
    fprintf(stderr,"Error: could not find o1.\n");
    exit(1);
  }
  if ((ps_fffloat(*par,"o2",&o2))) {
    fprintf(stderr,"Error: could not find o2.\n");
    exit(1);
  }
  //find n1,n2 or e1,e2
  if ((ps_ffint(*par,"n1",&n1))) {
    if((ps_fffloat(*par,"e1",&e1))){
      fprintf(stderr,"Error: could not find n1 or e1.\n");
      exit(1);
    } else n1 = EMPTY;
  } else e1 = (float)EMPTY;
  if ((ps_ffint(*par,"n2",&n2))) {
    if((ps_fffloat(*par,"e2",&e2))){
      fprintf(stderr,"Error: could not find n2 or e2.\n");
      exit(1);
    } else n2 = EMPTY;
  } else e2 = (float)EMPTY;

  //Computing some values ----------------------------//
  
  nx = 1+(int)(lx/src_dx);
  nz = 1+(int)(lz/src_dz);

  //compute d1 d2
  wavelength = CMIN/(fpeak/1000.0);
  if (mdl_dx==EMPTY) mdl_dx = wavelength/GPC*(horder/2.0);
  if (mdl_dz==EMPTY) mdl_dz = wavelength/GPC*(horder/2.0);
  float d1=mdl_dz, d2=mdl_dx;
  
  //compute e1,e2 or n1,n2
  if (n1==EMPTY) n1 = (int)((e1-o1)/d1)+1;
  else e1 = o1 + (n1-1)*d1;
  if (n2==EMPTY) n2 = (int)((e2-o2)/d2)+1;
  else e2 = o2 + (n2-1)*d2;  

  //other model variables
  if (model==2) {
    omv = (float*)emalloc(4*sizeof(float));
    
    if(ps_fffloat(*par,"cx1",omv  )) {
	omv[0]=(e1-o1)*.5; 
	fprintf(stdout,"WARNING: "
	               "cx1 is set to the default value!\n"); 
    } else omv[0] = cx1;
    
    if(ps_fffloat(*par,"cx2",omv+1)) {
	omv[1]=(e2-o2)*.5; 
	fprintf(stdout,"WARNING: "
	               "cx2 is set to the default value!\n"); 
    } else omv[1] = cx2;

    if(ps_fffloat(*par,"relsigmax1",omv+2)) {
	omv[2]=relsigmax1_def; 
	fprintf(stdout,"WARNING: "
	               "relsigmax1 is set to the default value!\n"); 
    } else omv[2] = relsigmax1;

    if(ps_fffloat(*par,"relsigmax2",omv+3)) {
	omv[3]=relsigmax2_def; 
	fprintf(stdout,"WARNING: "
	               "relsigmax2 is set to the default value!\n"); 
    } else omv[3] = relsigmax2;
  }
  if (model==3) {
    omv = (float*) emalloc(1*sizeof(float));
    if(ps_fffloat(*par,"cycles",omv)) {
	omv[0]=cycles_def; 
	fprintf(stdout,"WARNING: "
	               "cycles is set to the default value!\n"); 
    } else omv[0] = cycles;
  }

  // computing x,z coor of traces --------------------//
  //Assumption: x-coordinate is the fast moving coordinate.
  size_t nxnz = nx*nz;
  float *xcoor = (float*)calloc(nxnz,sizeof(float));
  float *zcoor = (float*)calloc(nxnz,sizeof(float));
  
  for(int iz; iz<nz; iz++){
    for(int ix; ix<nx; ix++){
      xcoor[iz*nx+ix] = xcent-lx/2.0 + ix*src_dx;
      zcoor[iz*nx+ix] = zcent-lz/2.0 + iz*src_dz;
    }
  }
  // computing bulkmod at trace coor -----------------//
  float *bulkmod_val = (float*)calloc(nxnz,sizeof(float));
  int choose = 1;
  
  for(int i=0; i<nxnz; i++)
    bulkmod_val[i] = get_zvalue( zcoor[i],xcoor[i],model,choose,omv );

  // read in traces from binary ----------------------//
  fp = fopen(src_bin,"rb");
  if(fp==NULL) {
    fprintf(stderr,"Error, from weightsource: "
		    "could not open file %s\n",src_bin);
    exit(1);
  }

  //Allocating room for data
  float *buff = (float*)calloc(nxnz*nts,sizeof(float)); 
  if(buff==NULL){
    fprintf(stderr,"Error, from weightsource: "
		   "could not allocate memory for buffer.\n");
  }  

  int final_nts = (int)(2.0/(fpeak*src_dt)) + nts;
  //reading in
  err = fread(buff,sizeof(float),nxnz*final_nts,fp); 
  if (err!=nxnz*final_nts){
    fprintf(stderr,"Error, from weightsource: "
	           "reading in binary data.\n");
    exit(1);
  }
  fclose(fp);

  // weight traces -----------------------------------//
  for(int i=0; i<nxnz; i++){
    for(int it=0; it<final_nts; it++){
      buff[it+i*final_nts] *= bulkmod_val[i];
    }
  }
  
  // write out traces --------------------------------//
  fp = fopen(out_bin,"wb");
  fwrite(buff,sizeof(float),nxnz*final_nts,fp);
  fclose(fp);

  free(omv);
  free(xcoor);
  free(zcoor);
  free(bulkmod_val);
  free(buff);
}
