#include <par.h>
#include <parser.h>
#include <math.h>

#define DEBUG

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


//----------------------------------------------------------------------------//
// model 0: Homogeneous model
static inline float homogeneous(int choose, float CMIN, float CMAX, float rho)
{
    float vel = (CMAX+CMIN)/2;
    
    if (choose==1) return vel*vel*rho; 
    if (choose==2) return rho;
    if (choose==3) return 1./rho;
    
    return vel;
}
//----------------------------------------------------------------------------//
// model 1: linear depth velocity model with constant density
static inline float lineardepth(float x1, int choose, 
				float CMIN, float CMAX, 
				float rho, float * modelpars)
{
    float top = modelpars[0];
    float bot = modelpars[1];
    
    float vel = (bot-x1)/(bot-top)*CMIN + (x1-top)/(bot-top)*CMAX;
    
    if (choose==1) return vel*vel*rho;
    if (choose==2) return rho;
    if (choose==3) return 1./rho;
        
    return vel;
}
//----------------------------------------------------------------------------//
// model 2: negative gaussian lense in velocity, with constant density
static inline float gaussianlense2D(float x1, float x2, int choose, 
				    float CMIN, float CMAX, float rho, float * modelpars)
{
    float cx1 = modelpars[0];
    float cx2 = modelpars[1];
    float sigmax1 = modelpars[2];
    float sigmax2 = modelpars[3];
    
    float gaussx1 = exp( -(x1-cx1)*(x1-cx1)/(2.*sigmax1*sigmax1) );
    float gaussx2 = exp( -(x2-cx2)*(x2-cx2)/(2.*sigmax2*sigmax2) );
    
    float vel = CMAX - (CMAX-CMIN)*gaussx1*gaussx2;
    
    if (choose==1) return vel*vel*rho;
    if (choose==2) return rho;
    if (choose==3) return 1./rho;
    
    return vel;
}
//----------------------------------------------------------------------------//
// model 3: Sinusodial-depth velocity with constant density
static inline float Sindepth(float x1, int choose, 
			     float CMIN, float CMAX, 
			     float rho, float * modelpars)
{
    float top = modelpars[0];
    float bot = modelpars[1];
    float cycles = modelpars[2];
    
    float wnum = cycles*2*PI;
    float vel = (CMAX-CMIN)*.5*sin(wnum*(x1-bot)/(bot-top)) + (CMAX+CMIN)*.5;
    
    if (choose==1) return vel*vel*rho;
    if (choose==2) return rho;
    if (choose==3) return 1./rho;
    
    return vel;
}




/******************************************************************************/
float get_zvalue(float x1, float x2, 
                 int model, int choose, 
		 float CMIN, float CMAX, float rho,
		 float * modelpars)
{
    float v;
	
    switch (model){
    case 0: v = homogeneous(choose,CMIN,CMAX,rho); break;
    case 1: v = lineardepth(x1,choose,CMIN,CMAX,rho,modelpars); break;
    case 2: v = gaussianlense2D(x1,x2,choose,CMIN,CMAX,rho,modelpars); break;
    case 3: v = Sindepth(x1,choose,CMIN,CMAX,rho,modelpars); break;
    default: v = homogeneous(choose,CMIN,CMAX,rho);
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
  float    *modelpars;
  float     top,bot;
  float     cx1,cx2;
  float     sigmax1,sigmax2;
  float     cycles;
  float     rho;
  
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

  #ifdef DEBUG
  fprintf(stderr,"After extracting in/out file names from command line\n");
  #endif

  //Extract parameters -------------------------------//
  par = ps_new();
  if (ps_createfile(par,par_file))
  {
    fprintf(stderr,"ERROR. Could not process par_file.\n");
    exit(1);
  }

  #ifdef DEBUG
  fprintf(stderr,"After extracting pars\n");
  #endif
  
  xargc=argc; xargv=argv;
  requestdoc(0);

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
    nts = EMPTY;
    fprintf(stderr,"Warning: could not find nts. Using other value.\n");
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
    fprintf(stderr,"Error: src_dx not found!.\n");
    exit(1);
  }
  if ((ps_fffloat(*par,"src_dz",&src_dz))) {
    fprintf(stderr,"Error: src_dz not found!\n");
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

  //find auxiliary model parameters
  if ((ps_fffloat(*par,"cx1",&cx1)))
    cx1 = (float)EMPTY;
  if ((ps_fffloat(*par,"cx2",&cx2)))
    cx2 = (float)EMPTY;
  if ((ps_fffloat(*par,"sigmax1",&sigmax1)))
    sigmax1 = (float)EMPTY;
  if ((ps_fffloat(*par,"sigmax2",&sigmax2)))
    sigmax2 = (float)EMPTY;
  if ((ps_fffloat(*par,"cycles",&cycles)))
    cycles = (float)EMPTY;
  if ((ps_fffloat(*par,"mdl_top",&top)))
    top = (float)EMPTY;
  if ((ps_fffloat(*par,"mdl_bot",&bot)))
    bot = (float)EMPTY;
  if ((ps_fffloat(*par,"mdl_rho",&rho)))
    rho = (float)EMPTY;

  #ifdef DEBUG
  fprintf(stderr,"After finding parameters\n");
  #endif

  //Computing some values ----------------------------//
  
  nx = 1+(int)(lx/src_dx);
  nz = 1+(int)(lz/src_dz);

  //compute d1 d2
  wavelength = CMIN/(fpeak/1000.0);
  if (mdl_dx==EMPTY) mdl_dx = wavelength/GPC*(horder/2.0);
  if (mdl_dz==EMPTY) mdl_dz = wavelength/GPC*(horder/2.0);
  //NOTE: d1<0 since it refers to depth
  float d1=-mdl_dz, d2=mdl_dx;
  
  //compute e1,e2 or n1,n2
  //NOTE: converting e1 into a negative quantity
  if (e1!=EMPTY && e1>0) e1*=-1;

  if (n1==EMPTY) n1 = (int)((e1-o1)/d1);
  else e1 = o1 + (float)(n1)*d1;
  if (n2==EMPTY) n2 = (int)((e2-o2)/d2);
  else e2 = o2 + (float)(n2)*d2;  

  if (n1<0) n1*=-1;
  if (n2<0) n2*=-1;

  if (nts==EMPTY){
    nts=(int)(2.0/(fpeak*src_dt))+1; 
    fprintf(stderr,"Warning: nts set to %d.\n",nts);
  }

  // Other model variables -----------------------//
    if(ps_fffloat(*par,"CMAX",&CMAX)) {CMAX=cmax_def; fprintf(stdout,"WARNING: cmax is set to default value!\n");}
    if(ps_fffloat(*par,"CMIN",&CMIN)) {CMIN=cmin_def; fprintf(stdout,"WARNING: CMIN is set to default value!\n");}
    if(ps_fffloat(*par,"mdl_rho",&rho)) {rho=rho_def; fprintf(stdout,"WARNING: rho is set to default value!\n");}

    if (model==1) {
      modelpars = (float*)emalloc(2*sizeof(float));
      if(ps_fffloat(*par,"top",modelpars  )) {modelpars[0]=o1; fprintf(stdout,"WARNING: top is set to o1!\n");}
      if(ps_fffloat(*par,"bot",modelpars+1)) {modelpars[1]=e1; fprintf(stdout,"WARNING: bot is set to e1!\n");}
    }
    if (model==2) {
        modelpars = (float*)emalloc(4*sizeof(float));
        if(ps_fffloat(*par,"cx1",modelpars  )) {
	  modelpars[0]=(e1-o1)*.5; 
	  fprintf(stdout, "WARNING: cx1 is set to the default value!\n"); 
	}
        if(ps_fffloat(*par,"cx2",modelpars+1)) {
	  modelpars[1]=(e2-o2)*.5; 
	  fprintf(stdout, "WARNING: cx2 is set to the default value!\n"); 
	}
        if(ps_fffloat(*par,"sigmax1",modelpars+2)) {
	  modelpars[2]=(e1-o1)*.1; 
	  fprintf(stdout, "WARNING: sigmax1 is set to the default value!\n"); 
	}
        if(ps_fffloat(*par,"sigmax2",modelpars+3)) {
	  modelpars[3]=(e2-o2)*.1; 
	  fprintf(stdout, "WARNING: sigmax2 is set to the default value!\n"); 
	}
    }
    if (model==3) {
        modelpars = (float*) emalloc(3*sizeof(float));
	if(ps_fffloat(*par,"top",modelpars  )) {modelpars[0]=o1; fprintf(stdout,"WARNING: top is set to o1!\n");}
	if(ps_fffloat(*par,"bot",modelpars+1)) {modelpars[1]=e1; fprintf(stdout,"WARNING: bot is set to e1!\n");}
        if(ps_fffloat(*par,"cycles",modelpars+2)) {
	  modelpars[2]=cycles_def; 
	  fprintf(stdout, "WARNING: cycles is set to the default value!\n"); 
	}
    }
  
  #ifdef DEBUG
  fprintf(stderr,"After computing some parameters.\n");
  #endif

  // computing x,z coor of traces --------------------//
  //Assumption: x-coordinate is the fast moving coordinate.
  size_t nxnz = nx*nz;
  float *xcoor = (float*)calloc(nxnz,sizeof(float));
  float *zcoor = (float*)calloc(nxnz,sizeof(float));

  for(int iz=0; iz<nz; iz++){
    for(int ix=0; ix<nx; ix++){
      xcoor[iz*nx+ix] = xcent-lx/2.0 + ix*src_dx;
      zcoor[iz*nx+ix] = zcent-lz/2.0 + iz*src_dz;
    }
  }

  #ifdef DEBUG
  fprintf(stderr,"After computing x,z coor\n");
  #endif

  // computing bulkmod at trace coor -----------------//
  float *bulkmod_val = (float*)calloc(nxnz,sizeof(float));
  int choose = 1;
  
  for(int i=0; i<nxnz; i++)
    bulkmod_val[i] = get_zvalue( zcoor[i],xcoor[i],
				 model,choose,CMIN,CMAX,rho,modelpars );

  #ifdef DEBUG
  fprintf(stderr,"After computing bulkmod\n");
  #endif

  //Allocating room for data
  int final_nts = (int)(2.0/(fpeak*src_dt)) + nts + 1;
  fprintf(stderr,"Final nts = %d\n",final_nts);
  fprintf(stderr,"nxnz      = %d\n",nxnz);  

  float *buff = (float*)calloc(nxnz*final_nts,sizeof(float)); 
  if(buff==NULL){
    fprintf(stderr,"Error, from weightsource: "
		   "could not allocate memory for buffer.\n");
  }  

  #ifdef DEBUG
  fprintf(stderr,"After allocating mem for buff\n");
  #endif

  // read in traces from binary ----------------------//
  fprintf(stderr,"WEIGHTSOURCE: reading from %s.\n",src_bin);  
  fp = fopen(src_bin,"rb");
  if(fp==NULL) {
    fprintf(stderr,"Error, from weightsource: "
		    "could not open file %s\n",src_bin);
    exit(1);
  }
    
  err = fread(buff,sizeof(float),nxnz*final_nts,fp); 
  if (err!=nxnz*final_nts){
    fprintf(stderr,"Error, from weightsource: "
	           "reading in binary data.\n");
    exit(1);
  }
  fprintf(stderr,"After reading in binary.\n");

  fclose(fp);

  // weight traces -----------------------------------//
  for(int i=0; i<nxnz; i++){
    for(int it=0; it<final_nts; it++){
      buff[it+i*final_nts] *= bulkmod_val[i]*src_dx*src_dz;
      //buff[it*nxnz+i] *= bulkmod_val[i]*src_dx*src_dz;
    }
    //fprintf(stderr,"xcoor[%d]=%f, zcoor[%d]=%f\n",i,xcoor[i],i,zcoor[i]);
    //fprintf(stderr,">>>> bulkmod_val[%d]=%f\n",i,bulkmod_val[i]);
  }

  // write out traces --------------------------------//
  fprintf(stderr,"WEIGHTSOURCE: writing to %s.\n",out_bin);
  fp = fopen(out_bin,"wb");
  fwrite(buff,sizeof(float),nxnz*final_nts,fp);
  fclose(fp);

  fprintf(stderr,"After writing out binary.\n");

  if (model>=1) free(modelpars);
  free(xcoor);
  free(zcoor);
  free(bulkmod_val);
  free(buff);
}
