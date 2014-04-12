#include <par.h>
#include <parser.h>

#define EMPTY -999

/******************************************************************************/
/* This code generates .par files and python scripts required to build models */
/* and source/hrd files repesctively, which will be needed in SConstruct in   */
/* order to run IWAVE simulations.                                            */
/******************************************************************************/
int main(int argc, char **argv) {

  PARARRAY *par;    
  char     *par_file;    //input file name

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

  float     gx_o,gx_e;   //rec. start and end [m]
  float     gx_dx;       //dist. between rec. [m]
  float     gelev;       //rec. elevation     [m]
  int       ntr;         //# of traces, for hdr
  float     hdr_T;       //final time for hdr [s]
  float     hdr_dt;      //time step for hdr  [s]
  int       hdr_nt;      //# of time steps for hdr

  int       choose1;
  int       choose2;
  int       model;
  float     o1,o2;
  float     e1,e2;
  int       n1,n2;      

  float     CFL;         //CFL condition for FD code
  float     dt;          //time step for simulation
  float     GPC;         //Gridcells per cycle
  float     wavelength;  //min. wavelength    [m/cycle]
  float     mdl_dx;      //spatial discretization for generating model
  float     mdl_dz;       

  float     nl1, nr1, nl2, nr2;

  char     *mdl_hfile1;
  char     *mdl_pfile1;
  char     *mdl_path1;

  char     *mdl_hfile2;
  char     *mdl_pfile2;
  char     *mdl_path2;

  char     *src_file;
  char     *src_sc;

  char     *hdr_file;
  char     *hdr_sc;

  char     *data_file;
  char     *test_par;

  if (argc!=2){
    fprintf(stderr,"Proper use:\nparserParam par=name.par\n");
    exit(1);
  }
  
  //Extract par=name.par
  if ( *(argv[1])  !='p' ||
       *(argv[1]+1)!='a' ||
       *(argv[1]+2)!='r' ||
       *(argv[1]+3)!='=' ){
    fprintf(stderr,"Proper use:\nparserParam par=name.par\n");
    exit(1);
  }
   par_file = argv[1]+4;
  
  //Extract parameters -------------------------------//
  par = ps_new();
  if (ps_createfile(par,par_file)){
    fprintf(stderr,"ERROR. could not process command line\n");
    exit(1);
  }

  //Parsing ------------------------------------------//
  //find horder
  if ((ps_ffint(*par,"horder",&horder))) {
    fprintf(stderr,"Error: could not find horder.\n");
    exit(1);
  }
  //find CMIN and CMAX
  if ((ps_fffloat(*par,"CMIN",&CMIN))) {
    fprintf(stderr,"Error: could not find CMIN.\n");
    exit(1);
  }
  if ((ps_fffloat(*par,"CMAX",&CMAX))) {
    fprintf(stderr,"Error: could not find CMAX.\n");
    exit(1);
  }
  //find fpeak
  if ((ps_fffloat(*par,"fpeak",&fpeak))) {
    fprintf(stderr,"Error: could not find fpeak.\n");
    exit(1);
  }
  //find GPC or mdl_dx,mdl_dz
  if ((ps_fffloat(*par,"GPC",&GPC))) {
    if ((ps_fffloat(*par,"mdl_dx",&mdl_dx))&&
	(ps_fffloat(*par,"mdl_dz",&mdl_dz))) {
      fprintf(stderr,"Error: could not find either GPC or mdl_dx,mdl_dz.\n");
      exit(1);
    } else GPC = (float)EMPTY;
  }
  else {
    mdl_dx = (float)EMPTY;
    mdl_dz = (float)EMPTY;
  }
  //find CFL
  if ((ps_fffloat(*par,"CFL",&CFL))) {
    fprintf(stderr,"Error: could not find CFL.\n");
     exit(1);
  }   
  //find dt if given
  if ((ps_fffloat(*par,"dt",&dt)))
    dt = (float)EMPTY; 

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
  //find gx_o, gx_e, gx_dx
  if ((ps_fffloat(*par,"gx_o",&gx_o))) {
    fprintf(stderr,"Error: could not find gx_o.\n");
    exit(1);
  }
  if ((ps_fffloat(*par,"gx_e",&gx_e))) {
    fprintf(stderr,"Error: could not find gx_e.\n");
    exit(1);
  }
  if ((ps_fffloat(*par,"gx_dx",&gx_dx))) {
    if((ps_ffint(*par,"ntr",&ntr))) {
      fprintf(stderr,"Error: could not find gx_dx or ntr.\n");
      exit(1);
    } else gx_dx = (float)EMPTY;
  } else ntr = EMPTY;
  //find gelev
  if ((ps_fffloat(*par,"gelev",&gelev))) {
    fprintf(stderr,"Error: could not find gelev.\n");
    exit(1);
  }
  //find hdr_T,hdr_dt,hdr_nt
  if ((ps_fffloat(*par,"hdr_T",&hdr_T))) {
    fprintf(stderr,"Error: could not find hdr_T.\n");
    exit(1);
  }
  if ((ps_fffloat(*par,"hdr_dt",&hdr_dt))) {
    if ((ps_ffint(*par,"hdr_nt",&hdr_nt))) {
      fprintf(stderr,"Error: could not find hdr_dt or hdr_nt.\n");
      exit(1);
    } else hdr_dt = (float)EMPTY;
  } else hdr_nt = EMPTY;
  //find choose1, choose2, model
  if ((ps_ffint(*par,"choose1",&choose1))) {
    fprintf(stderr,"Error: could not find choose1.\n");
    exit(1);
  }  
  if ((ps_ffint(*par,"choose2",&choose2))) {
    fprintf(stderr,"Error: could not find choose2.\n");
    exit(1);
  } 
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

  ps_fffloat(*par,"nl1",&nl1);
  ps_fffloat(*par,"nr1",&nr1);
  ps_fffloat(*par,"nl2",&nl2);
  ps_fffloat(*par,"nr2",&nr2);

  //find mdl_hfile1, mdl_pfile1, mdl_path1
  if ((ps_ffcstring(*par,"mdl_hfile1",&mdl_hfile1))) {
    fprintf(stderr,"Error: could not find mdl_hfile1.\n");
    exit(1);
  }  
  if ((ps_ffcstring(*par,"mdl_pfile1",&mdl_pfile1))) {
    fprintf(stderr,"Error: could not find mdl_pfile1.\n");
    exit(1);
  }  
  if ((ps_ffcstring(*par,"mdl_path1",&mdl_path1))) {
    fprintf(stderr,"Error: could not find mdl_path1.\n");
    exit(1);
  }  
  //find mdl_hfile2, mdl_pfile2, mdl_path2
  if ((ps_ffcstring(*par,"mdl_hfile2",&mdl_hfile2))) {
    fprintf(stderr,"Error: could not find mdl_hfile2.\n");
    exit(1);
  }  
  if ((ps_ffcstring(*par,"mdl_pfile2",&mdl_pfile2))) {
    fprintf(stderr,"Error: could not find mdl_pfile2.\n");
    exit(1);
  }  
  if ((ps_ffcstring(*par,"mdl_path2",&mdl_path2))) {
    fprintf(stderr,"Error: could not find mdl_path2.\n");
    exit(1);
  }  
  //find hdr_file, hdr_sc, hdr_path
  if ((ps_ffcstring(*par,"hdr_file",&hdr_file))) {
    fprintf(stderr,"Error: could not find hdr_file.\n");
    exit(1);
  }  
  if ((ps_ffcstring(*par,"hdr_sc",&hdr_sc))) {
    fprintf(stderr,"Error: could not find hdr_sc.\n");
    exit(1);
  }  

  //find src_file, src_sc, src_path
  if ((ps_ffcstring(*par,"src_file",&src_file))) {
    fprintf(stderr,"Error: could not find src_file.\n");
    exit(1);
  }  
  if ((ps_ffcstring(*par,"src_sc",&src_sc))) {
    fprintf(stderr,"Error: could not find src_sc.\n");
    exit(1);
  }  

  //find data_file and test_par
  if ((ps_ffcstring(*par,"data_file",&data_file))) {
    fprintf(stderr,"Error: could not find data_file.\n");
    exit(1);
  }
  if ((ps_ffcstring(*par,"test_par",&test_par))) {
    fprintf(stderr,"Error: could not find test_par.\n");
    exit(1);
  }

  //Computing some values ----------------------------//
  wavelength = CMIN/(fpeak/1000.0);

  if (mdl_dx==EMPTY) mdl_dx = wavelength/GPC*(horder/2.0);
  if (mdl_dz==EMPTY) mdl_dz = wavelength/GPC*(horder/2.0);

  nx = 1+(int)(lx/src_dx);
  nz = 1+(int)(lz/src_dz);

  if (ntr==EMPTY) ntr = (int)((gx_e-gx_o)/gx_dx)+1;
  else gx_dx = (gx_e-gx_o)/((float)ntr);

  if (hdr_nt==EMPTY) hdr_nt = (int)(hdr_T/hdr_dt)+2;
  else hdr_dt = hdr_T/((float)hdr_nt);

  float d1=mdl_dz, d2=mdl_dx;
  if (n1==EMPTY) n1 = (int)((e1-o1)/d1)+1;
  else e1 = o1 + (float)(n1)*d1;
  if (n2==EMPTY) n2 = (int)((e2-o2)/d2)+1;
  else e2 = o2 + (float)(n2)*d2;

  //Writing out mdl par files ------------------------//
  fprintf(stderr,"Writing out mdl1 par files:\n"
	         "mdl_pfile1 = %s\n"
	         "mdl_hfile1 = %s\n"
   	         "mdl_path1  = %s\n",
	  mdl_pfile1,mdl_hfile1,mdl_path1);

  FILE * fp = fopen(mdl_pfile1,"w");
  fprintf(fp,"model=%d ",model);
  fprintf(fp,"choose=%d ",choose1);
  fprintf(fp,"o1=%f n1=%d d1=%f ",o1,n1,d1);
  fprintf(fp,"o2=%f n2=%d d2=%f ",o2,n2,d2);
  fprintf(fp,"hfile=%s ",mdl_hfile1);
  fclose(fp);

  fprintf(stderr,"Writing out mdl2 par files:\n"
	         "mdl_pfile2 = %s\n"
	         "mdl_hfile2 = %s\n"
	         "mdl_path2  = %s\n",
	  mdl_pfile2,mdl_hfile2,mdl_path2);
  
  fp = fopen(mdl_pfile2,"w");
  fprintf(fp,"model=%d ",model);
  fprintf(fp,"choose=%d ",choose2);
  fprintf(fp,"o1=%f n1=%d d1=%f ",o1,n1,d1);
  fprintf(fp,"o2=%f n2=%d d2=%f ",o2,n2,d2);
  fprintf(fp,"hfile=%s ",mdl_hfile2);
  fclose(fp);

  //Writing out script for hdr file ------------------//
  fprintf(stderr,"Writing out hdr script:\n"
	         "hdr_sc   = %s\n"
  	         "hdr_file = %s\n",
	  hdr_sc,hdr_file);
  fp = fopen(hdr_sc,"w");
  
  //importing commands
  fprintf(fp,"import os\n\n"
             "CWPROOT      = os.getenv('CWPROOT')\n"
	     "sunull       = os.path.join(CWPROOT,'bin/sunull')\n"
	     "sushw        = os.path.join(CWPROOT,'bin/sushw')\n\n");

  fprintf(fp,"cmd = sunull + ' nt=%d ntr=%d dt=%f | '",   hdr_nt,ntr,hdr_dt);
  fprintf(fp,    "+ sushw  + ' key=sx a=%f c=0 j=%d | '", xcent,ntr);
  fprintf(fp,    "+ sushw  + ' key=gx a=%f b=%f j=%d | '",gx_o,gx_dx,ntr);
  fprintf(fp,    "+ sushw  + ' key=delrt a=0 | '"
	         "+ sushw  + ' key=selev a=%f | '",       zcent);
  fprintf(fp,    "+ sushw  + ' key=gelev a=%f > %s'\n\n", gelev,hdr_file);

  //Execute cmd
  fprintf(fp,"os.system(cmd)\n");
  fclose(fp);

  //Writing out script for initial src file ----------//
  fprintf(stderr,"Writing out src script:\n"
	         "src_sc   = %s\n"
 	         "src_file = %s\n",
      	  src_sc,src_file);
  fp = fopen(src_sc,"w");

  //importing commands
  fprintf(fp,"import os\n\n"
             "CWPROOT      = os.getenv('CWPROOT')\n"
	     "sunull       = os.path.join(CWPROOT,'bin/sunull')\n"
	     "suconv       = os.path.join(CWPROOT,'bin/suconv')\n"
	     "suplane      = os.path.join(CWPROOT,'bin/suplane')\n"
	     "sutxtaper    = os.path.join(CWPROOT,'bin/sutxtaper')\n"
	     "sushw        = os.path.join(CWPROOT,'bin/sushw')\n"
	     "suwaveform   = os.path.join(CWPROOT,'bin/suwaveform')\n"
	     "sustrip      = os.path.join(CWPROOT,'bin/sustrip')\n"
	     "supaste      = os.path.join(CWPROOT,'bin/supaste')\n"
	     "IWAVE        = os.path.join(os.getenv('RSFSRC'),'trip/iwave')\n"
	     "towed_array  = os.path.join(IWAVE,'trace/main/towed_array.x')\n"
	     "weightsource = os.path.join(IWAVE,'demo/main/weightsource.x')\n"
	  "\n");
  
  //generate ricker
  fprintf(fp,"cmd1 = suwaveform+' type=ricker1 fpeak=%f dt=%f>ricker.su'\n\n",
	  fpeak,src_dt);

  //convolve with delta in time, and taper
  fprintf(fp,"cmd2 = "
	     "suplane +"
	     "' npl=1 nt=%d dt=%lf ntr=%d len1=%lf | '",
	  nts,src_dt,nz*nx,lz*lx);
  fprintf(fp,"+ sushw +" 
	     "' key=gelev a=%lf c=%lf b=0 j=%d | '",
	  -lz/2.0,src_dz,nz);
  fprintf(fp,"+ sushw +"
	     "' key=gx a=%lf b=%lf c=0 j=%d | '",
	  -lx/2.0,src_dx,nz);
  fprintf(fp,"+ suconv +"
	     "' sufile=ricker.su | '");
  fprintf(fp,"+ sutxtaper +"
	     "' taper=3 key=gelev min=0 max=0 dx=%lf | '",
	  lz/2.0);
  fprintf(fp,"+ sutxtaper +"
	     "' taper=3 key=gx min=0 max=0 dx=%lf > source_base.su'\n\n",
	  lx/2.0);

  //apply towed_array command
  fprintf(fp,"cmd3 = "
	     "towed_array+"
	     "' data=%s src=source_base.su towed=%s'\n\n",
	  hdr_file,src_file);

  //striping bin and hdr from su file 
  fprintf(fp,"cmd4 = "
	     "sustrip + "
	     "' head=temp.h ftn=0 <%s >temp.bin'\n\n",
	  src_file);

  //applying weights
  fprintf(fp,"cmd5 = "
   	     "weightsource + "
	     "' %s temp.bin new_temp.bin'\n\n",
	  par_file);

  //The number of time stps in the source file is initially nts,
  //but then gets altered after convolving with a ricker wavelet.
  int final_nts = (int)(2.0/(fpeak*src_dt)) +1 + nts;

  //pasting back bin 
  fprintf(fp,"cmd6 = "
	     "supaste + "
	     "' <new_temp.bin >%s ns=%d head=temp.h'\n\n",
	  src_file,final_nts);

  //executing commands
  fprintf(fp,"os.system(cmd1)\n"
	     "os.system(cmd2)\n"
	     "os.system(cmd3)\n"
	     "os.system(cmd4)\n"
	     "os.system(cmd5)\n"
	     "os.system(cmd6)\n\n");

  //cleanup
  fprintf(fp,"os.system('rm ricker.su; "
	                "rm source_base.su; "
                        "rm temp.bin; "
                        "rm new_temp.bin; "
	                "rm temp.h')\n");
  fclose(fp);

  //Writing par file to run simulation ---------------//
  fp = fopen(test_par,"w");

  fprintf(fp,"order = %d\n"
             "cfl   = %f\n"
             "cmin  = %f\n"
             "cmax  = %f\n"
             "dmin  = 0.5\n"
	     "dmax  = 5.0\n"
	     "fpeak = %f\n\n",
	  horder, CFL, CMIN-0.5, CMAX+0.5, fpeak/1000.0);

  if (dt!=EMPTY) 
    fprintf(fp,"dt = %f\n",dt);
  
  fprintf(fp,"nl1   = %f\n"
	     "nr1   = %f\n"
	     "nl2   = %f\n"
	     "nr2   = %f\n\n",
	  nl1, nr1, nl2, nr2);

  fprintf(fp,"srctype = array\n"
	     "sampord = 1\n"
	     "source  = ../%s\n\n",
	  src_file);

  fprintf(fp,"hdrfile  = ../%s\n"
      	     "datafile = %s\n\n",
	  hdr_file, data_file);
  
  if (choose1==0)
    fprintf(fp,"velocity = ");
  else if(choose1==1)
    fprintf(fp,"bulkmodulus = ");
  else if(choose1==2)
    fprintf(fp,"density = ");
  else
    fprintf(fp,"buoyancy = ");
  fprintf(fp,"../%s%s\n",mdl_path1,mdl_hfile1);

  if (choose2==0)
    fprintf(fp,"velocity = ");
  else if(choose2==1)
    fprintf(fp,"bulkmodulus = ");
  else if(choose2==2)
    fprintf(fp,"density = ");
  else
    fprintf(fp,"buoyancy = ");
  fprintf(fp,"../%s%s\n\n",mdl_path2,mdl_hfile2);

  fprintf(fp,"dump_lda  =1\n"
             "dump_ldc  =1\n"
      	     "dump_term =1\n\n");

  fprintf(fp,"data_file = %s",data_file);

  fclose(fp);
  ps_delete(&par);
}
