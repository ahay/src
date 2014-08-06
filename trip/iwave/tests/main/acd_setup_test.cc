#include "gtest/gtest.h"
#include "acd_defn.hh"
#include "grid.h"
#include "traceio.h"
#include "istate.hh"

//#define GTEST_VERBOSE

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",    0, true,  true },
  {"data",   1, false, true },
  {"source", 1, true,  false},
  {"movie",  1, false, false},
  {"init",   1, true,  false},
  {"",       0, false, false}
};

namespace {

  using RVL::parse;
  using RVL::RVLException;
  using TSOpt::IWaveEnvironment;
  using TSOpt::IWaveTree;
  using TSOpt::IWaveSampler;
  using TSOpt::IWaveSim;
  using TSOpt::TASK_RELN;
  using TSOpt::IOTask;

  void create_fixed_2D_data(string fnm, 
			    int nt,
			    float dt,
			    float ot,
			    int nrec,
			    int ntr_in_rec,
			    float sx0,
			    float dsx,
			    float sz,
			    float rx0,
			    float drx,
			    float rz,
			    int scalel,
			    int scalco) {
			   
    string CWPROOT = getenv("CWPROOT");
    stringstream cmd;
    cmd << CWPROOT + "/bin/sunull nt=" <<nt<<" ntr="<<ntr_in_rec*nrec<<" dt="<<dt;
    cmd << "| " + CWPROOT + "/bin/sushw key=delrt a="<<ot;
    cmd << "| " + CWPROOT + "/bin/sushw key=sx a="<<sx0<<" c="<<dsx<<" j="<<ntr_in_rec;
    cmd << "| " + CWPROOT + "/bin/sushw key=selev a=" << (-1.0f)*sz;
    cmd << "| " + CWPROOT + "/bin/sushw key=gx a="<<rx0<<" b="<<drx<<" j="<<ntr_in_rec;
    cmd <<" | " + CWPROOT + "/bin/sushw key=gelev a="<<(-1.0f)*rz;
    cmd << "| " + CWPROOT + "/bin/sushw key=scalel a="<<scalel;
    cmd << "| " + CWPROOT + "/bin/sushw key=scalco a="<<scalco;
    cmd << "| " + CWPROOT + "/bin/suchw key1=offset key2=gx key3=sx b=1 c=-1 > "<<fnm<<"\n";

    if (system(cmd.str().c_str())) {
      RVLException e;
      e<<"Error: create_fixed_2D_data\n";
      e<<"  failed system call on command \n";
      e<<cmd.str();
      throw e;
    }

  }

  void create_hfile(string hfile, string dfile, grid g, float val, bool var=false) {
    
    if (hfile.size()>120) {
      RVLException e;
      e<<"Error: create_hfile\n";
      e<<"filename "<<hfile<<" longer than 120 chars\n";
      throw e;
    }
    
    char * fname=(char *)malloc(128*sizeof(char));
    FILE * fp = NULL;
      
    strcpy(fname,hfile.c_str());
    
    fp = iwave_fopen(&fname,"w",NULL,stderr);
    if (!fp) {
      RVLException e;
      e<<"Error: create_hfile\n";
      e<<"file testgrid.rsf not opened\n";
      throw e;
    }
    
    fprint_grid(fp,g);
    
    fprintf(fp,"data_format=native_float\n");
    
    fprintf(fp,"data_type = csq\n");
    fprintf(fp,"in=%s\n",dfile.c_str());
    
    fflush(fp);
    
    iwave_fclose(fp);
    
    strcpy(fname,dfile.c_str());
    float * buf = (float *)malloc(get_datasize_grid(g)*sizeof(float));
    if (var) {
      for (int i=0;i<get_datasize_grid(g);i++) buf[i]=val/((float)(i+1));
    }
    else {
      for (int i=0;i<get_datasize_grid(g);i++) buf[i]=val;
    }
    //      FILE * 
    fp = iwave_fopen(&fname,"w",NULL,stderr);
    if (!fp) {
      RVLException e;
      e<<"Error: create_hfile\n";
      e<<"file "<<fname<<" not opened\n";
      throw e;
    }
    for (int i=0;i<get_panelnum_grid(g);i++) 
      fwrite(buf,sizeof(float),get_datasize_grid(g),fp);
    
    free(buf);
    fflush(fp);
    
    iwave_fclose(fp);    
    free(fname);
  }

  void create_gauss(grid g, float lam) {
    
    char * fname=(char *)malloc(128*sizeof(char));
    FILE * fp = NULL;

    strcpy(fname,"gauss.rsf");
      
    fp = iwave_fopen(&fname,"w",NULL,stderr);
    if (!fp) {
      RVLException e;
      e<<"Error: create_hfile\n";
      e<<"file gauss.rsf not opened\n";
      throw e;
    }
    
    fprint_grid(fp,g);
    
    fprintf(fp,"data_format=native_float\n");
    
    fprintf(fp,"data_type = csq\n");
    fprintf(fp,"in=%s\n","gauss.rsf@");
    
    fflush(fp);
    
    iwave_fclose(fp);
    
    strcpy(fname,"gauss.rsf@");
    
    // compute center
    RPNT ctr;
    float radsq=0.0f;
    
    for (int i=0;i<g.dim;i++) {
      ctr[i] = g.axes[i].o + 0.5 * g.axes[i].d * g.axes[i].n;
      radsq += (0.5 * g.axes[i].d * g.axes[i].n)*(0.5 * g.axes[i].d * g.axes[i].n);
    }
    float x;
    float y; 
    float z;
    float rsq;
    float * buf = (float *)malloc(get_datasize_grid(g)*sizeof(float));
    if (g.dim==2) {
      //	cerr<<"g.axes[0].n="<<g.axes[0].n<<" g.axes[1].n="<<g.axes[1].n<<endl;
      //	cerr<<"g.axes[0].o="<<g.axes[0].o<<" g.axes[1].o="<<g.axes[1].o<<endl;
      //	cerr<<"g.axes[0].d="<<g.axes[0].d<<" g.axes[1].d="<<g.axes[1].d<<endl;
      //	cerr<<"radsq="<<radsq<<endl;
      for (int j=0;j<g.axes[1].n;j++) {
	for (int i=0;i<g.axes[0].n;i++) {
	  z=g.axes[0].o + i*g.axes[0].d;
	  x=g.axes[1].o + j*g.axes[1].d;
	  rsq=(z-ctr[0])*(z-ctr[0])+(x-ctr[1])*(x-ctr[1]);
	  //	    cerr<<"rsq/radsq="<<rsq/radsq<<endl;
	  buf[i+j*g.axes[0].n]=exp(-lam*lam*rsq/radsq);
	}
      }
    }
    else if (g.dim==3) {
      for (int k=0;k<g.axes[2].n;k++) {
	for (int j=0;j<g.axes[1].n;j++) {
	  for (int i=0;i<g.axes[0].n;i++) {
	    z=g.axes[0].o + i*g.axes[0].d;
	    x=g.axes[1].o + j*g.axes[1].d;
	    y=g.axes[2].o + k*g.axes[2].d;
	    buf[i+j*g.axes[0].n]=exp(-lam*lam*((z-ctr[0])*(z-ctr[0])+(x-ctr[1])*(x-ctr[1])+(y-ctr[2])*(y-ctr[2])));
	  }
	}
      }  
    }
    else {
      RVLException e;
      e<<"Error: create_hfile\n";
      e<<"  wake up and smell the roses\n";
      throw e;
    }
    //      FILE * 
    fp = iwave_fopen(&fname,"w",NULL,stderr);
    if (!fp) {
      RVLException e;
      e<<"Error: create_hfile\n";
      e<<"file "<<fname<<" not opened\n";
      throw e;
    }
    for (int i=0;i<get_panelnum_grid(g);i++) 
      fwrite(buf,sizeof(float),get_datasize_grid(g),fp);
    
    free(buf);
    fflush(fp);
    
    iwave_fclose(fp);
    free(fname);
  }

  class ACDSimTest : public ::testing::Test {
  public:

    string hfile;
    string dfile;
    string hfile1;
    string dfile1;
    string hfile2;
    string dfile2;
    string hfile3;
    string dfile3;
    string hfiles;
    string dfiles;
    string hfilem;
    string dfilem;
    string hfilei;
    string dfilei;

    IWaveInfo ic;

    ACDSimTest(): ic() {

	int nt;
	float dt;
	float ot;
	int nrec;
	int ntr_in_rec;
	float sx0;
	float dsx;
	float sz;
	float rx0;
	float drx;
	float rz;
	int scalel; 
	int scalco;

	// first segy - single gather
	nt = 101;
	dt = 0.004;
	ot = 0.0;
	nrec = 3;
	ntr_in_rec = 101;
	sx0=3300;
	dsx=100;
	sz=40;
	rx0=100;
	drx=20;
	rz=50;
	scalel=0;
	scalco=0;
	
	string dn = "data.su";
	
	create_fixed_2D_data(dn,nt,dt,ot,
			     nrec,ntr_in_rec,
			     sx0,dsx,sz,
			     rx0,drx,rz,
			     scalel,scalco);

	// segcond segy - source pulses
	nt = 51;
	dt = 0.004;
	ot = -100.0;
	nrec = 3;
	ntr_in_rec = 1;
	sx0=3300;
	dsx=100;
	sz=40;
	rx0=0.0;
	drx=0.0;
	rz=20;
	scalel=0;
	scalco=0;
	
	string wn = "wavelet.su";
	
	create_fixed_2D_data(wn,nt,dt,ot,
			     nrec,ntr_in_rec,
			     sx0,dsx,sz,
			     rx0,drx,rz,
			     scalel,scalco);

      // hfile etc. = simple 2D
      // hfile1 etc. = x, z reversed
      // hfile2      = gdim=3, n3=3 etc.
      // hfile3      = gdim=3, h=1, z=2, x=3
      hfile = "csq_4layer.rsf";
      dfile = "csq_4layer.rsf@";
      hfile1= "csq1.rsf";
      dfile1= "csq1.rsf@";
      hfile2= "csq2.rsf";
      dfile2= "csq2.rsf@";
      hfile3= "csq3.rsf";
      dfile3= "csq3.rsf@";
      hfiles= "csqsmall.rsf";
      dfiles= "csqsmall.rsf@";
      hfilem= "movie.rsf";
      dfilem= "movie.rsf@";
      hfilei= "init.rsf";
      dfilei= "init.rsf@";

      // create simple par file
      ofstream ps("parfile");
      ps<<"INPUT DATA FOR IWAVE\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"FD:\n";
      ps<<"\n";
      ps<<"         order = 4           scheme half-order\n";
      ps<<"           cfl = 0.5        cfl number - frac of max stable\n";
      ps<<"          cmin = 1.0         min velocity - checked\n";
      ps<<"          cmax = 5.0         max velocity - checked\n";
      ps<<"      max_step = 0           1 = set adaptively, 0 = use standard cfl from cmax\n";
      ps<<"         fpeak = 0.010       nominal central frequency \n";
      ps<<"\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"Model info:\n";
      ps<<"\n";
      ps<<"           csq = csq_4layer.rsf\n";
      ps<<"        csq_d1 = csq_4layer.rsf\n";
      ps<<"        csq_d2 = csq_4layer.rsf\n";
      ps<<"        csq_b1 = csq_4layer.rsf\n";
      ps<<"        csq_b2 = csq_4layer.rsf\n";
      ps<<"\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"MPI info:\n";
      ps<<"\n";
      ps<<"       mpi_np1 = 1      n_doms along axis 1\n";
      ps<<"       mpi_np2 = 1      n_doms along axis 2\n";
      ps<<"       mpi_np3 = 1      n_doms along axis 3\n";
      ps<<"       partask = 1      task parallelization\n";
      ps<<"\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"Source info:\n";
      ps<<"\n";
      ps<<"        source = wavelet.su\n";
      ps<<"       sampord = 1             sampling order\n";
      ps<<"\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"Trace info:\n";
      ps<<"\n";
      ps<<"            data = data.su    output data file\n";
      ps<<"\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"Output info:\n";
      ps<<"\n";
      ps<<"     printact = 1           per-time-step verbosity level\n";
      ps<<"                            0 - none\n";
      ps<<"                            1 - time step index\n";
      ps<<"                            2 - internal time step info\n";
      ps<<"                            > 5: dump everything\n";
      ps<<"      dump_pi = 1           dump parallel/dom. decomp info\n";
      ps<<"     dump_lda = 1           dump grid data for allocated arrays\n";
      ps<<"     dump_ldc = 1           dump grid data for computational arrays\n";
      ps<<"     dump_lds = 1           dump grid data for send arrays\n";
      ps<<"     dump_ldr = 1           dump grid data for receive arrays\n";
      ps<<"    dump_term = 0           dump terminator data\n";
      ps<<"    dump_pars = 0           print parameter table in IWaveOp\n";
      ps<<"   dump_steps = 0           print major steps in IWaveOp\n";
      ps.flush();
      ps.close();

      ofstream qs("initfile");
      qs<<"INPUT DATA FOR IWAVE\n";
      qs<<"------------------------------------------------------------------------\n";
      qs<<"FD:\n";
      qs<<"\n";
      qs<<"         order = 2           scheme half-order\n";
      qs<<"           cfl = 0.5        cfl number - frac of max stable\n";
      qs<<"          cmin = 1.0         min velocity - checked\n";
      qs<<"          cmax = 2.0         max velocity - checked\n";
      qs<<"      max_step = 0           1 = set adaptively, 0 = use standard cfl from cmax\n";
      qs<<"         fpeak = 0.010       nominal central frequency \n";
      qs<<"            dt = 2.0";
      qs<<"\n";
      qs<<"------------------------------------------------------------------------\n";
      qs<<"Model info:\n";
      qs<<"\n";
      qs<<"           csq = csq_4layer.rsf\n";
      qs<<"          init = gauss.rsf\n";
      qs<<"         movie = movie.rsf\n";
      qs<<"\n";
      qs<<"------------------------------------------------------------------------\n";
      qs<<"MPI info:\n";
      qs<<"\n";
      qs<<"       mpi_np1 = 1      n_doms along axis 1\n";
      qs<<"       mpi_np2 = 1      n_doms along axis 2\n";
      qs<<"       mpi_np3 = 1      n_doms along axis 3\n";
      qs<<"       partask = 1      task parallelization\n";
      qs<<"\n";
      qs<<"------------------------------------------------------------------------\n";
      qs<<"Output info:\n";
      qs<<"\n";
      qs<<"     printact = 1           per-time-step verbosity level\n";
      qs<<"                            0 - none\n";
      qs<<"                            1 - time step index\n";
      qs<<"                            2 - internal time step info\n";
      qs<<"                            > 5: dump everything\n";
      qs<<"      dump_pi = 1           dump parallel/dom. decomp info\n";
      qs<<"     dump_lda = 1           dump grid data for allocated arrays\n";
      qs<<"     dump_ldc = 1           dump grid data for computational arrays\n";
      qs<<"     dump_lds = 1           dump grid data for send arrays\n";
      qs<<"     dump_ldr = 1           dump grid data for receive arrays\n";
      qs<<"    dump_term = 0           dump terminator data\n";
      qs<<"    dump_pars = 0           print parameter table in IWaveOp\n";
      qs<<"   dump_steqs = 0           print major steqs in IWaveOp\n";
      qs.flush();
      qs.close();

      grid g;
      g.dim=2;
      g.gdim=2;
      g.axes[0].n=416;
      g.axes[1].n=800;
      g.axes[2].n=3;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[2].d=100.0;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[2].o=-100.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 3;

      float val = 1.5;
      create_hfile(hfile, dfile, g, val);

      g.dim=2;
      g.gdim=2;
      g.axes[0].n=416;
      g.axes[1].n=800;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;

      val = 1.0;
      create_hfile(hfilei, dfilei, g, val);

      g.axes[0].n=800;
      g.axes[1].n=416;      
      g.axes[0].id = 1;
      g.axes[1].id = 0;

      val=1.5;
      create_hfile(hfile1, dfile1, g, val);

      g.gdim=3;
      g.axes[0].n=416;
      g.axes[1].n=800;      
      g.axes[2].n=3;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[2].d=1.0;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[2].o=0.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 3;
      create_hfile(hfile2, dfile2, g, val);

      g.axes[0].n=3;
      g.axes[0].d=100.0;
      g.axes[0].o=-100.0;
      g.axes[1].n=416;
      g.axes[1].d=25.0;
      g.axes[1].o=0.0;
      g.axes[2].n=800;
      g.axes[2].d=25.0;
      g.axes[2].o=0.0;
      g.axes[0].id = 101;
      g.axes[1].id = 0;
      g.axes[2].id = 1;
      create_hfile(hfile3, dfile3, g, val);

      g.dim=2;
      g.gdim=2;
      g.axes[0].n=5;
      g.axes[1].n=10;
      g.axes[2].n=1;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[2].d=1.0;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[2].o=0.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 2;
      create_hfile(hfiles, dfiles, g, 1.0, true);

      g.dim=2;
      g.gdim=3;
      g.axes[0].n=416;
      g.axes[1].n=800;
      g.axes[2].n=33;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[2].d=200.0;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[2].o=2.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 2;

      val = 0.0;
      create_hfile(hfilem, dfilem, g, val);

      g.dim=2;
      g.gdim=3;
      g.axes[0].n=416;
      g.axes[1].n=800;
      g.axes[2].n=1;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[2].d=2.0;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[2].o=0.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 2;

      float lam = 10;
      create_gauss(g,lam);

    }
  };

  TEST_F(ACDSimTest, setup_environment) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;
#ifdef GTEST_VERBOSE
      ps_printall(*par,stderr);
#endif

      // check one of the main par entries
      string csq="";
      parse(*par,"csq",csq);
      EXPECT_EQ("csq_4layer.rsf",csq);
      
      ps_delete(&par);
      fclose(stream);

#ifndef GTEST_VERBOSE
      //      unlink("parfile");
#endif
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, iwavetree_2D_serial_order0) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=0;
      IWaveTree * wt = new IWaveTree(*par,stream,ic,order);
      IWAVE & w = *(wt->getStateArray()[0]);
      iwave_printf(&w, par, stream);
      IPNT cgs0, cge0, gs0, ge0;
      ra_a_gse(&((w.model).ld_c._s[0]),gs0,ge0);
      ra_gse(&((w.model).ld_c._s[0]),cgs0,cge0);
      EXPECT_EQ(1,gs0[0]);
      EXPECT_EQ(414,ge0[0]);
      EXPECT_EQ(1,gs0[1]);
      EXPECT_EQ(798,ge0[1]);
      EXPECT_EQ(1,cgs0[0]);
      EXPECT_EQ(414,cge0[0]);
      EXPECT_EQ(1,cgs0[1]);
      EXPECT_EQ(798,cge0[1]);
      IPNT cgs1, cge1, gs1, ge1;
      ra_a_gse(&((w.model).ld_c._s[1]),gs1,ge1);
      ra_gse(&((w.model).ld_c._s[1]),cgs1,cge1);
      EXPECT_EQ(-3,gs1[0]);
      EXPECT_EQ(418,ge1[0]);
      EXPECT_EQ(-3,gs1[1]);
      EXPECT_EQ(802,ge1[1]);
      EXPECT_EQ(1,cgs1[0]);
      EXPECT_EQ(414,cge1[0]);
      EXPECT_EQ(1,cgs1[1]);
      EXPECT_EQ(798,cge1[1]);
      IPNT cgs2, cge2, gs2, ge2;
      ra_a_gse(&((w.model).ld_c._s[2]),gs2,ge2);
      ra_gse(&((w.model).ld_c._s[2]),cgs2,cge2);
      EXPECT_EQ(-3,gs2[0]);
      EXPECT_EQ(418,ge2[0]);
      EXPECT_EQ(-3,gs2[1]);
      EXPECT_EQ(802,ge2[1]);
      EXPECT_EQ(1,cgs2[0]);
      EXPECT_EQ(414,cge2[0]);
      EXPECT_EQ(1,cgs2[1]);
      EXPECT_EQ(798,cge2[1]);

      delete wt;
      ps_delete(&par);
      fclose(stream);
#ifndef GTEST_VERBOSE

#endif
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }    

  TEST_F(ACDSimTest, iwavetree_2D_serial_order0_axis_swap) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // modify - set csq=csq1.rsf
      ps_slcstring(*par,"csq","csq1.rsf");

      // build order zero IWaveTree, check 
      int order=0;
      IWaveTree * wt = new IWaveTree(*par,stream,ic,order);
      IWAVE & w = *(wt->getStateArray()[0]);
      iwave_printf(&w, par, stream);
      IPNT cgs0, cge0, gs0, ge0;
      ra_a_gse(&((w.model).ld_c._s[0]),gs0,ge0);
      ra_gse(&((w.model).ld_c._s[0]),cgs0,cge0);
      EXPECT_EQ(1,gs0[1]);
      EXPECT_EQ(414,ge0[1]);
      EXPECT_EQ(1,gs0[0]);
      EXPECT_EQ(798,ge0[0]);
      EXPECT_EQ(1,cgs0[1]);
      EXPECT_EQ(414,cge0[1]);
      EXPECT_EQ(1,cgs0[0]);
      EXPECT_EQ(798,cge0[0]);
      IPNT cgs1, cge1, gs1, ge1;
      ra_a_gse(&((w.model).ld_c._s[1]),gs1,ge1);
      ra_gse(&((w.model).ld_c._s[1]),cgs1,cge1);
      EXPECT_EQ(-3,gs1[1]);
      EXPECT_EQ(418,ge1[1]);
      EXPECT_EQ(-3,gs1[0]);
      EXPECT_EQ(802,ge1[0]);
      EXPECT_EQ(1,cgs1[1]);
      EXPECT_EQ(414,cge1[1]);
      EXPECT_EQ(1,cgs1[0]);
      EXPECT_EQ(798,cge1[0]);
      IPNT cgs2, cge2, gs2, ge2;
      ra_a_gse(&((w.model).ld_c._s[2]),gs2,ge2);
      ra_gse(&((w.model).ld_c._s[2]),cgs2,cge2);
      EXPECT_EQ(-3,gs2[1]);
      EXPECT_EQ(418,ge2[1]);
      EXPECT_EQ(-3,gs2[0]);
      EXPECT_EQ(802,ge2[0]);
      EXPECT_EQ(1,cgs2[1]);
      EXPECT_EQ(414,cge2[1]);
      EXPECT_EQ(1,cgs2[0]);
      EXPECT_EQ(798,cge2[0]);

      delete wt;
      ps_delete(&par);
      fclose(stream);
#ifndef GTEST_VERBOSE

#endif
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }    

  TEST_F(ACDSimTest, iwavetree_2D_serial_order0_ext_after) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // modify - set csq=csq2.rsf
      ps_slcstring(*par,"csq","csq2.rsf");

      // build order zero IWaveTree, check 
      int order=0;
      IWaveTree * wt = new IWaveTree(*par,stream,ic,order);
      IWAVE & w = *(wt->getStateArray()[0]);
      iwave_printf(&w, par, stream);
      IPNT cgs0, cge0, gs0, ge0;
      ra_a_gse(&((w.model).ld_c._s[0]),gs0,ge0);
      ra_gse(&((w.model).ld_c._s[0]),cgs0,cge0);
      EXPECT_EQ(1,gs0[0]);
      EXPECT_EQ(414,ge0[0]);
      EXPECT_EQ(1,gs0[1]);
      EXPECT_EQ(798,ge0[1]);
      EXPECT_EQ(1,cgs0[0]);
      EXPECT_EQ(414,cge0[0]);
      EXPECT_EQ(1,cgs0[1]);
      EXPECT_EQ(798,cge0[1]);
      IPNT cgs1, cge1, gs1, ge1;
      ra_a_gse(&((w.model).ld_c._s[1]),gs1,ge1);
      ra_gse(&((w.model).ld_c._s[1]),cgs1,cge1);
      EXPECT_EQ(-3,gs1[0]);
      EXPECT_EQ(418,ge1[0]);
      EXPECT_EQ(-3,gs1[1]);
      EXPECT_EQ(802,ge1[1]);
      EXPECT_EQ(1,cgs1[0]);
      EXPECT_EQ(414,cge1[0]);
      EXPECT_EQ(1,cgs1[1]);
      EXPECT_EQ(798,cge1[1]);
      IPNT cgs2, cge2, gs2, ge2;
      ra_a_gse(&((w.model).ld_c._s[2]),gs2,ge2);
      ra_gse(&((w.model).ld_c._s[2]),cgs2,cge2);
      EXPECT_EQ(-3,gs2[0]);
      EXPECT_EQ(418,ge2[0]);
      EXPECT_EQ(-3,gs2[1]);
      EXPECT_EQ(802,ge2[1]);
      EXPECT_EQ(1,cgs2[0]);
      EXPECT_EQ(414,cge2[0]);
      EXPECT_EQ(1,cgs2[1]);
      EXPECT_EQ(798,cge2[1]);

      delete wt;
      ps_delete(&par);
      fclose(stream);
#ifndef GTEST_VERBOSE

#endif
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }    
  TEST_F(ACDSimTest, iwavetree_2D_serial_order0_ext_before) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // modify - set csq=csq3.rsf
      ps_slcstring(*par,"csq","csq3.rsf");

      // build order zero IWaveTree, check 
      int order=0;
      IWaveTree * wt = new IWaveTree(*par,stream,ic,order);
      IWAVE & w = *(wt->getStateArray()[0]);
      iwave_printf(&w, par, stream);
      IPNT cgs0, cge0, gs0, ge0;
      ra_a_gse(&((w.model).ld_c._s[0]),gs0,ge0);
      ra_gse(&((w.model).ld_c._s[0]),cgs0,cge0);
      EXPECT_EQ(1,gs0[0]);
      EXPECT_EQ(414,ge0[0]);
      EXPECT_EQ(1,gs0[1]);
      EXPECT_EQ(798,ge0[1]);
      EXPECT_EQ(1,cgs0[0]);
      EXPECT_EQ(414,cge0[0]);
      EXPECT_EQ(1,cgs0[1]);
      EXPECT_EQ(798,cge0[1]);
      IPNT cgs1, cge1, gs1, ge1;
      ra_a_gse(&((w.model).ld_c._s[1]),gs1,ge1);
      ra_gse(&((w.model).ld_c._s[1]),cgs1,cge1);
      EXPECT_EQ(-3,gs1[0]);
      EXPECT_EQ(418,ge1[0]);
      EXPECT_EQ(-3,gs1[1]);
      EXPECT_EQ(802,ge1[1]);
      EXPECT_EQ(1,cgs1[0]);
      EXPECT_EQ(414,cge1[0]);
      EXPECT_EQ(1,cgs1[1]);
      EXPECT_EQ(798,cge1[1]);
      IPNT cgs2, cge2, gs2, ge2;
      ra_a_gse(&((w.model).ld_c._s[2]),gs2,ge2);
      ra_gse(&((w.model).ld_c._s[2]),cgs2,cge2);
      EXPECT_EQ(-3,gs2[0]);
      EXPECT_EQ(418,ge2[0]);
      EXPECT_EQ(-3,gs2[1]);
      EXPECT_EQ(802,ge2[1]);
      EXPECT_EQ(1,cgs2[0]);
      EXPECT_EQ(414,cge2[0]);
      EXPECT_EQ(1,cgs2[1]);
      EXPECT_EQ(798,cge2[1]);

      delete wt;
      ps_delete(&par);
      fclose(stream);
#ifndef GTEST_VERBOSE

#endif
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }    

  TEST_F(ACDSimTest, iwavetree_2D_serial_order1) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=1;
      IWaveTree * wt = new IWaveTree(*par,stream,ic,order);
      IWAVE & w = *(wt->getStateArray()[1]);
      iwave_printf(&w, par, stream);
      IPNT cgs0, cge0, gs0, ge0;
      ra_a_gse(&((w.model).ld_c._s[0]),gs0,ge0);
      ra_gse(&((w.model).ld_c._s[0]),cgs0,cge0);
      EXPECT_EQ(1,gs0[0]);
      EXPECT_EQ(414,ge0[0]);
      EXPECT_EQ(1,gs0[1]);
      EXPECT_EQ(798,ge0[1]);
      EXPECT_EQ(1,cgs0[0]);
      EXPECT_EQ(414,cge0[0]);
      EXPECT_EQ(1,cgs0[1]);
      EXPECT_EQ(798,cge0[1]);
      IPNT cgs1, cge1, gs1, ge1;
      ra_a_gse(&((w.model).ld_c._s[1]),gs1,ge1);
      ra_gse(&((w.model).ld_c._s[1]),cgs1,cge1);
      EXPECT_EQ(-3,gs1[0]);
      EXPECT_EQ(418,ge1[0]);
      EXPECT_EQ(-3,gs1[1]);
      EXPECT_EQ(802,ge1[1]);
      EXPECT_EQ(1,cgs1[0]);
      EXPECT_EQ(414,cge1[0]);
      EXPECT_EQ(1,cgs1[1]);
      EXPECT_EQ(798,cge1[1]);
      IPNT cgs2, cge2, gs2, ge2;
      ra_a_gse(&((w.model).ld_c._s[2]),gs2,ge2);
      ra_gse(&((w.model).ld_c._s[2]),cgs2,cge2);
      EXPECT_EQ(-3,gs2[0]);
      EXPECT_EQ(418,ge2[0]);
      EXPECT_EQ(-3,gs2[1]);
      EXPECT_EQ(802,ge2[1]);
      EXPECT_EQ(1,cgs2[0]);
      EXPECT_EQ(414,cge2[0]);
      EXPECT_EQ(1,cgs2[1]);
      EXPECT_EQ(798,cge2[1]);

      delete wt;
      ps_delete(&par);
      fclose(stream);
#ifndef GTEST_VERBOSE

#endif
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }    

  TEST_F(ACDSimTest, isample_construct_rsf_input) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      IWAVE w;
      if (int err = iwave_construct(&w,par,stream,ic)) {
	RVLException e;
	e<<"Error: ACDSimTest - isample_construct_rsf_input\n";
	e<<"  error from iwave_construct = "<<err<<"\n";
	throw e;
      }
      string key="csq";
      
      IWaveSampler s(&w,key,*par,stream);
      EXPECT_EQ(2,s.getNumAxes());
      EXPECT_EQ(416,s.getAxis(0).n);
      EXPECT_EQ(25.0,s.getAxis(0).d);
      EXPECT_EQ(0.0,s.getAxis(0).o);
      EXPECT_EQ(0,s.getAxis(0).id);
      EXPECT_EQ(800,s.getAxis(1).n);
      EXPECT_EQ(25.0,s.getAxis(1).d);
      EXPECT_EQ(0.0,s.getAxis(1).o);
      EXPECT_EQ(1,s.getAxis(1).id);

      iwave_destroy(&w,ic.get_mdest());
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, isample_construct_rsf_input_small) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // modify - set csq=csqsmall.rsf
      ps_slcstring(*par,"csq","csqsmall.rsf");

      IWAVE w;
      if (int err = iwave_construct(&w,par,stream,ic)) {
	RVLException e;
	e<<"Error: ACDSimTest - isample_construct_rsf_input\n";
	e<<"  error from iwave_construct = "<<err<<"\n";
	throw e;
      }
      string key="csq";
      
      IWaveSampler s(&w,key,*par,stream);
      EXPECT_EQ(2,s.getNumAxes());
      EXPECT_EQ(5,s.getAxis(0).n);
      EXPECT_EQ(25.0,s.getAxis(0).d);
      EXPECT_EQ(0.0,s.getAxis(0).o);
      EXPECT_EQ(0,s.getAxis(0).id);
      EXPECT_EQ(10,s.getAxis(1).n);
      EXPECT_EQ(25.0,s.getAxis(1).d);
      EXPECT_EQ(0.0,s.getAxis(1).o);
      EXPECT_EQ(1,s.getAxis(1).id);

      iwave_destroy(&w,ic.get_mdest());
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, isample_sample_rsf_fwd_input_small) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // modify - set csq=csqsmall.rsf
      ps_slcstring(*par,"csq","csqsmall.rsf");

      IWAVE w;
      if (int err = iwave_construct(&w,par,stream,ic)) {
	RVLException e;
	e<<"Error: ACDSimTest - isample_construct_rsf_input\n";
	e<<"  error from iwave_construct = "<<err<<"\n";
	throw e;
      }
      string key="csq";
      
      IWaveSampler s(&w,key,*par,stream);

      grid gsim;
      init_default_grid(&gsim);
      gsim.dim=2;
      gsim.gdim=4;
      gsim.axes[0].n=5;
      gsim.axes[1].n=10;
      gsim.axes[2].n=2;
      gsim.axes[3].n=3;
      gsim.axes[0].d=25.0;
      gsim.axes[1].d=25.0;
      gsim.axes[2].d=10.0;
      gsim.axes[3].d=100.0;
      gsim.axes[0].o=0.0;
      gsim.axes[1].o=0.0;
      gsim.axes[2].o=-10.0;
      gsim.axes[3].o=-100.0;
      gsim.axes[0].id=0;
      gsim.axes[1].id=1;
      gsim.axes[2].id=2;
      gsim.axes[3].id=3;
      
      IPNT step;
      IASN(step,IPNT_0);
      step[2]=-1; // should be istart
      step[3]=-1;

      // must zero first - sample acts in update mode!
      rd_a_zero(&(w.model.ld_a));
      s.sample(gsim,step,true,true,&w,0,0,stream);

      IPNT rags; 
      IPNT rage;
      rd_gse(&(w.model.ld_a),0,rags,rage);
      /*
      cerr<<"rags[0]="<<rags[0]<<" rage[0]="<<rage[0]<<"\n";
      cerr<<"rags[1]="<<rags[1]<<" rage[1]="<<rage[1]<<"\n";
      for (int i1=rags[1]; i1<=rage[1]; i1++) {
	for (int i0=rags[0]; i0<=rage[0]; i0++) {
	  cerr<<"f["<<i1<<"]["<<i0<<"]="<<w.model.ld_a._s[0]._s2[i1][i0]<<"\n";
	}
      }
      */

      EXPECT_EQ(rags[0],1);
      EXPECT_EQ(rags[1],1);
      EXPECT_EQ(rage[0],3);
      EXPECT_EQ(rage[1],8);
      for (int i1=rags[1]; i1<=rage[1]; i1++) {
	for (int i0=rags[0]; i0<=rage[0]; i0++) {
	  float x = fabs(1.0/(1.0+i1*gsim.axes[0].n+i0) - w.model.ld_a._s[0]._s2[i1][i0]);
	  float e = 1.e-06 * fabs(1.0/(1.0+i1*gsim.axes[0].n+i0));
	  EXPECT_LT(x,e);
	}
      }
      iwave_destroy(&w,ic.get_mdest());
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, isample_construct_su_input) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      IWAVE w;
      if (int err = iwave_construct(&w,par,stream,ic)) {
	RVLException e;
	e<<"Error: ACDSimTest - isample_construct_su_input\n";
	e<<"  error from iwave_construct = "<<err<<"\n";
	throw e;
      }
      string key="data";
      
      IWaveSampler s(&w,key,*par,stream);
      EXPECT_EQ(2,s.getNumAxes());
      EXPECT_EQ(227,s.getAxis(0).n);
      EXPECT_GT(1.e-6,fabs(s.getAxis(0).d-1.767767));
      EXPECT_EQ(0.0,s.getAxis(0).o);
      EXPECT_EQ(2,s.getAxis(0).id);
      EXPECT_EQ(3,s.getAxis(1).n);
      EXPECT_EQ(1.0,s.getAxis(1).d);
      EXPECT_EQ(0.0,s.getAxis(1).o);
      EXPECT_EQ(3,s.getAxis(1).id);

      iwave_destroy(&w,ic.get_mdest());
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, spacetime_grid_build) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;

      // initialize grid
      grid g;
      init_default_grid(&g);

      // iotask list
      std::vector<TASK_RELN *> t;
      // step 1: create list of i/o tasks
      IOTask(t,order,fwd,ic);

      IWaveTree * wt = new IWaveTree(*par,stream,ic,order);

      for (int it=0; it<t.size(); it++) {
	IWaveSampler * tmp = NULL;
	//	cerr<<"constructing sampler on "<<t[it]->keyword<<endl;
	tmp = new IWaveSampler(wt->getStateArray()[0], t[it]->keyword, *par, stream);
	// note these are internal simulation axes, not archival, 
	// which will require some re-arranging in traceio
	for (int i=0;i<tmp->getNumAxes();i++) {
	  //	  cerr<<"axis "<<i<<endl;
	  //	  fprint_axis(stderr,tmp->getAxis(i));
	  if (!grid_union(&g,&(tmp->getAxis(i)))) {
	    RVLException e;
	    e<<"Error: IWaveSim constructor from grid_union\n";
	    fprintf(stream,"Error: IWaveSim constructor from grid_union\n");
	    fprintf(stream,"  failed to add these axes:\n");
	    for (int j=0;j<tmp->getNumAxes(); j++) 
	      fprint_axis(stream,tmp->getAxis(j));
	    fprintf(stream,"  to grid:\n");
	    fprint_grid(stream,g);
	    throw e;
	  }
	}
	// for this test no need to build sampler vector
	//	s.push_back(tmp);
	delete tmp;
      }

      EXPECT_EQ(4,g.gdim);
      EXPECT_EQ(416,g.axes[0].n);
      EXPECT_EQ(25.0,g.axes[0].d);
      EXPECT_EQ(0.0,g.axes[0].o);
      EXPECT_EQ(0,g.axes[0].id);
      EXPECT_EQ(800,g.axes[1].n);
      EXPECT_EQ(25.0,g.axes[1].d);
      EXPECT_EQ(0.0,g.axes[1].o);
      EXPECT_EQ(1,g.axes[1].id);
      EXPECT_EQ(284,g.axes[2].n);
      EXPECT_GT(1.e-6,fabs(g.axes[2].d-1.767767));
      EXPECT_GT(1.e-4,fabs(-100.7627-g.axes[2].o));
      EXPECT_EQ(2,g.axes[2].id);
      EXPECT_EQ(3,g.axes[3].n);
      EXPECT_EQ(1.0,g.axes[3].d);
      EXPECT_EQ(0.0,g.axes[3].o);
      EXPECT_EQ(3,g.axes[3].id);

      delete wt;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, spacetime_grid_build_axis_swap) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // modify - set csq=csq1.rsf
      ps_slcstring(*par,"csq","csq1.rsf");

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;

      // initialize grid
      grid g;
      init_default_grid(&g);

      // iotask list
      std::vector<TASK_RELN *> t;
      // step 1: create list of i/o tasks
      IOTask(t,order,fwd,ic);

      IWaveTree * wt = new IWaveTree(*par,stream,ic,order);

      for (int it=0; it<t.size(); it++) {
	IWaveSampler * tmp = NULL;
	//	cerr<<"constructing sampler on "<<t[it]->keyword<<endl;
	tmp = new IWaveSampler(wt->getStateArray()[0], t[it]->keyword, *par, stream);
	// note these are internal simulation axes, not archival, 
	// which will require some re-arranging in traceio
	for (int i=0;i<tmp->getNumAxes();i++) {
	  //	  cerr<<"axis "<<i<<endl;
	  //	  fprint_axis(stderr,tmp->getAxis(i));
	  if (!grid_union(&g,&(tmp->getAxis(i)))) {
	    RVLException e;
	    e<<"Error: IWaveSim constructor from grid_union\n";
	    fprintf(stream,"Error: IWaveSim constructor from grid_union\n");
	    fprintf(stream,"  failed to add these axes:\n");
	    for (int j=0;j<tmp->getNumAxes(); j++) 
	      fprint_axis(stream,tmp->getAxis(j));
	    fprintf(stream,"  to grid:\n");
	    fprint_grid(stream,g);
	    throw e;
	  }
	}
	//	cerr<<"grid"<<endl;
	//	fprint_grid(stderr,g);
	// for this test no need to build sampler vector
	//	s.push_back(tmp);
	delete tmp;
      }

      EXPECT_EQ(4,g.gdim);
      EXPECT_EQ(416,g.axes[0].n);
      EXPECT_EQ(25.0,g.axes[0].d);
      EXPECT_EQ(0.0,g.axes[0].o);
      EXPECT_EQ(0,g.axes[0].id);
      EXPECT_EQ(800,g.axes[1].n);
      EXPECT_EQ(25.0,g.axes[1].d);
      EXPECT_EQ(0.0,g.axes[1].o);
      EXPECT_EQ(1,g.axes[1].id);
      EXPECT_EQ(284,g.axes[2].n);
      EXPECT_GT(1.e-6,fabs(g.axes[2].d-1.767767));
      EXPECT_GT(1.e-4,fabs(-100.7627-g.axes[2].o));
      EXPECT_EQ(2,g.axes[2].id);
      EXPECT_EQ(3,g.axes[3].n);
      EXPECT_EQ(1.0,g.axes[3].d);
      EXPECT_EQ(0.0,g.axes[3].o);
      EXPECT_EQ(3,g.axes[3].id);

      delete wt;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, spacetime_grid_build_ext_after) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // modify - set csq=csq2.rsf
      ps_slcstring(*par,"csq","csq2.rsf");

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;

      // initialize grid
      grid g;
      init_default_grid(&g);

      // iotask list
      std::vector<TASK_RELN *> t;
      // step 1: create list of i/o tasks
      IOTask(t,order,fwd,ic);
      IWaveTree * wt = new IWaveTree(*par,stream,ic,order);

      for (int it=0; it<t.size(); it++) {
	IWaveSampler * tmp = NULL;
	//	cerr<<"constructing sampler on "<<t[it]->keyword<<endl;
	tmp = new IWaveSampler(wt->getStateArray()[0], t[it]->keyword, *par, stream);
	// note these are internal simulation axes, not archival, 
	// which will require some re-arranging in traceio
	for (int i=0;i<tmp->getNumAxes();i++) {
	  //	  cerr<<"axis "<<i<<endl;
	  //	  fprint_axis(stderr,tmp->getAxis(i));
	  if (!grid_union(&g,&(tmp->getAxis(i)))) {
	    RVLException e;
	    e<<"Error: IWaveSim constructor from grid_union\n";
	    fprintf(stream,"Error: IWaveSim constructor from grid_union\n");
	    fprintf(stream,"  failed to add these axes:\n");
	    for (int j=0;j<tmp->getNumAxes(); j++) 
	      fprint_axis(stream,tmp->getAxis(j));
	    fprintf(stream,"  to grid:\n");
	    fprint_grid(stream,g);
	    throw e;
	  }
	}
	// for this test no need to build sampler vector
	//	s.push_back(tmp);
	delete tmp;
      }

      EXPECT_EQ(4,g.gdim);
      EXPECT_EQ(416,g.axes[0].n);
      EXPECT_EQ(25.0,g.axes[0].d);
      EXPECT_EQ(0.0,g.axes[0].o);
      EXPECT_EQ(0,g.axes[0].id);
      EXPECT_EQ(800,g.axes[1].n);
      EXPECT_EQ(25.0,g.axes[1].d);
      EXPECT_EQ(0.0,g.axes[1].o);
      EXPECT_EQ(1,g.axes[1].id);
      EXPECT_EQ(284,g.axes[2].n);
      EXPECT_GT(1.e-6,fabs(g.axes[2].d-1.767767));
      EXPECT_GT(1.e-4,fabs(-100.7627-g.axes[2].o));
      EXPECT_EQ(2,g.axes[2].id);
      EXPECT_EQ(3,g.axes[3].n);
      EXPECT_EQ(1.0,g.axes[3].d);
      EXPECT_EQ(0.0,g.axes[3].o);
      EXPECT_EQ(3,g.axes[3].id);

      delete wt;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, spacetime_grid_build_ext_before) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // modify - set csq=csq3.rsf
      ps_slcstring(*par,"csq","csq3.rsf");

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;

      // initialize grid
      grid g;
      init_default_grid(&g);

      // iotask list
      std::vector<TASK_RELN *> t;
      // step 1: create list of i/o tasks

      IOTask(t,order,fwd,ic);

      IWaveTree * wt = new IWaveTree(*par,stream,ic,order);

      for (int it=0; it<t.size(); it++) {
	IWaveSampler * tmp = NULL;
	//	cerr<<"constructing sampler on "<<t[it]->keyword<<endl;
	tmp = new IWaveSampler(wt->getStateArray()[0], t[it]->keyword, *par, stream);
	// note these are internal simulation axes, not archival, 
	// which will require some re-arranging in traceio
	for (int i=0;i<tmp->getNumAxes();i++) {
	  if (!grid_union(&g,&(tmp->getAxis(i)))) {
	    RVLException e;
	    e<<"Error: IWaveSim constructor from grid_union\n";
	    fprintf(stream,"Error: IWaveSim constructor from grid_union\n");
	    fprintf(stream,"  failed to add these axes:\n");
	    for (int j=0;j<tmp->getNumAxes(); j++) 
	      fprint_axis(stream,tmp->getAxis(j));
	    fprintf(stream,"  to grid:\n");
	    fprint_grid(stream,g);
	    throw e;
	  }
	}
	// for this test no need to build sampler vector
	//	s.push_back(tmp);
	delete tmp;
      }

      //      cerr<<"after grid build:\n";
      //      fprint_grid(stderr,g);
      EXPECT_EQ(4,g.gdim);
      EXPECT_EQ(416,g.axes[0].n);
      EXPECT_EQ(25.0,g.axes[0].d);
      EXPECT_EQ(0.0,g.axes[0].o);
      EXPECT_EQ(0,g.axes[0].id);
      EXPECT_EQ(800,g.axes[1].n);
      EXPECT_EQ(25.0,g.axes[1].d);
      EXPECT_EQ(0.0,g.axes[1].o);
      EXPECT_EQ(1,g.axes[1].id);
      EXPECT_EQ(284,g.axes[2].n);
      EXPECT_GT(1.e-6,fabs(g.axes[2].d-1.767767));
      EXPECT_GT(1.e-4,fabs(-100.7627-g.axes[2].o));
      EXPECT_EQ(2,g.axes[2].id);
      EXPECT_EQ(3,g.axes[3].n);
      EXPECT_EQ(1.0,g.axes[3].d);
      EXPECT_EQ(0.0,g.axes[3].o);
      EXPECT_EQ(3,g.axes[3].id);

      delete wt;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, construct_sim_fwd_order0) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check  
      int order=0;
      bool fwd = true;
      int printact=0; int snaps=0;
      
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps);
#ifdef GTEST_VERBOSE
      sim->write(cerr);
#endif
      string tstr = "IWaveSim forward simulator\n  derivative order = 0\n  number of dynamic arrays in each RDOM = 2\n  number of i/o tasks = 3\n";
      stringstream ostr;
      sim->write(ostr);
      EXPECT_EQ(tstr,ostr.str());

      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, construct_sim_fwd_order0_movie) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=initfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;
      
      //      ps_slcstring(*par,"movie","movie.rsf");

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;
      int printact=0; int snaps=0;
      
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps);
#ifdef GTEST_VERBOSE
      sim->write(cerr);
#endif
      string tstr = "IWaveSim forward simulator\n  derivative order = 0\n  number of dynamic arrays in each RDOM = 2\n  number of i/o tasks = 3\n";
      stringstream ostr;
      sim->write(ostr);
      EXPECT_EQ(tstr,ostr.str());

      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, construct_sim_fwd_order1) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=1;
      bool fwd = true;
      int printact=0; int snaps=0;
      
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps);
#ifdef GTEST_VERBOSE
      sim->write(cerr);
#endif
      string tstr = "IWaveSim forward simulator\n  derivative order = 1\n  number of dynamic arrays in each RDOM = 2\n  number of i/o tasks = 4\n";
      stringstream ostr;
      sim->write(ostr);
      EXPECT_EQ(tstr,ostr.str());
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, construct_sim_adj_order01) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 

      bool fwd = false;
      stringstream err1;
      stringstream err2;
      stringstream err3;
      
      try {
	int order=0;
	int printact=0; int snaps=0;
	IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps);
	delete sim;
      }
      catch (RVLException & e) {
	e.write(err1);
      }

      try {
	int order=0;
	int printact=0; int snaps=1;
	IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps);
	delete sim;
      }
      catch (RVLException & e) {
	e.write(err2);
      }

      try {
	int order=1;
	int printact=0; int snaps=0;
	IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps);
	delete sim;
      }
      catch (RVLException & e) {
	e.write(err3);
      }

      string terr1 = "Error: IWaveSim constructor\n  adjoint mode (fwd=false):\n  must provide positive number of \n  storage units for checkpoints (argument \"snaps\")\n  however snaps = 0\n  does not make sense for reference simulation (order=0)\n  order supplied in argument list = 0\n\ncalled from IWaveSim constructor\n\n";
      string terr2 = "Error: IWaveSim constructor\n  adjoint mode (fwd=false):\n  does not make sense for reference simulation (order=0)\n  order supplied in argument list = 0\n\ncalled from IWaveSim constructor\n\n";
      string terr3 = "Error: IWaveSim constructor\n  adjoint mode (fwd=false):\n  must provide positive number of \n  storage units for checkpoints (argument \"snaps\")\n  however snaps = 0\n\ncalled from IWaveSim constructor\n\n";
      EXPECT_EQ(terr1,err1.str());
      EXPECT_EQ(terr2,err2.str());
      EXPECT_EQ(terr3,err3.str());

      int order=1;
      int printact=0; int snaps=5;
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps);
#ifdef GTEST_VERBOSE
      sim->write(cerr);
#endif
      string tstr = "IWaveSim adjoint simulator\n  derivative order = 1\n  number of dynamic arrays in each RDOM = 2\n  number of i/o tasks = 4\n  number of checkpoint states = 5\n  number of checkpoint RARRs  = 10\n";
      stringstream ostr;
      sim->write(ostr);
      EXPECT_EQ(tstr,ostr.str());
      delete sim;

      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, construct_sim_fwd_order2) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=2;
      bool fwd = true;
      int printact=0; int snaps=0;
      
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps);
#ifdef GTEST_VERBOSE
      sim->write(cerr);
#endif
      string tstr = "IWaveSim forward simulator\n  derivative order = 2\n  number of dynamic arrays in each RDOM = 2\n  number of i/o tasks = 5\n";
      stringstream ostr;
      sim->write(ostr);
      EXPECT_EQ(tstr,ostr.str());
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, construct_sim_adj_order2) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=2;
      bool fwd = false;
      int printact=0; int snaps=5;
      
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps);
#ifdef GTEST_VERBOSE
      sim->write(cerr);
#endif
      string tstr = "IWaveSim adjoint simulator\n  derivative order = 2\n  number of dynamic arrays in each RDOM = 2\n  number of i/o tasks = 5\n  number of checkpoint states = 5\n  number of checkpoint RARRs  = 20\n";
      stringstream ostr;
      sim->write(ostr);
      EXPECT_EQ(tstr,ostr.str());
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_sim_fwd_ord0) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;
      int printact=0; int snaps=0;
      
      bool dryrun=true;
      ofstream drystr("dryrun_sim_fwd_ord0");
      cerr<<"1\n";
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr,cerr);
      cerr<<"2\n";
      sim->run();
      cerr<<"3\n";
      drystr.close();
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_sim_fwd_ord0_swapxz) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      ps_slcstring(*par,"csq","csq1.rsf");

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;
      int printact=0; int snaps=0;
      
      bool dryrun=true;
      ofstream drystr("dryrun_sim_fwd_ord0_swapxz");
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr);
      sim->run();
      drystr.close();
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_sim_fwd_ord0_ext_sx) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      ps_slcstring(*par,"csq","csq2.rsf");

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;
      int printact=0; int snaps=0;
      
      bool dryrun=true;
      ofstream drystr("dryrun_sim_fwd_ord0_ext_sx");
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr);
      sim->run();
      drystr.close();
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_sim_fwd_ord0_ext_hx) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      ps_slcstring(*par,"csq","csq3.rsf");

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;
      int printact=0; int snaps=0;
      
      bool dryrun=true;
      ofstream drystr("dryrun_sim_fwd_ord0_hx");
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr);
      sim->run();
      drystr.close();
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_sim_fwd_ord0_movie) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=initfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;
      int printact=0; int snaps=0;
      
      bool dryrun=true;
      ofstream drystr("dryrun_sim_fwd_ord0_movie");
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr);
      sim->run();
      drystr.close();
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_sim_fwd_ord0_init) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=initfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=0;
      bool fwd = true;
      int printact=0; int snaps=0;
      
      bool dryrun=true;
      ofstream drystr("dryrun_sim_fwd_ord0_init");
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr);
      sim->run();
      drystr.close();
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_sim_fwd_ord1) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=1;
      bool fwd = true;
      int printact=0; int snaps=0;
      
      bool dryrun=true;
      ofstream drystr("dryrun_sim_fwd_ord1");
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr);
      sim->run();
      drystr.close();
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_sim_adj_ord1) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=1;
      bool fwd = false;
      int printact=0; int snaps=5;
      
      bool dryrun=true;
      ofstream drystr("dryrun_sim_adj_ord1");
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr);
      sim->run();
      drystr.close();
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_sim_fwd_ord2) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=2;
      bool fwd = true;
      int printact=0; int snaps=0;
      
      bool dryrun=true;
      ofstream drystr("dryrun_sim_fwd_ord2");
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr);
      sim->run();
      drystr.close();
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_sim_adj_ord2) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=parfile");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      // build order zero IWaveTree, check 
      int order=2;
      bool fwd = false;
      int printact=0; int snaps=5;
      
      bool dryrun=true;
      ofstream drystr("dryrun_sim_adj_ord2");
      IWaveSim * sim = new IWaveSim(order,fwd,*par,stream,ic,printact,snaps,dryrun,drystr);
      sim->run();
      drystr.close();
      delete sim;
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

}
int xargc;
char **xargv;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  int err = RUN_ALL_TESTS();
  iwave_fdestroy();
  return err;
}
