#include "gtest/gtest.h"
#include "acd_defn.hh"
#include "traceio.h"
#include "grid.h"
#include "except.hh"
#include "istate.hh"
#include "par.h"

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

  //  using RVL::parse;
  using RVL::RVLException;
  using TSOpt::IWaveEnvironment;

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
    cmd << "| " + CWPROOT + "/bin/sushw key=gelev a="<<(-1.0f)*rz;
    cmd << "| " + CWPROOT + "/bin/sushw key=scalel a="<<scalel;
    cmd	<< "| " + CWPROOT + "/bin/sushw key=scalco a="<<scalco;
    cmd << "| " + CWPROOT + "/bin/suchw key1=offset key2=gx key3=sx b=1 c=-1 > "<<fnm<<"\n";
    //    cerr<<"create_fixed_2D_data: cmd = "<<cmd.str()<<endl;
    if (system(cmd.str().c_str())) {
      RVLException e;
      e<<"Error: create_fixed_2D_data\n";
      e<<"  failed system call on command \n";
      e<<cmd.str();
      throw e;
    }

    //    cerr<<"create_fixed_2D_data: exit"<<endl;    
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

  class SEGY_Test : public ::testing::Test {
  public:

    IWaveInfo ic;

    SEGY_Test() {
      
      try {

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

	// first data set - single gather
	nt = 1501;
	dt = 0.002;
	ot = 0.0;
	nrec = 1;
	ntr_in_rec = 301;
	sx0=3300;
	dsx=0;
	sz=40;
	rx0=100;
	drx=20;
	rz=20;
	scalel=0;
	scalco=0;
	
	string fn1 = "data1.su";
	
	create_fixed_2D_data(fn1,nt,dt,ot,
			     nrec,ntr_in_rec,
			     sx0,dsx,sz,
			     rx0,drx,rz,
			     scalel,scalco);

	grid g;
	g.dim=2;
	g.gdim=2;
	g.axes[0].n=91;
	g.axes[1].n=391;
	g.axes[0].d=20.0;
	g.axes[1].d=20.0;
	g.axes[0].o=0.0;
	g.axes[1].o=0.0;
	g.axes[0].id = 0;
	g.axes[1].id = 1;

	string hfile = "vp2d_20m.rsf";
	string dfile = "vp2d_20m.rsf@";
	float val = 1.5;
	  
	create_hfile(hfile,dfile,g,val);

	// create simple par file
	ofstream ps("segy_test_pars");
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
	ps<<"           csq = vp2d_20m.rsf\n";
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
	ps<<"     \"   source = wavelet1.su\n\"";
	ps<<"       sampord = 1             sampling order\n";
	ps<<"\n";
	ps<<"------------------------------------------------------------------------\n";
	ps<<"Trace info:\n";
	ps<<"\n";
	ps<<"            data = data1.su    output data file\n";
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


      }
      catch (RVLException & e) {
	e<<"\ncalled from SEGY_Test constructor\n";
	throw e;
      }
    }
  };

  TEST_F(SEGY_Test, construct_nrec1) {
    try {

      tracegeom * tg = new tracegeom;
      setnull_tracegeom(tg);
      string dataname = "data1.su";
      float model_dt = 2.0;
      float src_tol = 0.001;
      int err = construct_tracegeom(tg,
				    dataname.c_str(),
				    model_dt,
				    src_tol,
				    stderr);
      if (err) {
	RVLException e;
	e<<"Error: TEST_F(SEGY_Test, construct_nrec1)\n";
	e<<"  from construct_tracegeom, err = "<<err<<"\n";
	throw e;
      }
      // see what we've got
      FILE * fp = fopen("dryrun_segy_test0","w");
      fprint_tracegeom(tg,fp);
      fclose(fp);
    }
    catch (RVLException & e) {
      e<<"\ncalled from TEST_F(SEGY_Test, construct_nrec1)\n";
      throw e;
    }
  }

  TEST_F(SEGY_Test, init_nrec1) {
    try {
      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=segy_test_pars");
      PARARRAY * par = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &par, &stream);
      delete [] argv;
      delete [] argvv;

      IWAVE w;
      IWAVE * state = &w;
      if (int err = iwave_construct(state,par,stream,ic)) {
	RVLException e;
	e<<"Error: SEGY_Test - isample_construct_rsf_input\n";
	e<<"  error from iwave_construct = "<<err<<"\n";
	throw e;
      }

      tracegeom * tg = new tracegeom;
      setnull_tracegeom(tg);
      string dataname = "data1.su";
      float src_tol = 0.001;
      int err = construct_tracegeom(tg,
				    dataname.c_str(),
				    (state->model).tsind.dt,
				    src_tol,
				    stderr);
      if (err) {
	RVLException e;
	e<<"Error: TEST_F(SEGY_Test, construct_nrec1)\n";
	e<<"  from construct_tracegeom, err = "<<err<<"\n";
	throw e;
      }

      IPNT m_n;                       /* axis lengths, local grid */
      RPNT m_o;                       /* axis origins, local grid */
      RPNT m_og;                      /* axis origins, global grid */
      RPNT m_d;                       /* axis steps, local grid */
      IPNT m_axord;                   /* axis order array */
      
      /* extract grid params - space and time */
      get_n(m_n,(state->model).gl);
      get_o(m_o,(state->model).gl);
      get_d(m_d,(state->model).gl);
      get_o(m_og,(state->model).g);
      get_ord(m_axord,(state->model).g);
      
      int sampord = 1;
      int initbuf = 0;
      int irec    = 0;

      /* initialize tracegeom */
      err=init_tracegeom(tg,
			 irec,
			 m_og,m_n,m_d,m_o,m_axord,
			 sampord,
			 (state->model).g.dim,
			 initbuf,
			 stream);
      
      // see what we've got
      FILE * fp = fopen("dryrun_segy_test0","w");
      fprint_tracegeom(tg,fp);
      fclose(fp);
    }
    catch (RVLException & e) {
      e<<"\ncalled from TEST_F(SEGY_Test, construct_nrec1)\n";
      throw e;
    }
  }
}

int xargc;
char ** xargv;
int main(int argc, char **argv) {
  //  xargc = argc; xargv = argv;
  try {
    ::testing::InitGoogleTest(&argc, argv);
    int err = RUN_ALL_TESTS();
    iwave_fdestroy();
    return err;
  }
  catch (RVLException &e) {
    e.write(cerr);
    exit(1);
  }
}

				    
      
