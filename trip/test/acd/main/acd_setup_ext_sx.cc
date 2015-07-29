#include "gtest/gtest.h"
#include "acd_defn.hh"
#include "grid.h"
#include "traceio.h"
#include "iwsim.hh"

//#define GTEST_VERBOSE

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",    0, true,  true },
  {"data",   1, false, true },
  {"source", 1, true,  false},
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

      hfile = "csq_ext_sx.rsf";
      dfile = "csq_ext_sx.rsf@";

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
      ps<<"           csq = "<<hfile<<"\n";
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
      ps<<"        source = "<<wn<<"\n";
      ps<<"       sampord = 1             sampling order\n";
      ps<<"\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"Trace info:\n";
      ps<<"\n";
      ps<<"            data = "<<dn<<"    output data file\n";
      ps<<"\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"Output info:\n";
      ps<<"\n";
      ps<<"     printact = 1           per-time-step verbosity level\n";
      ps<<"                            0 - none\n";
      ps<<"                            1 - time step index\n";
      ps<<"                            2 - internal time step info\n";
      ps<<"                            > 5: dump everything\n";
      ps<<"      dump_pi = 0           dump parallel/dom. decomp info\n";
      ps<<"     dump_lda = 1           dump grid data for allocated arrays\n";
      ps<<"     dump_ldc = 1           dump grid data for computational arrays\n";
      ps<<"     dump_lds = 0           dump grid data for send arrays\n";
      ps<<"     dump_ldr = 0           dump grid data for receive arrays\n";
      ps<<"    dump_term = 1           dump terminator data\n";
      ps<<"    dump_pars = 0           print parameter table in IWaveOp\n";
      ps<<"   dump_steps = 0           print major steps in IWaveOp\n";
      ps.flush();
      ps.close();

      grid g;

      g.gdim=3;
      g.dim=2;
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
      float val = 1.5;
      create_hfile(hfile, dfile, g, val);

    }
  };

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


}
int xargc;
char **xargv;

int main(int argc, char **argv) {
  int ts=0;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
  try {
    ::testing::InitGoogleTest(&argc, argv);
    int res = RUN_ALL_TESTS();
    MPI_Finalize();
    return res;
  }
  catch (RVLException & e) {
    e.write(cerr);
    MPI_Abort(MPI_COMM_WORLD,0);
    MPI_Finalize();
    exit(1);
  }
}

