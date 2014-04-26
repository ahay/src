#include "gtest/gtest.h"
#include "acd_defn.hh"
#include "grid.h"
#include "istate.hh"
#include "iwop.hh"
#include "gridpp.hh"
#include "functions.hh"

//#define GTEST_VERBOSE
IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",    0, true,  true },
  {"uc_in",  1, true, true},
  {"up_in",  2, true, true},
  {"uc_out",   1, false,  true},
  {"up_out",   2, false,  true},
  {"",       0, false, false}
};

namespace {

  using RVL::parse;
  using RVL::RVLException;
  using RVL::Vector;
  using RVL::Components;
  using RVL::Operator;
  using RVL::OperatorEvaluation;
  using RVL::LinearOp;
  using RVL::SymmetricBilinearOp;
  using RVL::AssignFilename;
  using RVL::AssignParams;
  using RVL::RVLRandomize;
  using TSOpt::IWaveEnvironment;
  using TSOpt::IWaveTree;
  using TSOpt::IWaveSampler;
  using TSOpt::IWaveSim;
  using TSOpt::TASK_RELN;
  using TSOpt::IOTask;
  using TSOpt::IWaveOp;
  using TSOpt::GridSpace;

  void create_hfile(string hfile, grid g, float val, bool var=false) {
    
    if (hfile.size()>120) {
      RVLException e;
      e<<"Error: create_hfile\n";
      e<<"filename "<<hfile<<" longer than 120 chars\n";
      throw e;
    }
    
    string dfile = hfile + "@";

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

    IWaveInfo ic;

    ACDSimTest(): ic() {

      ofstream qs("movie_onestep.par");
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
      qs<<"           csq = csq.rsf\n";
      qs<<"          uc_in = gauss.rsf\n";
      qs<<"          up_in = gauss.rsf\n";
      qs<<"         uc_out = movie_onestep_space.rsf\n";
      qs<<"         up_out = movie_onestep_space.rsf\n";
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
      qs<<"      dump_pi = 0           dump parallel/dom. decomp info\n";
      qs<<"     dump_lda = 0           dump grid data for allocated arrays\n";
      qs<<"     dump_ldc = 0           dump grid data for computational arrays\n";
      qs<<"     dump_lds = 0           dump grid data for send arrays\n";
      qs<<"     dump_ldr = 0           dump grid data for receive arrays\n";
      qs<<"    dump_term = 0           dump terminator data\n";
      qs<<"    dump_pars = 1           print parameter table in IWaveOp\n";
      qs<<"   dump_steps = 0           print major steps in IWaveOp\n";
      qs.flush();
      qs.close();

      ofstream ps("movie.par");
      ps<<"INPUT DATA FOR IWAVE\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"FD:\n";
      ps<<"\n";
      ps<<"         order = 2           scheme half-order\n";
      ps<<"           cfl = 0.5        cfl number - frac of max stable\n";
      ps<<"          cmin = 1.0         min velocity - checked\n";
      ps<<"          cmax = 2.0         max velocity - checked\n";
      ps<<"      max_step = 0           1 = set adaptively, 0 = use standard cfl from cmax\n";
      ps<<"         fpeak = 0.010       nominal central frequency \n";
      ps<<"            dt = 2.0";
      ps<<"\n";
      ps<<"------------------------------------------------------------------------\n";
      ps<<"Model info:\n";
      ps<<"\n";
      ps<<"           csq = csq.rsf\n";
      ps<<"         uc_in = gauss.rsf\n";
      ps<<"         up_in = gauss.rsf\n";
      ps<<"        uc_out = movie_space.rsf\n";
      ps<<"        up_out = movie_space.rsf\n";
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
      ps<<"Output info:\n";
      ps<<"\n";
      ps<<"     printact = 1           per-time-step verbosity level\n";
      ps<<"                            0 - none\n";
      ps<<"                            1 - time step index\n";
      ps<<"                            2 - internal time step info\n";
      ps<<"                            > 5: dump everything\n";
      ps<<"      dump_pi = 0           dump parallel/dom. decomp info\n";
      ps<<"     dump_lda = 0           dump grid data for allocated arrays\n";
      ps<<"     dump_ldc = 0           dump grid data for computational arrays\n";
      ps<<"     dump_lds = 0           dump grid data for send arrays\n";
      ps<<"     dump_ldr = 0           dump grid data for receive arrays\n";
      ps<<"    dump_term = 0           dump terminator data\n";
      ps<<"    dump_pars = 1           print parameter table in IWaveOp\n";
      ps<<"   dump_steps = 0           print major steps in IWaveOp\n";
      ps.flush();
      ps.close();

      grid g;
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

      float val = 1.5;
      create_hfile("csq.rsf", g, val);

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
      g.axes[2].o=2.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 2;

      val = 0.0;
      create_hfile("movie_onestep_space.rsf", g, 1.0, true);

      g.dim=2;
      g.gdim=3;
      g.axes[0].n=416;
      g.axes[1].n=800;
      g.axes[2].n=21;
      g.axes[0].d=25.0;
      g.axes[1].d=25.0;
      g.axes[2].d=200.0;
      g.axes[0].o=0.0;
      g.axes[1].o=0.0;
      g.axes[2].o=0.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 2;

      val = 0.0;
      create_hfile("movie_space.rsf", g, 1.0, true);

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

  TEST_F(ACDSimTest, dryrun_op_fwd_onestep_ord0) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie_onestep.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;

      bool dryrun=true;
      ofstream drystr("dryrun_op_fwd_onestep_ord0");

      IWaveOp op(*pars,stream,dryrun,drystr);

      Vector<ireal> m(op.getDomain());
      Vector<ireal> d(op.getRange());
      
      Components<ireal> cm(m);
      Components<ireal> cd(d);

      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);
      AssignFilename mfn1("gauss.rsf");
      cm[1].eval(mfn1);
      cm[2].zero();
      
      AssignFilename dfn0("uc_out_onestep.rsf");
      cd[0].eval(dfn0);
      AssignFilename dfn1("up_out_onestep.rsf");
      cd[1].eval(dfn1);

      OperatorEvaluation<ireal> opeval(op,m);
      d.copy(opeval.getValue());

      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_op_fwd_onestep_ord1) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie_onestep.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;

      bool dryrun=true;
      ofstream drystr("dryrun_op_fwd_onestep_ord1");

      IWaveOp op(*pars,stream,dryrun,drystr);

      Vector<ireal> m(op.getDomain());
      Vector<ireal> dm(op.getDomain());
      Vector<ireal> dd(op.getRange());
      
      Components<ireal> cm(m);

      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);
      AssignFilename mfn1("gauss.rsf");
      cm[1].eval(mfn1);
      cm[2].zero();

      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
      dm.eval(rnd);
      dd.zero();

      OperatorEvaluation<ireal> opeval(op,m);
      opeval.getDeriv().applyOp(dm,dd);

      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_op_adj_onestep_ord1) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie_onestep.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;

      bool dryrun=true;
      ofstream drystr("dryrun_op_adj_onestep_ord1");

      IWaveOp op(*pars,stream,dryrun,drystr);

      Vector<ireal> m(op.getDomain());
      Vector<ireal> dm(op.getDomain());
      Vector<ireal> dd(op.getRange());
      
      Components<ireal> cm(m);
      Components<ireal> cdm(dm);
      Components<ireal> cdd(dd);

      string csqname = "csq.rsf";
      string dcsqname = "migcsq.rsf";
      string bmoviename = "bmovie_onestep.rsf";
      
      AssignFilename mfn(csqname);
      cm[0].eval(mfn);
      AssignFilename dmfn(dcsqname);
      cdm[0].eval(dmfn);
      AssignFilename ddfn(bmoviename);
      cdd[0].eval(ddfn);

      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
      dd.eval(rnd);
      dm.zero();

      OperatorEvaluation<ireal> opeval(op,m);
      opeval.getDeriv().applyAdjOp(dd,dm);

      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_op_fwd_ord0) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;

      bool dryrun=true;
      ofstream drystr("dryrun_op_fwd_ord0");

      string csqname = "csq.rsf";
      string moviename = "movie.rsf";
      
      IWaveOp op(*pars,stream,dryrun,drystr);

      Vector<ireal> m(op.getDomain());
      Vector<ireal> d(op.getRange());
      
      Components<ireal> cm(m);
      Components<ireal> cd(d);

      AssignFilename mfn(csqname);
      cm[0].eval(mfn);
      AssignFilename dfn(moviename);
      cd[0].eval(dfn);

      OperatorEvaluation<ireal> opeval(op,m);
      d.copy(opeval.getValue());

      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_op_fwd_ord1) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;

      bool dryrun=true;
      ofstream drystr("dryrun_op_fwd_ord1");

      IWaveOp op(*pars,stream,dryrun,drystr);

      Vector<ireal> m(op.getDomain());
      Vector<ireal> dm(op.getDomain());
      Vector<ireal> dd(op.getRange());
      
      Components<ireal> cm(m);
      Components<ireal> cdm(dm);
      Components<ireal> cdd(dd);

      string csqname = "csq.rsf";
      string dcsqname = "randcsq.rsf";
      string dmoviename = "dmovie.rsf";
      
      AssignFilename mfn(csqname);
      cm[0].eval(mfn);
      AssignFilename dmfn(dcsqname);
      cdm[0].eval(dmfn);
      AssignFilename ddfn(dmoviename);
      cdd[0].eval(ddfn);

      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
      dm.eval(rnd);
      dd.zero();

      OperatorEvaluation<ireal> opeval(op,m);
      opeval.getDeriv().applyOp(dm,dd);

      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDSimTest, dryrun_op_adj_ord1) {
    try {

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;

      bool dryrun=true;
      ofstream drystr("dryrun_op_adj_ord1");

      IWaveOp op(*pars,stream,dryrun,drystr);

      Vector<ireal> m(op.getDomain());
      Vector<ireal> dm(op.getDomain());
      Vector<ireal> dd(op.getRange());
      
      Components<ireal> cm(m);
      Components<ireal> cdm(dm);
      Components<ireal> cdd(dd);

      string csqname = "csq.rsf";
      string dcsqname = "migcsq.rsf";
      string bmoviename = "bmovie.rsf";
      
      AssignFilename mfn(csqname);
      cm[0].eval(mfn);
      AssignFilename dmfn(dcsqname);
      cdm[0].eval(dmfn);
      AssignFilename ddfn(bmoviename);
      cdd[0].eval(ddfn);

      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
      dd.eval(rnd);
      dm.zero();

      OperatorEvaluation<ireal> opeval(op,m);
      opeval.getDeriv().applyAdjOp(dd,dm);

      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }
}

int xargc; 
char ** xargv;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  int err = RUN_ALL_TESTS();
  iwave_fdestroy();
  return err;
}
