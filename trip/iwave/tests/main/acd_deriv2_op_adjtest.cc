#include "gtest/gtest.h"
#include "acd_defn.hh"
#include "grid.h"
#include "istate.hh"
#include "iwop.hh"
#include "gridpp.hh"
#include "functions.hh"
#include "adjtest.hh"

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
  using RVL::LinearBilinearOp;
  using RVL::SymmetricBilinearOp;
  using RVL::AssignFilename;
  using RVL::AssignParams;
  using RVL::RVLRandomize;
  using RVL::AdjointTest;
  using TSOpt::IWaveEnvironment;
  using TSOpt::IWaveTree;
  using TSOpt::IWaveSampler;
  using TSOpt::IWaveSim;
  using TSOpt::TASK_RELN;
  using TSOpt::IOTask;
  using TSOpt::IWaveOp;
  using TSOpt::GridSpace;

  void constructargs(std::vector<std::string> args,
		     int * argc,
		     char ***argv) {
    if (args.size() < 2) {
      RVLException e;
      e<<"Error: constructargs\n";
      e<<"  input std::vector must have length > 1\n";
      throw e;
    }
    *argc = args.size()+1;
    *argv = new char*[*argc];
    (*argv)[0]=NULL;
    for (int i=0;i<*argc-1;i++) {
      (*argv)[i+1]=new char[args[i].size()+1];
      strcpy((*argv)[i+1],args[i].c_str());
      //      cerr<<"argv["<<i+1<<"]="<<(*argv)[i+1]<<"\n";
    }
  }

  void destroyargs(int argc, char ***argv) {
    for (int i=0;i<argc-1;i++) delete (*argv)[i+1];
    delete [] *argv;
  }

  void create_hfile(string hfile, grid g, float val, bool var=false, int position=0) {
    
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
      if(position==0){
    if (var) {
      for (int i=0;i<get_datasize_grid(g);i++) buf[i]=val/((float)(i+1));
    }
    else {
      for (int i=0;i<get_datasize_grid(g);i++) buf[i]=val;
    }}
      else{
          for (int i=0;i<get_datasize_grid(g);i++) buf[i]=val;
          buf[position] = 1;
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

  class ACDCauchyOpAdjTest : public ::testing::Test {
  public:

    IWaveInfo ic;
    ofstream rpt;

    ACDCauchyOpAdjTest(): ic() {

      rpt.open("detail.txt",ios_base::app);

      ofstream qs("movie_1step.par");
      qs<<"INPUT DATA FOR IWAVE\n";
      qs<<"------------------------------------------------------------------------\n";
      qs<<"FD:\n";
      qs<<"\n";
      qs<<"         order = 1           scheme half-order\n";
      qs<<"           cfl = 0.5        cfl number - frac of max stable\n";
      qs<<"          cmin = 1.0         min velocity - checked\n";
      qs<<"          cmax = 2.0         max velocity - checked\n";
      qs<<"            dt = 2.0";
      qs<<"\n";
      qs<<"------------------------------------------------------------------------\n";
      qs<<"Model info:\n";
      qs<<"\n";
      qs<<"           csq = csq.rsf\n";
      qs<<"          uc_in = u_in.rsf\n";
      qs<<"          up_in = u_in.rsf\n";
      qs<<"         uc_out = u_out.rsf\n";
      qs<<"         up_out = u_out.rsf\n";
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
      qs<<"     printact = 0           per-time-step verbosity level\n";
      qs<<"                            0 - none\n";
      qs<<"                            1 - time step index\n";
      qs<<"                            2 - internal time step info\n";
      qs<<"                            > 5: dump everything\n";
      qs<<"      dump_pi = 0           dump parallel/dom. decomp info\n";
      qs<<"     dump_lda = 1           dump grid data for allocated arrays\n";
      qs<<"     dump_ldc = 1           dump grid data for computational arrays\n";
      qs<<"     dump_lds = 0           dump grid data for send arrays\n";
      qs<<"     dump_ldr = 0           dump grid data for receive arrays\n";
      qs<<"    dump_term = 0           dump terminator data\n";
      qs<<"    dump_pars = 0           print parameter table in IWaveOp\n";
      qs<<"   dump_steps = 0           print major steps in IWaveOp\n";
      qs.flush();
      qs.close();

      qs.open("movie_3step.par");
      qs<<"INPUT DATA FOR IWAVE\n";
      qs<<"------------------------------------------------------------------------\n";
      qs<<"FD:\n";
      qs<<"\n";
      qs<<"         order = 1           scheme half-order\n";
      qs<<"           cfl = 0.5        cfl number - frac of max stable\n";
      qs<<"          cmin = 1.0         min velocity - checked\n";
      qs<<"          cmax = 2.0         max velocity - checked\n";
      qs<<"            dt = 2.0";
      qs<<"\n";
      qs<<"------------------------------------------------------------------------\n";
      qs<<"Model info:\n";
      qs<<"\n";
      qs<<"           csq = csq.rsf\n";
      qs<<"          uc_in = u_in.rsf\n";
      qs<<"          up_in = u_in.rsf\n";
      qs<<"         uc_out = u_out3.rsf\n";
      qs<<"         up_out = u_out3.rsf\n";
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
      qs<<"     printact = 0           per-time-step verbosity level\n";
      qs<<"                            0 - none\n";
      qs<<"                            1 - time step index\n";
      qs<<"                            2 - internal time step info\n";
      qs<<"                            > 5: dump everything\n";
      qs<<"      dump_pi = 0           dump parallel/dom. decomp info\n";
      qs<<"     dump_lda = 1           dump grid data for allocated arrays\n";
      qs<<"     dump_ldc = 1           dump grid data for computational arrays\n";
      qs<<"     dump_lds = 0           dump grid data for send arrays\n";
      qs<<"     dump_ldr = 0           dump grid data for receive arrays\n";
      qs<<"    dump_term = 0           dump terminator data\n";
      qs<<"    dump_pars = 0           print parameter table in IWaveOp\n";
      qs<<"   dump_steps = 0           print major steps in IWaveOp\n";
      qs.flush();
      qs.close();

      qs.open("movie_100step.par");
      qs<<"INPUT DATA FOR IWAVE\n";
      qs<<"------------------------------------------------------------------------\n";
      qs<<"FD:\n";
      qs<<"\n";
      qs<<"         order = 1           scheme half-order\n";
      qs<<"           cfl = 0.5        cfl number - frac of max stable\n";
      qs<<"          cmin = 1.0         min velocity - checked\n";
      qs<<"          cmax = 2.0         max velocity - checked\n";
      qs<<"            dt = 2.0";
      qs<<"\n";
      qs<<"------------------------------------------------------------------------\n";
      qs<<"Model info:\n";
      qs<<"\n";
      qs<<"           csq = csq.rsf\n";
      qs<<"          uc_in = u_in.rsf\n";
      qs<<"          up_in = u_in.rsf\n";
      qs<<"         uc_out = u_out100.rsf\n";
      qs<<"         up_out = u_out100.rsf\n";
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
      qs<<"     printact = 0           per-time-step verbosity level\n";
      qs<<"                            0 - none\n";
      qs<<"                            1 - time step index\n";
      qs<<"                            2 - internal time step info\n";
      qs<<"                            > 5: dump everything\n";
      qs<<"      dump_pi = 0           dump parallel/dom. decomp info\n";
      qs<<"     dump_lda = 1           dump grid data for allocated arrays\n";
      qs<<"     dump_ldc = 1           dump grid data for computational arrays\n";
      qs<<"     dump_lds = 0           dump grid data for send arrays\n";
      qs<<"     dump_ldr = 0           dump grid data for receive arrays\n";
      qs<<"    dump_term = 0           dump terminator data\n";
      qs<<"    dump_pars = 0           print parameter table in IWaveOp\n";
      qs<<"   dump_steps = 0           print major steps in IWaveOp\n";
      qs.flush();
      qs.close();

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
      g.axes[2].o=0.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 2;

      val = 0.0;
      create_hfile("u_in.rsf",g,val);

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
      create_hfile("u_out.rsf",g,val);

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
      g.axes[2].o=6.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 2;

      val = 0.0;
      create_hfile("u_out3.rsf",g,val);

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
      g.axes[2].o=200.0;
      g.axes[0].id = 0;
      g.axes[1].id = 1;
      g.axes[2].id = 2;

      val = 0.0;
      create_hfile("u_out100.rsf",g,val);

    }

    ~ACDCauchyOpAdjTest() { rpt.close(); }
  };

/*
  TEST_F(ACDCauchyOpAdjTest, dryrun_op_onestep_ord1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDCauchyOpAdjTest, dryrun_op_onestep_ord1\n\n";

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie_1step.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;
      
      bool dryrun=true;
      ofstream drystr("dryrun_op_adj_onestep_ord1");
      
      IWaveOp op(*pars,stream,dryrun,drystr,rpt);
      
      Vector<ireal> m(op.getDomain());
      
      Vector<ireal> dm(op.getDomain());
      Vector<ireal> bm(op.getDomain());
      
      Vector<ireal> dd(op.getRange());
      Vector<ireal> bd(op.getRange());
      
      int seed = 19490615;
      RVL::RVLRandomize<float> rnd(seed,-1.0f,1.0f);            
      Components<ireal> cm(m);
      Components<ireal> cdm(dm);
      Components<ireal> cdd(dd);
      Components<ireal> cbm(bm);
      Components<ireal> cbd(bd);
      
      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);
      
      cm[1].eval(rnd);
      cm[2].zero();
      
      cdm[0].eval(rnd);
      cdm[1].zero();
      cdm[2].zero();
      cbd[0].eval(rnd);
      cbd[1].zero();
      
      dd.zero();
      bm.zero();
      
      // apply Op and AdjOp
      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();
      lop.applyOp(dm,dd);
      lop.applyAdjOp(bd,bm);
      
      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDCauchyOpAdjTest, wetrun_op_uonly_onestep_ord1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDCauchyOpAdjTest, wetrun_op_uonly_onestep_ord1\n\n";
      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie_1step.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;
      
      bool dryrun=false;
      //      ofstream drystr("op_adj_onestep_ord1");
      
      IWaveOp op(*pars,stream,dryrun,rpt,rpt);
      
      Vector<ireal> m(op.getDomain());
      
      Vector<ireal> dm(op.getDomain());
      Vector<ireal> bm(op.getDomain());
      
      Vector<ireal> dd(op.getRange());
      Vector<ireal> bd(op.getRange());
      
      //      int seed = 19490615;
      int seed = getpid();
      RVL::RVLRandomize<float> rnd(seed,-1.0f,1.0f);            
      Components<ireal> cm(m);
      Components<ireal> cdm(dm);
      Components<ireal> cdd(dd);
      Components<ireal> cbm(bm);
      Components<ireal> cbd(bd);
      
      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);
      
      cm[1].eval(rnd);
      cm[2].zero();
      
      rpt << "m.norm = " << m.norm() << endl;
      rpt << "m[0].norm = " << cm[0].norm() << endl;
      rpt << "m[1].norm = " << cm[1].norm() << endl;
      rpt << "m[2].norm = " << cm[2].norm() << endl;
      
      cdm[0].zero();
      cdm[1].eval(rnd);
      cdm[2].eval(rnd);
      cbd[0].eval(rnd);
      cbd[1].eval(rnd);
      
      dd.zero();
      bm.zero();
      
      // apply Op and AdjOp
      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();
      lop.applyOp(dm,dd);
      lop.applyAdjOp(bd,bm);
      rpt << "\n=======   x ========================\n";
      rpt << "dm[0].norm = " << cdm[0].norm() << endl;
      rpt << "dm[1].norm = " << cdm[1].norm() << endl;
      rpt << "dm[2].norm = " << cdm[2].norm() << endl;
      rpt << "\n=======  Ax ========================\n";
      rpt << "dd[0].norm = " << cdd[0].norm() << endl;
      rpt << "dd[1].norm = " << cdd[1].norm() << endl;
      rpt << "\n=======   y [ucb, upb] =============\n";
      rpt << "bd[0].norm = " << cbd[0].norm() << endl;
      rpt << "bd[1].norm = " << cbd[1].norm() << endl;
      rpt << "\n=======A^Ty [csqb, ucb, upb] =======\n";
      rpt << "bm[0].norm = " << cbm[0].norm() << endl;
      rpt << "bm[1].norm = " << cbm[1].norm() << endl;
      rpt << "bm[2].norm = " << cbm[2].norm() << endl;
      
      
      rpt << "<csqd, csqb> = " << cdm[0].inner(cbm[0]) << endl;
      // compute norms and inner products
      ireal axnm, ynm, axy, xaty;
      axnm = dd.norm();
      ynm  = bd.norm();
      axy  = dd.inner(bd);
      xaty = bm.inner(dm);
      rpt << "<Ax,    y> = " << axy << endl;
      rpt << "< x, A^Ty> = " << xaty << endl;
      rpt << " |Ax| * |y| = " << axnm * ynm << endl;
      rpt << " adjvalue = " << abs(axy - xaty)/(axnm*ynm) << endl;
      
      EXPECT_GT(1.0e-6,abs(axy - xaty)/(axnm*ynm));

      //      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }
  
  TEST_F(ACDCauchyOpAdjTest, wetrun_op_conly_onestep_ord1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDCauchyOpAdjTest, wetrun_op_conly_onestep_ord1\n\n";
      
      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie_1step.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;
      
      bool dryrun=false;
      //      ofstream drystr("op_adj_onestep_ord1");
      
      IWaveOp op(*pars,stream,dryrun,rpt,rpt);
      
      Vector<ireal> m(op.getDomain());
      
      Vector<ireal> dm(op.getDomain());
      Vector<ireal> bm(op.getDomain());
      
      Vector<ireal> dd(op.getRange());
      Vector<ireal> bd(op.getRange());
      
      //      int seed = 19490615;
      int seed = getpid();
      RVL::RVLRandomize<float> rnd(seed,-1.0f,1.0f);            
      Components<ireal> cm(m);
      Components<ireal> cdm(dm);
      Components<ireal> cdd(dd);
      Components<ireal> cbm(bm);
      Components<ireal> cbd(bd);
      
      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);
      //      AssignFilename ucfn("../uc_spike.rsf");
      cm[1].eval(rnd);
      cm[2].zero();
      
      rpt << "m.norm = " << m.norm() << endl;
      rpt << "m[0].norm = " << cm[0].norm() << endl;
      rpt << "m[1].norm = " << cm[1].norm() << endl;
      rpt << "m[2].norm = " << cm[2].norm() << endl;

      //      AssignFilename csqdfn("../csq_spike.rsf");
      cdm[0].eval(rnd);
      cdm[1].zero();
      cdm[2].zero();
      //      AssignFilename ucbdfn("../ucbd_spike.rsf");
      cbd[0].eval(rnd);
      cbd[1].zero();
      
      dd.zero();
      bm.zero();
      
      // apply Op and AdjOp
      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();
      lop.applyOp(dm,dd);
      lop.applyAdjOp(bd,bm);
      rpt << "\n=======   x ========================\n";
      rpt << "dm[0].norm = " << cdm[0].norm() << endl;
      rpt << "dm[1].norm = " << cdm[1].norm() << endl;
      rpt << "dm[2].norm = " << cdm[2].norm() << endl;
      rpt << "\n=======  Ax ========================\n";
      rpt << "dd[0].norm = " << cdd[0].norm() << endl;
      rpt << "dd[1].norm = " << cdd[1].norm() << endl;
      rpt << "\n=======   y [ucb, upb] =============\n";
      rpt << "bd[0].norm = " << cbd[0].norm() << endl;
      rpt << "bd[1].norm = " << cbd[1].norm() << endl;
      rpt << "\n=======A^Ty [csqb, ucb, upb] =======\n";
      rpt << "bm[0].norm = " << cbm[0].norm() << endl;
      rpt << "bm[1].norm = " << cbm[1].norm() << endl;
      rpt << "bm[2].norm = " << cbm[2].norm() << endl;
      
      
      rpt << "<csqd, csqb> = " << cdm[0].inner(cbm[0]) << endl;
      // compute norms and inner products
      ireal axnm, ynm, axy, xaty;
      axnm = dd.norm();
      ynm  = bd.norm();
      axy  = dd.inner(bd);
      xaty = bm.inner(dm);
      rpt << "<Ax,    y> = " << axy << endl;
      rpt << "< x, A^Ty> = " << xaty << endl;
      rpt << " |Ax| * |y| = " << axnm * ynm << endl;
      rpt << " adjvalue = " << abs(axy - xaty)/(axnm*ynm) << endl;
      
      EXPECT_GT(1.0e-6,abs(axy - xaty)/(axnm*ynm));

      //      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_onestep_ord1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_onestep_ord1\n\n";

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie_1step.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;
      
      bool dryrun=false;
      //      ofstream drystr("op_adj_onestep_ord1");
      
      IWaveOp op(*pars,stream,dryrun,rpt,rpt);
      
      Vector<ireal> m(op.getDomain());
      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);            
      Components<ireal> cm(m);
      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);      
      cm[1].eval(rnd);
      cm[2].zero();

      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();

      bool res = AdjointTest(lop,rnd,rpt);

      EXPECT_EQ(res,true);

      //      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }
  
  TEST_F(ACDCauchyOpAdjTest, dryrun_op_3step_ord1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDCauchyOpAdjTest, dryrun_op_3step_ord1\n\n";

      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie_3step.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;
      
      bool dryrun=true;
      ofstream drystr("dryrun_op_adj_3step_ord1");
      
      IWaveOp op(*pars,stream,dryrun,drystr,rpt);
      
      Vector<ireal> m(op.getDomain());
      
      Vector<ireal> dm(op.getDomain());
      Vector<ireal> bm(op.getDomain());
      
      Vector<ireal> dd(op.getRange());
      Vector<ireal> bd(op.getRange());
      
      int seed = 19490615;
      RVL::RVLRandomize<float> rnd(seed,-1.0f,1.0f);            
      Components<ireal> cm(m);
      Components<ireal> cdm(dm);
      Components<ireal> cdd(dd);
      Components<ireal> cbm(bm);
      Components<ireal> cbd(bd);
      
      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);
      
      cm[1].eval(rnd);
      cm[2].zero();
      
      cdm[0].eval(rnd);
      cdm[1].zero();
      cdm[2].zero();
      cbd[0].eval(rnd);
      cbd[1].zero();
      
      dd.zero();
      bm.zero();
      
      // apply Op and AdjOp
      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();
      lop.applyOp(dm,dd);
      lop.applyAdjOp(bd,bm);
      
      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDCauchyOpAdjTest, wetrun_op_uonly_3step_ord1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDCauchyOpAdjTest, wetrun_op_uonly_3step_ord1\n\n";
      
      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie_3step.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;
      
      bool dryrun=false;
      //      ofstream drystr("op_adj_onestep_ord1");
      
      IWaveOp op(*pars,stream,dryrun,rpt,rpt);
      
      Vector<ireal> m(op.getDomain());
      
      Vector<ireal> dm(op.getDomain());
      Vector<ireal> bm(op.getDomain());
      
      Vector<ireal> dd(op.getRange());
      Vector<ireal> bd(op.getRange());
      
      //      int seed = 19490615;
      int seed = getpid();
      RVL::RVLRandomize<float> rnd(seed,-1.0f,1.0f);            
      Components<ireal> cm(m);
      Components<ireal> cdm(dm);
      Components<ireal> cdd(dd);
      Components<ireal> cbm(bm);
      Components<ireal> cbd(bd);
      
      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);
      
      cm[1].eval(rnd);
      cm[2].zero();
      
      rpt << "m.norm = " << m.norm() << endl;
      rpt << "m[0].norm = " << cm[0].norm() << endl;
      rpt << "m[1].norm = " << cm[1].norm() << endl;
      rpt << "m[2].norm = " << cm[2].norm() << endl;
      
      cdm[0].zero();
      cdm[1].eval(rnd);
      cdm[2].eval(rnd);
      cbd[0].eval(rnd);
      cbd[1].eval(rnd);
      
      dd.zero();
      bm.zero();
      
      // apply Op and AdjOp
      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();
      lop.applyOp(dm,dd);
      lop.applyAdjOp(bd,bm);
      rpt << "\n=======   x ========================\n";
      rpt << "dm[0].norm = " << cdm[0].norm() << endl;
      rpt << "dm[1].norm = " << cdm[1].norm() << endl;
      rpt << "dm[2].norm = " << cdm[2].norm() << endl;
      rpt << "\n=======  Ax ========================\n";
      rpt << "dd[0].norm = " << cdd[0].norm() << endl;
      rpt << "dd[1].norm = " << cdd[1].norm() << endl;
      rpt << "\n=======   y [ucb, upb] =============\n";
      rpt << "bd[0].norm = " << cbd[0].norm() << endl;
      rpt << "bd[1].norm = " << cbd[1].norm() << endl;
      rpt << "\n=======A^Ty [csqb, ucb, upb] =======\n";
      rpt << "bm[0].norm = " << cbm[0].norm() << endl;
      rpt << "bm[1].norm = " << cbm[1].norm() << endl;
      rpt << "bm[2].norm = " << cbm[2].norm() << endl;
      
      
      rpt << "<csqd, csqb> = " << cdm[0].inner(cbm[0]) << endl;
      // compute norms and inner products
      ireal axnm, ynm, axy, xaty;
      axnm = dd.norm();
      ynm  = bd.norm();
      axy  = dd.inner(bd);
      xaty = bm.inner(dm);
      rpt << "<Ax,    y> = " << axy << endl;
      rpt << "< x, A^Ty> = " << xaty << endl;
      rpt << " |Ax| * |y| = " << axnm * ynm << endl;
      rpt << " adjvalue = " << abs(axy - xaty)/(axnm*ynm) << endl;
      
      EXPECT_GT(1.0e-6,abs(axy - xaty)/(axnm*ynm));

      //      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }
  
  TEST_F(ACDCauchyOpAdjTest, wetrun_op_conly_3step_ord1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDCauchyOpAdjTest, wetrun_op_conly_3step_ord1\n\n";
      
      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie_3step.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;
      
      bool dryrun=false;
      //      ofstream drystr("op_adj_onestep_ord1");
      
      IWaveOp op(*pars,stream,dryrun,rpt,rpt);
      
      Vector<ireal> m(op.getDomain());
      
      Vector<ireal> dm(op.getDomain());
      Vector<ireal> bm(op.getDomain());
      
      Vector<ireal> dd(op.getRange());
      Vector<ireal> bd(op.getRange());
      
      //      int seed = 19490615;
      int seed = getpid();
      RVL::RVLRandomize<float> rnd(seed,-1.0f,1.0f);            
      Components<ireal> cm(m);
      Components<ireal> cdm(dm);
      Components<ireal> cdd(dd);
      Components<ireal> cbm(bm);
      Components<ireal> cbd(bd);
      
      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);
      //      AssignFilename ucfn("../uc_spike.rsf");
      cm[1].eval(rnd);
      cm[2].zero();
      
      rpt << "m.norm = " << m.norm() << endl;
      rpt << "m[0].norm = " << cm[0].norm() << endl;
      rpt << "m[1].norm = " << cm[1].norm() << endl;
      rpt << "m[2].norm = " << cm[2].norm() << endl;
      
      AssignFilename csqdfn("../csq_spike.rsf");
      cdm[0].eval(rnd);
      //      cdm[0].eval(csqdfn);
      cdm[1].zero();
      cdm[2].zero();
      //      AssignFilename ucbdfn("../ucbd_spike.rsf");
      cbd[0].eval(rnd);
      cbd[1].zero();
      
      dd.zero();
      bm.zero();
      
      // apply Op and AdjOp
      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();
      lop.applyOp(dm,dd);
      lop.applyAdjOp(bd,bm);
      rpt << "\n=======   x ========================\n";
      rpt << "dm[0].norm = " << cdm[0].norm() << endl;
      rpt << "dm[1].norm = " << cdm[1].norm() << endl;
      rpt << "dm[2].norm = " << cdm[2].norm() << endl;
      rpt << "\n=======  Ax ========================\n";
      rpt << "dd[0].norm = " << cdd[0].norm() << endl;
      rpt << "dd[1].norm = " << cdd[1].norm() << endl;
      rpt << "\n=======   y [ucb, upb] =============\n";
      rpt << "bd[0].norm = " << cbd[0].norm() << endl;
      rpt << "bd[1].norm = " << cbd[1].norm() << endl;
      rpt << "\n=======A^Ty [csqb, ucb, upb] =======\n";
      rpt << "bm[0].norm = " << cbm[0].norm() << endl;
      rpt << "bm[1].norm = " << cbm[1].norm() << endl;
      rpt << "bm[2].norm = " << cbm[2].norm() << endl;
      
      
      rpt << "<csqd, csqb> = " << cdm[0].inner(cbm[0]) << endl;
      // compute norms and inner products
      ireal axnm, ynm, axy, xaty;
      axnm = dd.norm();
      ynm  = bd.norm();
      axy  = dd.inner(bd);
      xaty = bm.inner(dm);
      rpt << "<Ax,    y> = " << axy << endl;
      rpt << "< x, A^Ty> = " << xaty << endl;
      rpt << " |Ax| * |y| = " << axnm * ynm << endl;
      rpt << " adjvalue = " << abs(axy - xaty)/(axnm*ynm) << endl;

      EXPECT_GT(1.0e-6,abs(axy - xaty)/(axnm*ynm));
      
      //      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }  

  TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_3step_ord1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_3step_ord1\n\n";
      
      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie_3step.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;
      
      bool dryrun=false;
      //      ofstream drystr("op_adj_onestep_ord1");
      
      IWaveOp op(*pars,stream,dryrun,rpt,rpt);
      
      Vector<ireal> m(op.getDomain());
      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);            
      Components<ireal> cm(m);
      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);      
      cm[1].eval(rnd);
      cm[2].zero();

      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();

      bool res = AdjointTest(lop,rnd,rpt);

      EXPECT_EQ(res,true);

      //      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }

  TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_100step_ord1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_100step_ord1\n\n";
      
      // fake command line environment
      int argc = 2;
      char * argvv = new char[128];
      char ** argv = new char*[2];
      argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
      strcpy(argv[1],"par=movie_100step.par");
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      delete [] argv;
      delete [] argvv;
      
      bool dryrun=false;
      //      ofstream drystr("op_adj_onestep_ord1");
      
      IWaveOp op(*pars,stream,dryrun,rpt,rpt);
      
      Vector<ireal> m(op.getDomain());
      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);            
      Components<ireal> cm(m);
      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);      
      cm[1].eval(rnd);
      cm[2].zero();

      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();

      bool res = AdjointTest(lop,rnd,rpt);

      EXPECT_EQ(res,true);

      //      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }
  
  TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_1step_fdord4_ord1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_1step_fdord4_ord1\n\n";

      // fake command line environment
      std::vector<std::string> args(0);
      args.push_back("par=movie_1step.par");
      args.push_back("order=2");

      int argc;
      char ** argv;
      constructargs(args,&argc,&argv);
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      destroyargs(argc,&argv);

      bool dryrun=false;
      //      ofstream drystr("op_adj_onestep_ord1");
      
      IWaveOp op(*pars,stream,dryrun,rpt,rpt);
      
      Vector<ireal> m(op.getDomain());
      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);            
      Components<ireal> cm(m);
      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);      
      cm[1].eval(rnd);
      cm[2].zero();

      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();

      bool res = AdjointTest(lop,rnd,rpt);

      EXPECT_EQ(res,true);

      //      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }
  
  TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_3step_fdord4_ord1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_3step_fdord4_ord1\n\n";

      // fake command line environment
      std::vector<std::string> args(0);
      args.push_back("par=movie_3step.par");
      args.push_back("order=2");

      int argc;
      char ** argv;
      constructargs(args,&argc,&argv);
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      destroyargs(argc,&argv);

      bool dryrun=false;
      //      ofstream drystr("op_adj_onestep_ord1");
      
      IWaveOp op(*pars,stream,dryrun,rpt,rpt);
      
      Vector<ireal> m(op.getDomain());
      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);            
      Components<ireal> cm(m);
      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);      
      cm[1].eval(rnd);
      cm[2].zero();

      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();

      bool res = AdjointTest(lop,rnd,rpt);

      EXPECT_EQ(res,true);

      //      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }
  
  TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_100step_fdord4_ord1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_100step_fdord4_ord1\n\n";

      // fake command line environment
      std::vector<std::string> args(0);
      args.push_back("par=movie_100step.par");
      args.push_back("order=2");

      int argc;
      char ** argv;
      constructargs(args,&argc,&argv);
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      destroyargs(argc,&argv);

      bool dryrun=false;
      //      ofstream drystr("op_adj_onestep_ord1");
      
      IWaveOp op(*pars,stream,dryrun,rpt,rpt);
      
      Vector<ireal> m(op.getDomain());
      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);            
      Components<ireal> cm(m);
      AssignFilename mfn0("csq.rsf");
      cm[0].eval(mfn0);      
      cm[1].eval(rnd);
      cm[2].zero();

      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();

      bool res = AdjointTest(lop,rnd,rpt);

      EXPECT_EQ(res,true);

      //      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }
    
    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_1step_fdord8_ord1) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_1step_fdord8_ord1\n\n";
            
            // fake command line environment
            std::vector<std::string> args(0);
            args.push_back("par=movie_1step.par");
            args.push_back("order=4");
            
            int argc;
            char ** argv;
            constructargs(args,&argc,&argv);
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            destroyargs(argc,&argv);
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
            Components<ireal> cm(m);
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            cm[1].eval(rnd);
            cm[2].zero();
            
            OperatorEvaluation<ireal> opeval(op,m);
            LinearOp<ireal> const & lop = opeval.getDeriv();
            
            bool res = AdjointTest(lop,rnd,rpt);
            
            EXPECT_EQ(res,true);
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }
    
    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_3step_fdord8_ord1) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_3step_fdord8_ord1\n\n";
            
            // fake command line environment
            std::vector<std::string> args(0);
            args.push_back("par=movie_3step.par");
            args.push_back("order=4");
            
            int argc;
            char ** argv;
            constructargs(args,&argc,&argv);
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            destroyargs(argc,&argv);
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
            Components<ireal> cm(m);
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            cm[1].eval(rnd);
            cm[2].zero();
            
            OperatorEvaluation<ireal> opeval(op,m);
            LinearOp<ireal> const & lop = opeval.getDeriv();
            
            bool res = AdjointTest(lop,rnd,rpt);
            
            EXPECT_EQ(res,true);
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }
    
    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_100step_fdord8_ord1) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_100step_fdord8_ord1\n\n";
            
            // fake command line environment
            std::vector<std::string> args(0);
            args.push_back("par=movie_100step.par");
            args.push_back("order=4");
            
            int argc;
            char ** argv;
            constructargs(args,&argc,&argv);
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            destroyargs(argc,&argv);
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
            Components<ireal> cm(m);
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            cm[1].eval(rnd);
            cm[2].zero();
            
            OperatorEvaluation<ireal> opeval(op,m);
            LinearOp<ireal> const & lop = opeval.getDeriv();
            
            bool res = AdjointTest(lop,rnd,rpt);
            
            EXPECT_EQ(res,true);
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }
*/
    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uonly_onestep_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uonly_onestep_ord2\n\n";
            // fake command line environment
            int argc = 2;
            char * argvv = new char[128];
            char ** argv = new char*[2];
            argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
            strcpy(argv[1],"par=movie_1step.par");
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            delete [] argv;
            delete [] argvv;
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            
            Vector<ireal> dm(op.getDomain());
            Vector<ireal> dm0(op.getDomain());
            Vector<ireal> bm(op.getDomain());
            
            Vector<ireal> dd(op.getRange());
            Vector<ireal> bd(op.getRange());
            
            //      int seed = 19490615;
            int seed = getpid();
            RVL::RVLRandomize<float> rnd(seed,-1.0f,1.0f);
            Components<ireal> cm(m);
            Components<ireal> cdm(dm);
            Components<ireal> cdd(dd);
            Components<ireal> cbm(bm);
            Components<ireal> cbd(bd);
            
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            cm[1].eval(rnd);
            cm[2].zero();
            dm0.eval(rnd);
            
            rpt << "m.norm = " << m.norm() << endl;
            rpt << "m[0].norm = " << cm[0].norm() << endl;
            rpt << "m[1].norm = " << cm[1].norm() << endl;
            rpt << "m[2].norm = " << cm[2].norm() << endl;
            
            cdm[0].zero();
            cdm[1].eval(rnd);
            cdm[2].eval(rnd);
            cbd[0].eval(rnd);
            cbd[1].eval(rnd);
            
            dd.zero();
            bm.zero();
            
            // apply Op and AdjOp
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm0);
            lbl.applyOp(dm,dd);
            lbl.applyAdjOp(bd,bm);
            rpt << "\n=======   x ========================\n";
            rpt << "dm[0].norm = " << cdm[0].norm() << endl;
            rpt << "dm[1].norm = " << cdm[1].norm() << endl;
            rpt << "dm[2].norm = " << cdm[2].norm() << endl;
            rpt << "\n=======  Ax ========================\n";
            rpt << "dd[0].norm = " << cdd[0].norm() << endl;
            rpt << "dd[1].norm = " << cdd[1].norm() << endl;
            rpt << "\n=======   y [ucdb, updb] =============\n";
            rpt << "bd[0].norm = " << cbd[0].norm() << endl;
            rpt << "bd[1].norm = " << cbd[1].norm() << endl;
            rpt << "\n=======A^Ty [csqb, ucdb, updb] =======\n";
            rpt << "bm[0].norm = " << cbm[0].norm() << endl;
            rpt << "bm[1].norm = " << cbm[1].norm() << endl;
            rpt << "bm[2].norm = " << cbm[2].norm() << endl;
            
            
            rpt << "<csqd, csqb> = " << cdm[0].inner(cbm[0]) << endl;
            // compute norms and inner products
            ireal axnm, ynm, axy, xaty;
            axnm = dd.norm();
            ynm  = bd.norm();
            axy  = dd.inner(bd);
            xaty = bm.inner(dm);
            rpt << "<Ax,    y> = " << axy << endl;
            rpt << "< x, A^Ty> = " << xaty << endl;
            rpt << " |Ax| * |y| = " << axnm * ynm << endl;
            rpt << " adjvalue = " << abs(axy - xaty)/(axnm*ynm) << endl;
            
            EXPECT_GT(1.0e-6,abs(axy - xaty)/(axnm*ynm));
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }

    TEST_F(ACDCauchyOpAdjTest, wetrun_op_conly_onestep_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_conly_onestep_ord2\n\n";
            
            // fake command line environment
            int argc = 2;
            char * argvv = new char[128];
            char ** argv = new char*[2];
            argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
            strcpy(argv[1],"par=movie_1step.par");
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            delete [] argv;
            delete [] argvv;
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            
            Vector<ireal> dm(op.getDomain());
            Vector<ireal> dm0(op.getDomain());
            Vector<ireal> bm(op.getDomain());
            
            Vector<ireal> dd(op.getRange());
            Vector<ireal> bd(op.getRange());
            
            //      int seed = 19490615;
            int seed = getpid();
            RVL::RVLRandomize<float> rnd(seed,-1.0f,1.0f);
            Components<ireal> cm(m);
            Components<ireal> cdm(dm);
            Components<ireal> cdm0(dm0);
            Components<ireal> cdd(dd);
            Components<ireal> cbm(bm);
            Components<ireal> cbd(bd);
            
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            //      AssignFilename ucfn("../uc_spike.rsf");
            cm[1].eval(rnd);
            cm[2].zero();
//            cdm0[0].eval(rnd);
//            cdm0[1].zero();
//            cdm0[2].zero();
            dm0.eval(rnd);
            rpt << "m.norm = " << m.norm() << endl;
            rpt << "m[0].norm = " << cm[0].norm() << endl;
            rpt << "m[1].norm = " << cm[1].norm() << endl;
            rpt << "m[2].norm = " << cm[2].norm() << endl;
 
            //      AssignFilename csqdfn("../csq_spike.rsf");
            cdm[0].eval(rnd);
            cdm[1].zero();
            cdm[2].zero();
            //      AssignFilename ucbdfn("../ucbd_spike.rsf");
            cbd[0].eval(rnd);
            cbd[1].zero();
            
            dd.zero();
            bm.zero();
            
            // apply Op and AdjOp
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm0);
            lbl.applyOp(dm,dd);
            lbl.applyAdjOp(bd,bm);
            rpt << "\n=======   x ========================\n";
            rpt << "dm[0].norm = " << cdm[0].norm() << endl;
            rpt << "dm[1].norm = " << cdm[1].norm() << endl;
            rpt << "dm[2].norm = " << cdm[2].norm() << endl;
            rpt << "\n=======  Ax ========================\n";
            rpt << "dd[0].norm = " << cdd[0].norm() << endl;
            rpt << "dd[1].norm = " << cdd[1].norm() << endl;
            rpt << "\n=======   y [ucb, upb] =============\n";
            rpt << "bd[0].norm = " << cbd[0].norm() << endl;
            rpt << "bd[1].norm = " << cbd[1].norm() << endl;
            rpt << "\n=======A^Ty [csqb, ucb, upb] =======\n";
            rpt << "bm[0].norm = " << cbm[0].norm() << endl;
            rpt << "bm[1].norm = " << cbm[1].norm() << endl;
            rpt << "bm[2].norm = " << cbm[2].norm() << endl;
            
            
            rpt << "<csqd, csqb> = " << cdm[0].inner(cbm[0]) << endl;
            // compute norms and inner products
            ireal axnm, ynm, axy, xaty;
            axnm = dd.norm();
            ynm  = bd.norm();
            axy  = dd.inner(bd);
            xaty = bm.inner(dm);
            rpt << "<Ax,    y> = " << axy << endl;
            rpt << "< x, A^Ty> = " << xaty << endl;
            rpt << " |Ax| * |y| = " << axnm * ynm << endl;
            rpt << " adjvalue = " << abs(axy - xaty)/(axnm*ynm) << endl;
            
            EXPECT_GT(1.0e-6,abs(axy - xaty)/(axnm*ynm));
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }
    
    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_onestep_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_onestep_ord2\n\n";
            
            // fake command line environment
            int argc = 2;
            char * argvv = new char[128];
            char ** argv = new char*[2];
            argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
            strcpy(argv[1],"par=movie_1step.par");
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            delete [] argv;
            delete [] argvv;
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            Vector<ireal> dm(op.getDomain());
            RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
            Components<ireal> cm(m);
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            cm[1].eval(rnd);
            cm[2].zero();
            dm.eval(rnd);
            
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm);

            bool res = AdjointTest(lbl,rnd,rpt);
            
            EXPECT_EQ(res,true);
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }

    TEST_F(ACDCauchyOpAdjTest, dryrun_op_3step_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, dryrun_op_3step_ord2\n\n";
            
            // fake command line environment
            int argc = 2;
            char * argvv = new char[128];
            char ** argv = new char*[2];
            argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
            strcpy(argv[1],"par=movie_3step.par");
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            delete [] argv;
            delete [] argvv;
            
            bool dryrun=true;
            ofstream drystr("dryrun_op_adj_3step_ord2");
            
            IWaveOp op(*pars,stream,dryrun,drystr,rpt);
            
            Vector<ireal> m(op.getDomain());
            
            Vector<ireal> dm(op.getDomain());
            Vector<ireal> dm0(op.getDomain());
            Vector<ireal> bm(op.getDomain());
            
            Vector<ireal> dd(op.getRange());
            Vector<ireal> bd(op.getRange());
            
            int seed = 19490615;
            RVL::RVLRandomize<float> rnd(seed,-1.0f,1.0f);
            Components<ireal> cm(m);
            Components<ireal> cdm(dm);
            Components<ireal> cdd(dd);
            Components<ireal> cbm(bm);
            Components<ireal> cbd(bd);
            
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            
            cm[1].eval(rnd);
            cm[2].zero();
            dm0.eval(rnd);
            cdm[0].eval(rnd);
            cdm[1].zero();
            cdm[2].zero();
            cbd[0].eval(rnd);
            cbd[1].zero();
            
            dd.zero();
            bm.zero();
            
            // apply Op and AdjOp
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm0);
            lbl.applyOp(dm,dd);
            lbl.applyAdjOp(bd,bm);
            
            drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }

    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uonly_3step_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uonly_3step_ord2\n\n";
            
            // fake command line environment
            int argc = 2;
            char * argvv = new char[128];
            char ** argv = new char*[2];
            argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
            strcpy(argv[1],"par=movie_3step.par");
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            delete [] argv;
            delete [] argvv;
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            
            Vector<ireal> dm(op.getDomain());
            Vector<ireal> dm0(op.getDomain());
            Vector<ireal> bm(op.getDomain());
            
            Vector<ireal> dd(op.getRange());
            Vector<ireal> bd(op.getRange());
            
            //      int seed = 19490615;
            int seed = getpid();
            RVL::RVLRandomize<float> rnd(seed,-1.0f,1.0f);
            Components<ireal> cm(m);
            Components<ireal> cdm(dm);
            Components<ireal> cdd(dd);
            Components<ireal> cbm(bm);
            Components<ireal> cbd(bd);
            
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            
            cm[1].eval(rnd);
            cm[2].zero();
            dm0.eval(rnd);
            rpt << "m.norm = " << m.norm() << endl;
            rpt << "m[0].norm = " << cm[0].norm() << endl;
            rpt << "m[1].norm = " << cm[1].norm() << endl;
            rpt << "m[2].norm = " << cm[2].norm() << endl;
            
            cdm[0].zero();
            cdm[1].eval(rnd);
            cdm[2].eval(rnd);
            cbd[0].eval(rnd);
            cbd[1].eval(rnd);
            
            dd.zero();
            bm.zero();
            
            // apply Op and AdjOp
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm0);
            lbl.applyOp(dm,dd);
            lbl.applyAdjOp(bd,bm);
            rpt << "\n=======   x ========================\n";
            rpt << "dm[0].norm = " << cdm[0].norm() << endl;
            rpt << "dm[1].norm = " << cdm[1].norm() << endl;
            rpt << "dm[2].norm = " << cdm[2].norm() << endl;
            rpt << "\n=======  Ax ========================\n";
            rpt << "dd[0].norm = " << cdd[0].norm() << endl;
            rpt << "dd[1].norm = " << cdd[1].norm() << endl;
            rpt << "\n=======   y [ucb, upb] =============\n";
            rpt << "bd[0].norm = " << cbd[0].norm() << endl;
            rpt << "bd[1].norm = " << cbd[1].norm() << endl;
            rpt << "\n=======A^Ty [csqb, ucb, upb] =======\n";
            rpt << "bm[0].norm = " << cbm[0].norm() << endl;
            rpt << "bm[1].norm = " << cbm[1].norm() << endl;
            rpt << "bm[2].norm = " << cbm[2].norm() << endl;
            
            
            rpt << "<csqd, csqb> = " << cdm[0].inner(cbm[0]) << endl;
            // compute norms and inner products
            ireal axnm, ynm, axy, xaty;
            axnm = dd.norm();
            ynm  = bd.norm();
            axy  = dd.inner(bd);
            xaty = bm.inner(dm);
            rpt << "<Ax,    y> = " << axy << endl;
            rpt << "< x, A^Ty> = " << xaty << endl;
            rpt << " |Ax| * |y| = " << axnm * ynm << endl;
            rpt << " adjvalue = " << abs(axy - xaty)/(axnm*ynm) << endl;
            
            EXPECT_GT(1.0e-6,abs(axy - xaty)/(axnm*ynm));
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }
    
    TEST_F(ACDCauchyOpAdjTest, wetrun_op_conly_3step_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_conly_3step_ord2\n\n";
            
            // fake command line environment
            int argc = 2;
            char * argvv = new char[128];
            char ** argv = new char*[2];
            argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
            strcpy(argv[1],"par=movie_3step.par");
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            delete [] argv;
            delete [] argvv;
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            
            Vector<ireal> dm(op.getDomain());
            Vector<ireal> dm0(op.getDomain());
            Vector<ireal> bm(op.getDomain());
            
            Vector<ireal> dd(op.getRange());
            Vector<ireal> bd(op.getRange());
            
            //      int seed = 19490615;
            int seed = getpid();
            RVL::RVLRandomize<float> rnd(seed,-1.0f,1.0f);
            Components<ireal> cm(m);
            Components<ireal> cdm(dm);
            Components<ireal> cdd(dd);
            Components<ireal> cbm(bm);
            Components<ireal> cbd(bd);
            
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            //      AssignFilename ucfn("../uc_spike.rsf");
            cm[1].eval(rnd);
            cm[2].zero();
            dm0.eval(rnd);
            
            rpt << "m.norm = " << m.norm() << endl;
            rpt << "m[0].norm = " << cm[0].norm() << endl;
            rpt << "m[1].norm = " << cm[1].norm() << endl;
            rpt << "m[2].norm = " << cm[2].norm() << endl;
            
            AssignFilename csqdfn("../csq_spike.rsf");
            cdm[0].eval(rnd);
            //      cdm[0].eval(csqdfn);
            cdm[1].zero();
            cdm[2].zero();
            //      AssignFilename ucbdfn("../ucbd_spike.rsf");
            cbd[0].eval(rnd);
            cbd[1].zero();
            
            dd.zero();
            bm.zero();
            
            // apply Op and AdjOp
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm0);
            lbl.applyOp(dm,dd);
            lbl.applyAdjOp(bd,bm);
            rpt << "\n=======   x ========================\n";
            rpt << "dm[0].norm = " << cdm[0].norm() << endl;
            rpt << "dm[1].norm = " << cdm[1].norm() << endl;
            rpt << "dm[2].norm = " << cdm[2].norm() << endl;
            rpt << "\n=======  Ax ========================\n";
            rpt << "dd[0].norm = " << cdd[0].norm() << endl;
            rpt << "dd[1].norm = " << cdd[1].norm() << endl;
            rpt << "\n=======   y [ucb, upb] =============\n";
            rpt << "bd[0].norm = " << cbd[0].norm() << endl;
            rpt << "bd[1].norm = " << cbd[1].norm() << endl;
            rpt << "\n=======A^Ty [csqb, ucb, upb] =======\n";
            rpt << "bm[0].norm = " << cbm[0].norm() << endl;
            rpt << "bm[1].norm = " << cbm[1].norm() << endl;
            rpt << "bm[2].norm = " << cbm[2].norm() << endl;
            
            
            rpt << "<csqd, csqb> = " << cdm[0].inner(cbm[0]) << endl;
            // compute norms and inner products
            ireal axnm, ynm, axy, xaty;
            axnm = dd.norm();
            ynm  = bd.norm();
            axy  = dd.inner(bd);
            xaty = bm.inner(dm);
            rpt << "<Ax,    y> = " << axy << endl;
            rpt << "< x, A^Ty> = " << xaty << endl;
            rpt << " |Ax| * |y| = " << axnm * ynm << endl;
            rpt << " adjvalue = " << abs(axy - xaty)/(axnm*ynm) << endl;
            
            EXPECT_GT(1.0e-6,abs(axy - xaty)/(axnm*ynm));
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }
    
    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_3step_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_3step_ord2\n\n";
            
            // fake command line environment
            int argc = 2;
            char * argvv = new char[128];
            char ** argv = new char*[2];
            argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
            strcpy(argv[1],"par=movie_3step.par");
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            delete [] argv;
            delete [] argvv;
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            Vector<ireal> dm0(op.getDomain());
            RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
            Components<ireal> cm(m);
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            cm[1].eval(rnd);
            cm[2].zero();
            dm0.eval(rnd);
            
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm0);
            bool res = AdjointTest(lbl,rnd,rpt);
            
            EXPECT_EQ(res,true);
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }
    
    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_100step_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_100step_ord2\n\n";
            
            // fake command line environment
            int argc = 2;
            char * argvv = new char[128];
            char ** argv = new char*[2];
            argv[0]=&(argvv[0]); argv[1]=&(argvv[65]);
            strcpy(argv[1],"par=movie_100step.par");
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            delete [] argv;
            delete [] argvv;
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            Vector<ireal> dm0(op.getDomain());
            RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
            Components<ireal> cm(m);
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            cm[1].eval(rnd);
            cm[2].zero();
            dm0.eval(rnd);
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm0);
            
            bool res = AdjointTest(lbl,rnd,rpt);
            
            EXPECT_EQ(res,true);
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }

    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_1step_fdord4_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_1step_fdord4_ord2\n\n";
            
            // fake command line environment
            std::vector<std::string> args(0);
            args.push_back("par=movie_1step.par");
            args.push_back("order=2");
            
            int argc;
            char ** argv;
            constructargs(args,&argc,&argv);
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            destroyargs(argc,&argv);
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            Vector<ireal> dm0(op.getDomain());
            RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
            Components<ireal> cm(m);
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            cm[1].eval(rnd);
            cm[2].zero();
            dm0.eval(rnd);
            
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm0);
            bool res = AdjointTest(lbl,rnd,rpt);
            
            EXPECT_EQ(res,true);
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }
    
    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_3step_fdord4_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_3step_fdord4_ord2\n\n";
            
            // fake command line environment
            std::vector<std::string> args(0);
            args.push_back("par=movie_3step.par");
            args.push_back("order=2");
            
            int argc;
            char ** argv;
            constructargs(args,&argc,&argv);
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            destroyargs(argc,&argv);
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            Vector<ireal> dm0(op.getDomain());
            RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
            Components<ireal> cm(m);
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            cm[1].eval(rnd);
            cm[2].zero();
            dm0.eval(rnd);
            
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm0);
            
            bool res = AdjointTest(lbl,rnd,rpt);
            
            EXPECT_EQ(res,true);
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }
    
    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_100step_fdord4_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_100step_fdord4_ord2\n\n";
            
            // fake command line environment
            std::vector<std::string> args(0);
            args.push_back("par=movie_100step.par");
            args.push_back("order=2");
            
            int argc;
            char ** argv;
            constructargs(args,&argc,&argv);
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            destroyargs(argc,&argv);
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            Vector<ireal> dm0(op.getDomain());
            RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
            Components<ireal> cm(m);
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            cm[1].eval(rnd);
            cm[2].zero();
            dm0.eval(rnd);
            
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm0);
            
            bool res = AdjointTest(lbl,rnd,rpt);
            
            EXPECT_EQ(res,true);
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }
    
    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_1step_fdord8_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_1step_fdord4_ord2\n\n";
            
            // fake command line environment
            std::vector<std::string> args(0);
            args.push_back("par=movie_1step.par");
            args.push_back("order=4");
            
            int argc;
            char ** argv;
            constructargs(args,&argc,&argv);
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            destroyargs(argc,&argv);
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            Vector<ireal> dm0(op.getDomain());
            RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
            Components<ireal> cm(m);
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            cm[1].eval(rnd);
            cm[2].zero();
            dm0.eval(rnd);
            
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm0);
            bool res = AdjointTest(lbl,rnd,rpt);
            
            EXPECT_EQ(res,true);
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }
    
    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_3step_fdord8_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_3step_fdord8_ord2\n\n";
            
            // fake command line environment
            std::vector<std::string> args(0);
            args.push_back("par=movie_3step.par");
            args.push_back("order=4");
            
            int argc;
            char ** argv;
            constructargs(args,&argc,&argv);
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            destroyargs(argc,&argv);
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            Vector<ireal> dm0(op.getDomain());
            RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
            Components<ireal> cm(m);
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            cm[1].eval(rnd);
            cm[2].zero();
            dm0.eval(rnd);
            
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm0);
            
            bool res = AdjointTest(lbl,rnd,rpt);
            
            EXPECT_EQ(res,true);
            
            //      drystr.close();
            ps_delete(&pars);
            fclose(stream);
        }
        catch (RVLException & e) {
            e.write(cerr);
            exit(1);
        }
    }
    
    TEST_F(ACDCauchyOpAdjTest, wetrun_op_uandc_100step_fdord8_ord2) {
        try {
            rpt<<"\n======================================================\n";
            rpt<<"ACDCauchyOpAdjTest, wetrun_op_uandc_100step_fdord8_ord2\n\n";
            
            // fake command line environment
            std::vector<std::string> args(0);
            args.push_back("par=movie_100step.par");
            args.push_back("order=4");
            
            int argc;
            char ** argv;
            constructargs(args,&argc,&argv);
            PARARRAY * pars = NULL;
            FILE * stream = NULL;
            IWaveEnvironment(argc, argv, 0, &pars, &stream);
            destroyargs(argc,&argv);
            
            bool dryrun=false;
            //      ofstream drystr("op_adj_onestep_ord1");
            
            IWaveOp op(*pars,stream,dryrun,rpt,rpt);
            
            Vector<ireal> m(op.getDomain());
            Vector<ireal> dm0(op.getDomain());
            RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);
            Components<ireal> cm(m);
            AssignFilename mfn0("csq.rsf");
            cm[0].eval(mfn0);
            cm[1].eval(rnd);
            cm[2].zero();
            dm0.eval(rnd);
            
            OperatorEvaluation<ireal> opeval(op,m);
            SymmetricBilinearOp<ireal> const & blop = opeval.getDeriv2();
            LinearBilinearOp<float> lbl(blop,dm0);
            
            bool res = AdjointTest(lbl,rnd,rpt);
            
            EXPECT_EQ(res,true);
            
            //      drystr.close();
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
