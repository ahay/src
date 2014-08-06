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
  {"data",   1, false, true },
  {"source", 1, true,  false},
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
    if (args.size() < 1) {
      RVLException e;
      e<<"Error: constructargs\n";
      e<<"  input std::vector must have length > 0\n";
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

  class ACDTraceOpAdjTest : public ::testing::Test {
  public:

    IWaveInfo ic;
    ofstream rpt;

    ACDTraceOpAdjTest(): ic() {

      rpt.open("details.txt",ios_base::app);

      ofstream qs("trace.par");
      qs<<"INPUT DATA FOR IWAVE\n";
      qs<<"------------------------------------------------------------------------\n";
      qs<<"FD:\n";
      qs<<"\n";
      qs<<"         order = 1           scheme half-order\n";
      qs<<"           cfl = 0.5        cfl number - frac of max stable\n";
      qs<<"          cmin = 1.0         min velocity - checked\n";
      qs<<"          cmax = 5.0         max velocity - checked\n";
      qs<<"\n";
      qs<<"------------------------------------------------------------------------\n";
      qs<<"Model info:\n";
      qs<<"\n";
      qs<<"           csq = ../csq.rsf\n";
      qs<<"          data = ../data.su\n";
      qs<<"        source = ../source.su\n";
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
      qs<<"    dump_term = 1           dump terminator data\n";
      qs<<"    dump_pars = 1           print parameter table in IWaveOp\n";
      qs<<"   dump_steps = 1           print major steps in IWaveOp\n";
      qs.flush();
      qs.close();
    }
  };

  TEST_F(ACDTraceOpAdjTest, dryrun_adjtest_fdord2_deriv1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDTraceOpAdjTest, dryrun_adjtest_fdord2_deriv1\n\n";

      // fake command line environment
      std::vector<std::string> args(0);
      args.push_back("par=trace.par");

      int argc;
      char ** argv;
      constructargs(args,&argc,&argv);
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      destroyargs(argc,&argv);

      bool dryrun=true;
      ofstream drystr("dryrun_fdord2_deriv1");
      
      IWaveOp op(*pars,stream,dryrun,drystr,rpt);
      
      Vector<ireal> m(op.getDomain());
      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);            
      Components<ireal> cm(m);
      AssignFilename mfn0("../csq.rsf");
      cm[0].eval(mfn0);      

      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();

      bool res = AdjointTest(lop,rnd,rpt);
      // result will be false because computations will not be done!!
      EXPECT_EQ(res,false);

      drystr.close();
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }
 
  TEST_F(ACDTraceOpAdjTest, byhand_adjtest_cv_fdord2_deriv1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDTraceOpAdjTest, byhand_adjtest_fdord2_deriv1\n\n";

      // fake command line environment
      std::vector<std::string> args(0);
      args.push_back("par=trace.par");

      int argc;
      char ** argv;
      constructargs(args,&argc,&argv);
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      destroyargs(argc,&argv);

      bool dryrun=false;
      
      IWaveOp op(*pars,stream,dryrun,rpt,rpt);
      
      Vector<ireal> m(op.getDomain());
      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);            
      Components<ireal> cm(m);
      AssignFilename mfn0("../csq.rsf");
      cm[0].eval(mfn0);      

      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();

      Vector<ireal> xin(op.getDomain());
      Vector<ireal> xout(op.getDomain());
      Vector<ireal> yin(op.getRange());
      Vector<ireal> yout(op.getRange());
      AssignFilename xinfn("../csq_spike.rsf");
      xin.eval(xinfn);
      AssignFilename youtfn("fwd_data.su");
      yout.eval(youtfn);
      yout.zero();
      xout.zero();
      yin.eval(rnd);

      lop.applyOp(xin,yout);
      lop.applyAdjOp(yin,xout);
      ireal nout=xout.norm();
      ireal nin=yin.norm();
      rpt<<"xin.norm="<<xin.norm()<<endl;
      rpt<<"yout.norm="<<yout.norm()<<endl;
      rpt<<"yin.norm="<<nin<<endl;
      rpt<<"xout.norm="<<nout<<endl;
      ireal ipout=xout.inner(xin);
      ireal ipin=yin.inner(yout);
      rpt<<"<xout,xin>="<<ipout<<endl;
      rpt<<"<yout,yin>="<<ipin<<endl;
      rpt<<"rel error = "<< (abs(ipout-ipin))/(nin*nout)<<endl;
      ps_delete(&pars);
      fclose(stream);
    }
    catch (RVLException & e) {
      e.write(cerr);
      exit(1);
    }
  }
  
  TEST_F(ACDTraceOpAdjTest, wetrun_adjtest_cv_fdord2_deriv1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDTraceOpAdjTest, wetrun_adjtest_cv_fdord2_deriv1\n\n";

      // fake command line environment
      std::vector<std::string> args(0);
      args.push_back("par=trace.par");

      int argc;
      char ** argv;
      constructargs(args,&argc,&argv);
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      destroyargs(argc,&argv);

      bool dryrun=false;
      //      ofstream drystr("dryrun_fdord2_deriv1");
      
      IWaveOp op(*pars,stream,dryrun,rpt,rpt);
      
      Vector<ireal> m(op.getDomain());
      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);            
      Components<ireal> cm(m);
      AssignFilename mfn0("../csq.rsf");
      cm[0].eval(mfn0);      

      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();

      bool res = AdjointTest(lop,rnd,rpt);
      // result will be false because computations will not be done!!
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
 
  TEST_F(ACDTraceOpAdjTest, wetrun_adjtest_dome_fdord2_deriv1) {
    try {
      rpt<<"\n======================================================\n";
      rpt<<"ACDTraceOpAdjTest, wetrun_adjtest_dome_fdord2_deriv1\n\n";

      // fake command line environment
      std::vector<std::string> args(0);
      args.push_back("par=trace.par");

      int argc;
      char ** argv;
      constructargs(args,&argc,&argv);
      PARARRAY * pars = NULL;
      FILE * stream = NULL;
      IWaveEnvironment(argc, argv, 0, &pars, &stream);
      destroyargs(argc,&argv);

      bool dryrun=false;
      //      ofstream drystr("dryrun_fdord2_deriv1");
      
      IWaveOp op(*pars,stream,dryrun,rpt,rpt);
      
      Vector<ireal> m(op.getDomain());
      RVL::RVLRandomize<float> rnd(getpid(),-1.0f,1.0f);            
      Components<ireal> cm(m);
      AssignFilename mfn0("../domecsq20m.rsf");
      cm[0].eval(mfn0);      

      OperatorEvaluation<ireal> opeval(op,m);
      LinearOp<ireal> const & lop = opeval.getDeriv();

      bool res = AdjointTest(lop,rnd,rpt);
      // result will be false because computations will not be done!!
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
