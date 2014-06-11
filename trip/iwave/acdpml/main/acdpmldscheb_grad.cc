#include "acdpml_defn.hh"
#include "grid.h"
#include "gridpp.hh"
#include "gridops.hh"
#include "istate.hh"
#include "iwop.hh"
#include "functions.hh"
#include "op.hh"
#include "blockop.hh"
#include "chebalg.hh"
#include "LBFGSBT.hh"
#include "LinFitLS.hh"
#include "acdpmlds_grad_selfdoc.hh"
#include <par.h>
#include "gradtest.hh"
#include "segyops.hh"

//#define GTEST_VERBOSE
IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csqext",    0, true,  true },
  {"data",   1, false, true },
  {"source", 1, true,  false},
  {"",       0, false, false}
};

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
using RVL::ScaleOpFwd;
using RVL::TensorOp;
using RVL::GradientTest;
using RVLUmin::ChebPolicy;
using RVLUmin::ChebPolicyData;
using RVLUmin::LBFGSBT;
using RVLUmin::LinFitLS;
using RVLUmin::ChebAlg;
using TSOpt::IWaveEnvironment;
using TSOpt::IWaveTree;
using TSOpt::IWaveSampler;
using TSOpt::IWaveSim;
using TSOpt::TASK_RELN;
using TSOpt::IOTask;
using TSOpt::IWaveOp;
using TSOpt::SEGYLinMute;
#ifdef IWAVE_USE_MPI
using TSOpt::MPIGridSpace;
using TSOpt::MPISEGYSpace;
#else
using TSOpt::GridSpace;
using TSOpt::SEGYSpace;
#endif
using TSOpt::GridExtendOp;
using TSOpt::GridDerivOp;

int xargc;
char **xargv;

/** this version requires that two model files be present:
    csqext - extended version of csq, with [d=spatial dimn]
    gdim=d+1
    dim=d
    n[d+1]=number of shot records, 
    d[d+1]=1
    o[d+1]=0
    id[d+1]=d+1
    csq - physical space version of csq, with dim=gdim=d
    Uses these keywords for geometry only - data content irrelevant

    [intended near-term update: require only physical space material
    params, generate extended version internally, so automatic 
    compatibility with data]

*/
int main(int argc, char ** argv) {

  try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif

    PARARRAY * pars = NULL;
    FILE * stream = NULL;
    IWaveEnvironment(argc, argv, 0, &pars, &stream);

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }

#ifdef IWAVE_USE_MPI
    if (retrieveGroupID() == MPI_UNDEFINED) {
      fprintf(stream,"NOTE: finalize MPI, cleanup, exit\n");
      fflush(stream);
    }
    else {
#endif

      // the Op
      IWaveOp iwop(*pars,stream);
      SEGYLinMute mute(valparse<float>(*pars,"mute_slope",0.0f),
                       valparse<float>(*pars,"mute_zotime",0.0f),
                       valparse<float>(*pars,"mute_width",0.0f));
        
      LinearOpFO<float> muteop(iwop.getRange(),iwop.getRange(),mute,mute);
      OpComp<float> op(iwop,muteop);
    
      /* generate physical model space */
#ifdef IWAVE_USE_MPI
      MPIGridSpace csqsp(valparse<std::string>(*pars,"csq"),"notype",true);
#else
      GridSpace csqsp(valparse<std::string>(*pars,"csq"),"notype",true);
#endif
      // make it a product, so it's compatible with domain of op
      StdProductSpace<ireal> dom(csqsp);

      // vel-squared model!
      Vector<ireal> m(dom);
      Vector<ireal> m_end(dom);
      Vector<ireal> dm(dom);
      Vector<ireal> dd(op.getRange());
      Vector<ireal> mdd(op.getRange());

      // read in start point and end point
      AssignFilename mfn(valparse<std::string>(*pars,"csq"));
      Components<ireal> cm(m);
      cm[0].eval(mfn);

      AssignFilename mfn_end(valparse<std::string>(*pars,"csq_end"));
      Components<ireal> cm_end(m_end);
      cm_end[0].eval(mfn_end);
        
      Components<ireal> cdm(dm);
      dm.copy(m_end);
      dm.linComb(-1.0f,m);
    
      AssignFilename ddfn(valparse<std::string>(*pars,"data"));
      Components<ireal> cdd(dd);
      cdd[0].eval(ddfn);

      std::string mddnm = valparse<std::string>(*pars,"datamut","");
      if (mddnm.size()>0) {
        AssignFilename mddfn(mddnm);
        mdd.eval(mddfn);
      }
      muteop.applyOp(dd,mdd);

      ChebPolicyData<float> pd(valparse<float>(*pars,"gamma",0.04f),
                               valparse<float>(*pars,"epsilon",0.1),
                               valparse<float>(*pars,"alpha",1.01),
                               valparse<float>(*pars,"lbd_est",1.01),
                               valparse<int>(*pars,"MaxIter",10),true);
    
      /* output stream */
      std::stringstream res,strgrad;
      res<<scientific;
      strgrad<<scientific;
    
      // choice of GridDerivOp for semblance op is placeholder - extended axis is dim-th, with id=dim+1
      // note default weight of zero!!!
      // Note added 10.03.14: getGrid is not usably implemented for MPIGridSpace at this time, so 
      // must fish the derivative index out by hand and bcast
      // THIS SHOULD ALL BE FIXED! (1) getGrid properly implemented in parallel, (2) proper
      // defn of DSOp to include internal extd axis case (space/time shift)
      int dsdir = INT_MAX;
      if (retrieveGlobalRank()==0) dsdir=csqsp.getGrid().dim;
#ifdef IWAVE_USE_MPI
      if (MPI_Bcast(&dsdir,1,MPI_INT,0,retrieveGlobalComm())) {
	RVLException e;
	e<<"Error: acdiva, rank="<<retrieveGlobalRank()<<"\n";
	e<<"  failed to bcast dsdir\n";
	throw e;
      }
#endif
      GridDerivOp dsop(op.getDomain(),dsdir,valparse<float>(*pars,"DSWeight",0.0f));
      TensorOp<float> top(op,dsop);
      // create RHS of block system
      Vector<float> td(top.getRange());
      Components<float> ctd(td);
      ctd[0].copy(mdd);
      ctd[1].zero();

      // choice of preop is placeholder
      ScaleOpFwd<float> preop(top.getDomain(),1.0f);

      LinFitLS<float, ChebPolicy<float>, ChebPolicyData<float> > f(top,preop,td,pd,valparse<bool>(*pars,"refine",false),res);
      GridExtendOp g(dom,op.getDomain());
      FcnlOpComp<float> gf(f,g);

      // lower, upper bds for csq
      Vector<float> lb(gf.getDomain());
      Vector<float> ub(gf.getDomain());
      float lbcsq=valparse<float>(*pars,"cmin");
      lbcsq=lbcsq*lbcsq;
      RVLAssignConst<float> asl(lbcsq);
      lb.eval(asl);
      float ubcsq=valparse<float>(*pars,"cmax");
      ubcsq=ubcsq*ubcsq;
      RVLAssignConst<float> asu(ubcsq);
      ub.eval(asu);
      RVLMin<float> mn;
#ifdef IWAVE_USE_MPI
      MPISerialFunctionObjectRedn<float,float> mpimn(mn);
      ULBoundsTest<float> ultest(lb,ub,mpimn);
#else
      ULBoundsTest<float> ultest(lb,ub,mn);
#endif

      FunctionalBd<float> fbd(gf,ultest);

      // read in scan test parameters
      int  nhval=valparse<int>(*pars,"Nhval",11);
      float hmin=valparse<float>(*pars,"HMin",0.1f);
      float hmax=valparse<float>(*pars,"HMax",1.0f);
        
      /* scan */
      GradientTest(fbd,m,dm,strgrad,nhval,hmin,hmax);

      if (retrieveRank() == 0) {
	std::string outfile = valparse<std::string>(*pars,"outfile","");
	if (outfile.size()>0) {
	  ofstream outf(outfile.c_str());
          outf<<strgrad.str();
          outf<<"\n ================================================= \n";
	  outf<<res.str();
	  outf.close();
	}
	else {
          cout<<strgrad.str();
          cout<<"\n ================================================= \n";
	  cout<<res.str();
	}
      }
    
#ifdef IWAVE_USE_MPI
    }
    MPI_Finalize();
#endif
  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
  
}
