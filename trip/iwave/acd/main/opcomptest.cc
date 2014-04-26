#include "acd_defn.hh"
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
#include "acdds_grad_selfdoc.hh"
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
using RVL::LinearOpFO;
using RVL::OpComp;
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
      Vector<ireal> m(op.getDomain());
      Vector<ireal> dm0(op.getDomain());
      Vector<ireal> dm1(op.getDomain());
      Vector<ireal> bm0(op.getDomain());
      Vector<ireal> bm1(op.getDomain());
        
      Vector<ireal> dd(op.getRange());
      Vector<ireal> d0(op.getRange());
      Vector<ireal> d1(op.getRange());

      // read in start point and end point
      AssignFilename mfn(valparse<std::string>(*pars,"csqsm"));
      Components<ireal> cm(m);
      cm[0].eval(mfn);

      AssignFilename dmfn0(valparse<std::string>(*pars,"csq_d1"));
      Components<ireal> cdm0(dm0);
      cdm0[0].eval(dmfn0);
        
      AssignFilename dmfn1(valparse<std::string>(*pars,"csq_d2"));
      Components<ireal> cdm1(dm1);
      cdm1[0].eval(dmfn1);
        
      AssignFilename bmfn0(valparse<std::string>(*pars,"csq_b1"));
      Components<ireal> cbm0(bm0);
      cbm0[0].eval(bmfn0);
        
      AssignFilename bmfn1(valparse<std::string>(*pars,"csq_b2"));
      Components<ireal> cbm1(bm1);
      cbm1[0].eval(bmfn1);
      
//      AssignFilename dmfn(valparse<std::string>(*pars,"reflectivity"));
//      Components<ireal> cdm(dm);
//      cdm[0].eval(dmfn);
//      dm.zero();
    
      AssignFilename ddfn(valparse<std::string>(*pars,"data"));
      Components<ireal> cdd(dd);
      cdd[0].eval(ddfn);
        
      AssignFilename d0fn(valparse<std::string>(*pars,"data1"));
      Components<ireal> cd0(d0);
      cd0[0].eval(d0fn);
        
      AssignFilename d1fn(valparse<std::string>(*pars,"data2"));
      Components<ireal> cd1(d1);
      cd1[0].eval(d1fn);
    
      /* output stream */
      std::stringstream res,strgrad;
      res<<scientific;
      //strgrad<<scientific;

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
	e<<"Error: acddscheb_grad, rank="<<retrieveGlobalRank()<<"\n";
	e<<"  failed to bcast dsdir\n";
	throw e;
      }
#endif

        OperatorEvaluation<float> opeval(op, m);
        dd.copy(opeval.getValue());

        opeval.getDeriv().applyOp(dm0,d0);
        opeval.getDeriv().applyAdjOp(d0,bm0);
        float iprng = d0.normsq();
        float ipdom = dm0.inner(bm0);
        float denom = iprng;
        res<< "\n AdjontTest: OpComp(iwop,muteop)\n";
        res<< "\n Deriv=1: \n";
        res<< "<Ax,y>  = " << iprng << endl;
        res<< "<x,A'y> = " << ipdom << endl;
        res<< "|<Ax,y>-<x,A'y>|/|Ax||y| = " << fabs(iprng - ipdom)/denom << endl;
        
        opeval.getDeriv2().applyOp(dm0,dm1,d1);
        opeval.getDeriv2().applyAdjOp(dm0,d1,bm1);
        iprng = d1.normsq();
        ipdom = dm1.inner(bm1);
        denom = iprng;
        res<< "\n Deriv=2: \n";
        res<< "<Ax,y>  = " << iprng << endl;
        res<< "<x,A'y> = " << ipdom << endl;
        res<< "|<Ax,y>-<x,A'y>|/|Ax||y| = " << fabs(iprng - ipdom)/denom << endl;
        
/*        OperatorEvaluation<float> iwopeval(iwop, m);

        iwopeval.getDeriv().applyOp(dm0,d0);
        iwopeval.getDeriv().applyAdjOp(d0,bm0);
        iprng = d0.normsq();
        ipdom = dm0.inner(bm0);
        denom = iprng;
        res<< "\n AdjontTest: IWaveOp\n";
        res<< "\n Deriv=1: \n";
        res<< "<Ax,y>  = " << iprng << endl;
        res<< "<x,A'y> = " << ipdom << endl;
        res<< "|<Ax,y>-<x,A'y>|/|Ax||y| = " << fabs(iprng - ipdom)/denom << endl;
        
        iwopeval.getDeriv2().applyOp(dm0,dm1,d1);
        iwopeval.getDeriv2().applyAdjOp(dm0,d1,bm1);
        iprng = d1.normsq();
        ipdom = dm1.inner(bm1);
        denom = iprng;
        res<< "\n Deriv=2: \n";
        res<< "<Ax,y>  = " << iprng << endl;
        res<< "<x,A'y> = " << ipdom << endl;
        res<< "|<Ax,y>-<x,A'y>|/|Ax||y| = " << fabs(iprng - ipdom)/denom << endl;
*/        
      
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
