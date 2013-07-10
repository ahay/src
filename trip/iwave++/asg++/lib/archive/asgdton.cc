#include "asgdton.hh"

namespace ASG {
  using TSOpt::IWaveEnvironment;
  using TSOpt::IWaveState;
  using TSOpt::IWaveLinState;
  using TSOpt::IWaveStep;
  using TSOpt::IWaveStaticInit;
  using TSOpt::IWaveDynamicInit;
  using TSOpt::IWaveOp;
  using TSOpt::Sampler;
  using TSOpt::LinSampler;
  using TSOpt::Sim;
  using TSOpt::StdSim;
  using TSOpt::StdSimData;
  using TSOpt::OpNewCreatePolicy;
  using TSOpt::PolicyBase;
  using TSOpt::CPSim;

#ifdef IWAVE_USE_MPI
  using TSOpt::MPIGridSpace;
  using TSOpt::MPISEGYSpace;
#else
  using TSOpt::GridSpace;
  using TSOpt::SEGYSpace;
#endif
  using TSOpt::segytrace;

  using RVL::LocalDataContainer;
  using RVL::StdProductSpace;
  using RVL::Vector;
  using RVL::Components;
  using RVL::Operator;
  using RVL::LinearOp;
  using RVL::OperatorEvaluation;
  using RVL::RVLException;
  using RVL::AssignFilename;
  using RVL::AssignTag;
  using RVL::AssignParams;
  using RVL::BinaryLocalFunctionObject;

  void FwdInt::operator()(LocalDataContainer<float> & out,
			  LocalDataContainer<float> const & in) {
    try {
      segytrace & trout = dynamic_cast<segytrace &>(out);
      segytrace const & trin = dynamic_cast<segytrace const &>(in);

      if (trout.getMetadata().ns != trin.getMetadata().ns) {
	RVLException e;
	e<<"Error: FwdInt::operator()\n";
	e<<"input, output traces of different lengths";
	throw e;
      }

      if (trout.getMetadata().delrt != trin.getMetadata().delrt) {
	RVLException e;
	e<<"Error: FwdInt::operator()\n";
	e<<"input, output traces have different time origins";
	throw e;
      }

      if (trout.getMetadata().dt != trin.getMetadata().dt) {
	RVLException e;
	e<<"Error: FwdInt::operator()\n";
	e<<"input, output traces have different time steps";
	throw e;
      }

      int ns=trout.getMetadata().ns;
      float hdt=0.0005f*((float)trout.getMetadata().dt);

      trout.getData()[0]=0.0f;

      for (int i=1;i<ns;i++) {
	trout.getData()[i]=trout.getData()[i-1]+
	  hdt*(trin.getData()[i-1]+trin.getData()[i]);
      }
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: AdjInt::operator()\n";
      e<<"either input or output LDC not segytrace\n";
      throw e;
    }
  }

  void AdjInt::operator()(LocalDataContainer<float> & out,
			  LocalDataContainer<float> const & in) {

    try {
      segytrace & trout = dynamic_cast<segytrace &>(out);
      segytrace const & trin = dynamic_cast<segytrace const &>(in);
      
      if (trout.getMetadata().ns != trin.getMetadata().ns) {
	RVLException e;
	e<<"Error: AdjInt::operator()\n";
	e<<"input, output traces of different lengths";
	throw e;
      }
      
      if (trout.getMetadata().delrt != trin.getMetadata().delrt) {
	RVLException e;
	e<<"Error: AdjInt::operator()\n";
	e<<"input, output traces have different time origins";
	throw e;
      }
      
      if (trout.getMetadata().dt != trin.getMetadata().dt) {
	RVLException e;
	e<<"Error: AdjInt::operator()\n";
	e<<"input, output traces have different time steps";
	throw e;
      }
      
      int ns=trout.getMetadata().ns;
      float hdt=0.0005f*((float)trout.getMetadata().dt);
      
      trout.getData()[ns-1]=hdt*trin.getData()[ns-1];
      
      for (int i=ns-2;i>0;i--) {
	trout.getData()[i]=trout.getData()[i+1]+
	  hdt*(trin.getData()[i+1]+trin.getData()[i]);
      }
      
      trout.getData()[0]=0.5*(hdt*trin.getData()[1]+trout.getData()[1]);
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: AdjInt::operator()\n";
      e<<"either input or output LDC not segytrace\n";
      throw e;
    }
  }


  void TimeRev::operator()(LocalDataContainer<float> & out,
			   LocalDataContainer<float> const & in) {
    try {

      segytrace & trout = dynamic_cast<segytrace &>(out);
      segytrace const & trin = dynamic_cast<segytrace const &>(in);

      if (trout.getMetadata().ns != trin.getMetadata().ns) {
	RVLException e;
	e<<"Error: TimeRev::operator()\n";
	e<<"input, output traces of different lengths";
	throw e;
      }

      if (trout.getMetadata().delrt != trin.getMetadata().delrt) {
	RVLException e;
	e<<"Error: TimeRev::operator()\n";
	e<<"input, output traces have different time origins";
	throw e;
      }

      if (trout.getMetadata().dt != trin.getMetadata().dt) {
	RVLException e;
	e<<"Error: TimeRev::operator()\n";
	e<<"input, output traces have different time steps";
	throw e;
      }

      int ns=trout.getMetadata().ns;

      for (int i=0;i<ns;i++) 
	trout.getData()[i]=trin.getData()[ns-i-1];
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: TimeRev::operator()\n";
      e<<"either input or output LDC not segytrace\n";
      throw e;
    }
  }

  ASGDtoN::ASGDtoN(Vector<float> const & _mdl,
#ifdef IWAVE_USE_MPI
		   MPISEGYSpace const & _dom,
#else
		   SEGYSpace const & _dom,
#endif	  
		   PARARRAY const & par,
		   FILE * _str)
    : mdl(_mdl), dom(_dom), w(dom),
      str(_str), integrate(false), 
      op(mdl.getSpace(),dom,par,str,&asg_gfdm) {

    /* sanity test - is srctype specified and = array*/
    char * tst;
    if (ps_ffcstring(par,"srctype",&tst)) {
      RVLException e;
      e<<"Error: ASGDtoN constructor\n";
      e<<"srctype key not set in param array\n";
      e<<"check parameter file\n";
      throw e;
    }
    if (strcmp(tst,"array")) {
      RVLException e;
      e<<"Error: ASGDtoN constructor\n";
      e<<"srctype != array\n";
      e<<"check parameter file\n";
      throw e;
    }
    free(tst);

    /* check for integration */
    int flute=0;
    ps_ffint(par,"integrate",&flute);
    if (flute) integrate=true;

    /* add filename for source workspace to param table */
    AssignParams ap(op.getPar(),"source");
    w.eval(ap);
    
  }

  ASGDtoN::ASGDtoN(ASGDtoN const & a) 
    : mdl(a.mdl), dom(a.dom), w(dom), 
      str(a.str), integrate(a.integrate),
      op(a.op) {
    /* add filename for source workspace to param table */
    AssignParams ap(op.getPar(),"source");
    w.eval(ap);
  }

  void ASGDtoN::apply(Vector<float> const & x,
		      Vector<float> & y) const {
    try { 

      /* integrate if flag set */
      if (integrate) {
	FwdInt fi;
	w.eval(fi,x);
      }
      else {
	w.copy(x);
      }

      OperatorEvaluation<float> opeval(op,mdl);

      y.copy(opeval.getValue());
      
    }
    catch (RVLException & e) {
      e<<"\ncalled from ASGDtoN::applyOp\n";
      throw e;
    }
  }
  
  void ASGDtoN::applyAdj(Vector<float> const & x,
			   Vector<float> & y) const {
    try { 

      /* time rev FO */
      TimeRev tr;

      /* time-reverse input */
      w.eval(tr,x);

      /* evaluate, copy results to output */
      OperatorEvaluation<float> opeval(op,mdl);
      w.eval(tr,opeval.getValue());
      
      if (integrate) {
	AdjInt ai;
	y.eval(ai,w);
      }
      else {
	y.copy(w);
      }
      
    }
    catch (RVLException & e) {
      e<<"\ncalled from ASGDtoN::applyAdjOp\n";
      throw e;
    }
  }

  ostream & ASGDtoN::write(ostream & str) const {
    str<<"ASGDtoN: Acoustic Staggered Grid Dirichlet-to-Neumann Op\n";
    return str;
  }
  
}

