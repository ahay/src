#include "iwenv.hh"

//#define IWAVE_VERBOSE

namespace TSOpt {

  void IWaveEnvironment(int argc, char ** argv, int ts,
			PARARRAY ** par, 
			FILE ** stream) {
    int err=0;
    //    cerr<<"IWaveEnvironment -> initparallel_global\n";
    initparallel_global(ts);
    //    cerr<<"IWaveEnvironment -> initinoutstream\n";
    if ((err=initoutstream(stream,retrieveGlobalRank(),retrieveGlobalSize()))) {
      RVLException e;
      e<<"Error: IWaveEnvironment from initoutstream, err="<<err<<"\n";
      e<<"failed to initialize output stream\n";
      throw e;
    }
    //    cerr<<"IWaveEnvironment -> readinput\n";
#ifdef IWAVE_VERBOSE
    for (int i=1;i<argc;i++) 
      cerr<<"IWaveEnvironment: argv["<<i<<"]="<<argv[i]<<endl;
#endif
    if ((err=readinput(par,*stream,argc,argv))) {
      RVLException e;
      e<<"Error: IWaveEnvironment from readinput, err="<<err<<"\n";
      e<<"failed to parse param file \n";
      throw e;
    }
#ifdef IWAVE_VERBOSE
    cerr<<"IWaveEnvironment: after construct par\n";
    ps_printall(*(*par),stderr);
#endif
    //    cerr<<"IWaveEnvironment -> initparallel_local\n";
    if ((err=initparallel_local(*(*par),*stream))) {
      RVLException e;
      e<<"Error: IWaveEnvironment from initparallel_local, err="<<err<<"\n";
      e<<"failed to create cartesian grid or local comm\n";
      throw e;
    }
    //    cerr<<"IWaveEnvironment -> exit\n";
  }

}

