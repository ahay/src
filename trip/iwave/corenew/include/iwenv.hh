#ifndef __IW_ENV__
#define __IW_ENV__

#include "except.hh"
#include "parser.h"
#include "parserdecl.hh"
#include "iwave.h"

namespace TSOpt {
  using RVL::parse;
  using RVL::valparse;
  using RVL::RVLException;

  /** Function to initialize parallel environment, output stream, and param table */
  void IWaveEnvironment(int argc, char ** argv, int ts, PARARRAY ** par, FILE ** stream);

}

#endif
