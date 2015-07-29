#ifndef __TSOPT_GRIDPP_TOP__
#define __TSOPT_GRIDPP_TOP__

#ifdef IWAVE_USE_MPI
#include "parserpp.hh"
#include "axispp.hh"
#include "axisfun.hh"
#include "mpigridpp.hh"
#include "gridfun.hh"
#else
#include "parserpp.hh"
#include "axispp.hh"
#include "axisfun.hh"
#include "gridpp.hh"
#include "gridfun.hh"
#endif

#endif
