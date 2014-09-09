#ifndef __ACDPML_DEFN__
#define __ACDPML_DEFN__

#include "acdpml.hh"
#include "iwinfo.hh"

std::string IWaveInfo::iwave_model = "acdpml";
FIELD IWaveInfo::iwave_fields[]
= {
  {"csq",    0,    0,  {0, 0, 0}},
  {"uc",     1,    0,  {0, 0, 0}},
  {"up",     1,    0,  {0, 0, 0}},
  {"phi1",   1,    0,  {1, 1, 1}},
  {"phi0",   1,    0,  {1, 1, 1}},
  {"",       0,    0,  {0, 0, 0}}
};

FD_MODELINIT IWaveInfo::minit = acdpml_modelinit;
FD_MODELDEST IWaveInfo::mdest = acdpml_modeldest;
FD_TIMESTEP IWaveInfo::timestep = acdpml_timestep;
FD_TIMEGRID IWaveInfo::timegrid = acdpml_timegrid;
FD_STENCIL IWaveInfo::createstencil = acdpml_create_sten;
FD_CHECK IWaveInfo::check = acdpml_check;

#endif

