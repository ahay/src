#ifndef __ACD_DEFN__
#define __ACD_DEFN__

#include "acd.hh"
#include "iwinfo.hh"

std::string IWaveInfo::iwave_model = "acd";
FIELD IWaveInfo::iwave_fields[]
= {
  {"csq",    0,    0,  {0, 0, 0}},
  {"uc",     1,    0,  {0, 0, 0}},
  {"up",     1,    0,  {0, 0, 0}},
  {"",       0,    0,  {0, 0, 0}}
};

FD_MODELINIT IWaveInfo::minit = acd_modelinit;
FD_MODELDEST IWaveInfo::mdest = acd_modeldest;
FD_TIMESTEP IWaveInfo::timestep = acd_timestep;
FD_TIMEGRID IWaveInfo::timegrid = acd_timegrid;
FD_STENCIL IWaveInfo::createstencil = acd_create_sten;
FD_CHECK IWaveInfo::check = acd_check;

#endif

