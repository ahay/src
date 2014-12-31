#ifndef __ASG_DEFN__
#define __ASG_DEFN__

#include "asg.hh"
#include "iwinfo.hh"

std::string IWaveInfo::iwave_model = "asg";
FIELD IWaveInfo::iwave_fields[]
= {
  {"bulkmod",    0,    0,  {0, 0, 0}},
  {"buoyancy",   0,    0,  {0, 0, 0}},
  {"p0",         1,    0,  {0, 0, 0}},
  {"p1",         1,    0,  {0, 0, 0}},
  {"p2",         1,    0,  {0, 0, 0}},
  {"v0",         1,    1,  {1, 0, 0}},
  {"v1",         1,    1,  {0, 1, 0}},
  {"v2",         1,    1,  {0, 0, 1}},
  {"",           0,    0,  {0, 0, 0}}
};

FD_MODELINIT IWaveInfo::minit = asg_modelinit;
FD_MODELDEST IWaveInfo::mdest = asg_modeldest;
FD_TIMESTEP IWaveInfo::timestep = asg_timestep;
FD_TIMEGRID IWaveInfo::timegrid = asg_timegrid;
FD_STENCIL IWaveInfo::createstencil = asg_create_sten;
FD_CHECK IWaveInfo::check = asg_check;

#endif

