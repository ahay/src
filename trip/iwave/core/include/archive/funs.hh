#ifndef __DMODEL_FUNS__
#define __DMODEL_FUNS__

#include "mtuple.hh"

int get_num_fields() {
  int num=0;
  while ((get_iwave_fields()[num].field != "") && num<IWAVEMAXDATA) num++;
  if (num >= IWAVEMAXDATA) {
    RVLException e;
    e<<"Error: get_num_fields\n";
    e<<"  over limit for number - probably left off last entry with field=\"\"\n";
    throw e;
  }
  return num;
}

int get_num_iokeys() {
  int num=0;
  while ((get_iwave_iokeys()[num].keyword != "") && num<IWAVEMAXDATA) num++;
  if (num >= IWAVEMAXDATA) {
    RVLException e;
    e<<"Error: get_num_iokeys\n";
    e<<"  over limit for number - probably left off last entry with keyword=\"\"\n";
    throw e;
  }
  return num;
}

ostream & write_iwave_fields(ostream & str) {
  str <<"Field Definition: name = "<<get_iwave_model()<<"\n";
  for (int i=0;i<get_num_fields();i++) {
    str<<"  field["<<i<<"]="<<get_iwave_fields()[i].field
       <<" dynamic="<<get_iwave_fields()[i].dynamic
       <<" substep="<<get_iwave_fields()[i].substep
       <<" gtype=[";
    for (int j=0; j<RARR_MAX_NDIM-1; j++) 
      str<<get_iwave_fields()[i].gtype[j]<<",";
    str<<get_iwave_fields()[i].gtype[RARR_MAX_NDIM-1]<<"\n";
  }
  return str;
}

ostream & write_iwave_iokeys(ostream & str) {
  str <<"IO Definition: name = "<<get_iwave_model()<<"\n";
  for (int i=0;i<get_num_iokeys();i++) {
    str<<"  keyword["<<i<<"]="<<get_iwave_iokeys()[i].keyword
       <<" index="<<get_iwave_iokeys()[i].rarrindex
       <<" input="<<get_iwave_iokeys()[i].input
       <<" active="<<get_iwave_iokeys()[i].active
       <<"\n";
  }
  return str;
}

#endif
