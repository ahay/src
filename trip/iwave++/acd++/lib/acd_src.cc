#include "acd_src.hh"

namespace ACD {

  ACDSource::ACDSource(IWaveState & _state, tracegeom & _tg)
    : state(_state), arr(NULL), ready(false), ans(true), tg(_tg) {

    if (!(arr=(SAMPLER *)usermalloc_(sizeof(SAMPLER)))) {
      RVLException e;
      e<<"ERROR: ACDSource constructor\n";
      e<<"failed to allocate array source. ABORT\n";
      throw e;
    }
    
    // change of 03.12: for sources, all weights are 1!
    int err = 0;
    // repetition means only one sampling - this must be changed!!!!
    IPNT tmp;
    for (int i=0;i<RARR_MAX_NDIM;i++) tmp[i]=D_UC;
    if (err=sampler_construct(arr,
			      &state.getPAR(),
			      tmp,
			      RPNT_1,
			      RPNT_0,
			      -1,
			      "source",
			      "source",
			      state.getStream())) {
      RVLException e;
      e<<"ERROR: ACDSource constructor from sampler_construct\n";
      e<<"returned code "<<err<<"\n";
      userfree_(arr);
      throw e;
    }  
  }

  ACDSource::~ACDSource() {
    if (arr) {
      sampler_destroy(arr);
      userfree_(arr);
    }
  }

  void ACDSource::init() {

    int err=0;
    
    /*
      int towed=0;
      void ACDAdjSampler::getCoords(RPNT & s) const  { 
      for (int i=0;i<RARR_MAX_NDIM; i++) s[i]=(trace.t.tg.src[trace.t.tg.xrec])[i]; 
      }
  

      ps_flint(state.getPAR(),"towed",&towed);
      if (towed) RASN(scoord,tg.src[tg.xrec]);
      
      IPNT mindex;
      for (int i=0;i<RARR_MAX_NDIM;i++) mindex[i]=D_MP0;
      err=sampler_construct(arr,&state.getPAR(),D_P,mindex,scoord,1,
			    NULL,"source",state.getStream());
      if (err) {
	RVLException e;
	e<<"ERROR: ACDSource from sampler_construct. ABORT\n";
	userfree_(arr);
	throw e;
      }
    */
    
    err=sampler_init(arr,
		     &(state.IWaveState::getIWAVE().model),
		     &(state.getPAR()),
		     state.getStream());
    if (err) {
      RVLException e;
      e<<"ERROR: ACDSource from sampler_init. ABORT\n";
      sampler_destroy(arr);
      userfree_(arr);
      throw e;
    }
    if (state.getIWAVE().printact > 1) {
      fprintf(state.getStream(),"------------ ACD_SRC: array source ----------\n");
      sampler_fprint(arr,state.getStream());
    }
    ready=true;
  }

  int ACDSource::getStartIndex() {
    if (!ready) {
      RVLException e;
      e<<"Error: ACDSource::getStartIndex\n";
      e<<"Neither point source nor array source objects initialized\n";
      throw e;
    }
    if (arr) return (arr->t).istart;
    RVLException e;
    e<<"Error: ACDSource::getStartIndex\n";
    e<<"Array source object not initialized\n";
    throw e;
  }

  void ACDSource::run() {

    if (!ready) {
      RVLException e;
      e<<"Error: ACDSource::run\n";
      e<<"called with ready flag unset\n";
      if (arr) {
	sampler_destroy(arr);
	userfree_(arr);
      }
      throw e;
    }

    int iv = state.IWaveState::getTSIndex().getCstruct().iv;
#if 0
    int printact = state.IWaveState::getIWAVE().printact;
    FILE *stream = state.getStream();
    int it = state.IWaveState::getTSIndex().getCstruct().it;
    if (printact>1){
      fprintf(stream,"\n**********enter ACDSource::run for ");
      state.cid(stream);
      fprintf(stream,"\n");
    }
    if (printact>2){
      fprintf(stream,"\n----ACDSource::run() calling pointsorce_run() --it=%d iv=%d ------\n",it,iv);
      if (iv==0){
	fprintf(stream,"\n----print, only when iv==0\n ");
	//fprintf(stream,"\n----ref-domain--pressure array---------\n ");
	//rd_print(&(state.IWaveState::getIWAVE().model.ld_a),0,stream);
	rd_a_print(&(state.IWaveState::getIWAVE().model.ld_a),stream);
	//fprintf(stream,"\n----pert-domain--------------------\n ");
	//rd_a_print(&(state.getIWAVE().model.ld_a),stream);
      }
    }
#endif
    if (arr && iv==0) {
      if (sampler_run(arr,
		      &(state.IWaveState::getIWAVE().model))) ans=false;
    }
#if 0
    if(printact>2){
      fprintf(stream,"\n----ACDSource::run() after pointsorce_run() --it=%d iv=%d ------\n",it,iv);
      if(iv==0){
	fprintf(stream,"\n----print, only when iv==0 \n");
	//fprintf(stream,"\n----ref-domain--pressure array---------------\n ");
	rd_print(&(state.IWaveState::getIWAVE().model.ld_a),0,stream);
	//rd_a_print(&(state.IWaveState::getIWAVE().model.ld_a),stream);
	//fprintf(stream,"\n----pert-domain--------------------\n ");
	//rd_a_print(&(state.getIWAVE().model.ld_a),stream); 
      }
    }
    if(printact >1){
      fprintf(stream,"\n********** exit ACDSource::run for ");
      state.cid(stream);
      fprintf(stream,"\n");
    }
#endif 
  }

  void ACDSource::fprint() {
    if (arr) sampler_fprint(arr,state.getStream());
  }
  
}

