#include "asg_src.hh"

namespace ASG {

  ASGSource::ASGSource(IWaveState & _state, tracegeom & _tg)
      : tg(_tg), src(NULL), arr(NULL), state(_state), ready(false), ans(true)  {

    RASN(scoord,RPNT_0);
   
    // read source type from param table - default is "point"
    char * srctype = NULL;
    if (ps_flcstring(state.getPAR(),"srctype",&srctype)) {
      srctype=(char *)usermalloc_(6*sizeof(char));
      strcpy(srctype,"point");
    }
    
    if (!strcmp(srctype,"point")) {
      if (!(src=(POINTSRC *)usermalloc_(sizeof(POINTSRC)))) {
	RVLException e;
	e<<"ERROR: ASGSource constructor\n";
	e<<"failed to allocate point source. ABORT\n";
	throw e;
      }
      // this is OK, as a POINTSRC consists only of elementary types
      memset(src,0,sizeof(POINTSRC));
    }
    else if (!strcmp(srctype,"array")) {
      if (!(arr=(SAMPLER *)usermalloc_(sizeof(SAMPLER)))) {
	RVLException e;
	e<<"ERROR: ASGSource constructor\n";
	e<<"failed to allocate array source. ABORT\n";
	throw e;
      }

      // change of 03.12: for sources, all weights are 1!
      int err = 0;
      if ((err=sampler_construct(arr,
				&state.getPAR(),
				D_P,
				RPNT_1,
				RPNT_0,
				1,
				"source",
				"source",
				 state.getStream()))) {
	RVLException e;
	e<<"ERROR: ASGSource constructor from sampler_construct\n";
	e<<"returned code "<<err<<"\n";
	userfree_(arr);
	throw e;
      }
    }
    else {
      RVLException e;
      if (srctype) e<<"ERROR: ASGSource constructor: unknown source option = "<<srctype<<"\n";
      else e<<"ERROR: ASGSource constructor: unknown source option\n";
      throw e;
    }

    // this works because srctype is never referenced again without being
    // re-initialized (should have done this in asg driver!)
    if (srctype) userfree_(srctype);
  }

  ASGSource::~ASGSource() {
    if (src) {
      pointsrc_destroy(src);
      userfree_(src);
    }
    if (arr) {
      sampler_destroy(arr);
      userfree_(arr);
    }
  }

  void ASGSource::init() {

    int err=0;
    
    if (src) {

      err=pointsrc_init(src,
			&(state.IWaveState::getIWAVE().model),
			&(state.getPAR()),
			&tg,
			state.getStream());
      if (err) {
	RVLException e;
	e<<"ERROR: ASGSource constructor from pointsrc_init. ABORT\n";
	userfree_(src);
	throw e;
      }
      if (state.getIWAVE().printact > 1) {
	fprintf(state.getStream(),"------------ ASG_SRC: point source ----------\n");
	pointsrc_fprint(src,state.getStream());
      }
    }
    else {
      /*
      int towed=0;
  void ASGAdjSampler::getCoords(RPNT & s) const  { 
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
	e<<"ERROR: ASGSource from sampler_construct. ABORT\n";
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
	e<<"ERROR: ASGSource from sampler_init. ABORT\n";
	sampler_destroy(arr);
	userfree_(arr);
	throw e;
      }
      if (state.getIWAVE().printact > 1) {
	fprintf(state.getStream(),"------------ ASG_SRC: array source ----------\n");
	sampler_fprint(arr,state.getStream());
      }
    }

    ready=true;
  }

  int ASGSource::getStartIndex() {
    if (!ready) {
      RVLException e;
      e<<"Error: ASGSource::getStartIndex\n";
      e<<"Neither point source nor array source objects initialized\n";
      throw e;
    }
    if (src) return src->istart;
    if (arr) return (arr->t).istart;
    RVLException e;
    e<<"Error: ASGSource::getStartIndex\n";
    e<<"Neither point source nor array source objects initialized\n";
    throw e;
  }

  void ASGSource::run() {

    if (!ready) {
      RVLException e;
      e<<"Error: ASGSource::run\n";
      e<<"called with ready flag unset\n";
      if (src) {
	pointsrc_destroy(src);
	userfree_(src);
      }
      if (arr) {
	sampler_destroy(arr);
	userfree_(arr);
      }
      throw e;
    }

    // cerr<<"pre()\n";
    int iv = state.IWaveState::getTSIndex().getCstruct().iv;
#if 0
    int printact = state.IWaveState::getIWAVE().printact;
    FILE *stream = state.getStream();
    int it = state.IWaveState::getTSIndex().getCstruct().it;
    if (printact>1){
      fprintf(stream,"\n**********enter ASGSource::run for ");
      state.cid(stream);
      fprintf(stream,"\n");
    }
    if (printact>2){
      fprintf(stream,"\n----ASGSource::run() calling pointsorce_run() --it=%d iv=%d ------\n",it,iv);
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
    if (src && iv==0) {
      if (pointsrc_run(src,
		       &(state.IWaveState::getIWAVE().model))) ans=false;
    }
    if (arr && iv==0) {
      if (sampler_run(arr,
		      &(state.IWaveState::getIWAVE().model))) ans=false;
    }
#if 0
    if(printact>2){
      fprintf(stream,"\n----ASGSource::run() after pointsorce_run() --it=%d iv=%d ------\n",it,iv);
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
      fprintf(stream,"\n********** exit ASGSource::run for ");
      state.cid(stream);
      fprintf(stream,"\n");
    }
#endif 
  }

  void ASGSource::fprint() {
    if (src) pointsrc_fprint(src,state.getStream());
    if (arr) sampler_fprint(arr,state.getStream());
  }
  
}

