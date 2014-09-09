#ifndef __TSOpt_GridRWFO
#define __TSOpt_GridRWFO

#include "mpigridpp.hh"
#include "griddec.hh"
#include "local.hh"
#include "localreduction.hh"

namespace TSOpt {

  using RVL::writeMeta;
  using RVL::ContentPackage;
  using RVL::RVLException;
  using RVL::Writeable; 
  using RVL::UnaryLocalFunctionObject;
  using RVL::UnaryLocalFunctionObjectConstEval;
  using RVL::LocalDataContainer;
  using RVL::MPISynchRoot;

  /** GridReader
   * read content on a subgrid into memory   
  */

  template<typename DataType> 
#ifdef IWAVE_USE_MPI
  class GridReader: public UnaryLocalFunctionObjectConstEval<DataType>,
		    public MPISynchRoot
#else
  class GridReader: public UnaryLocalFunctionObjectConstEval<DataType>
#endif
  {
  
  private:
    Grid<DataType> subgrid;
    DataType * buf;
    int buf_size;
    mutable int local;

  public:
#ifdef IWAVE_USE_MPI    
    void set() {}
    virtual void synch() {}
#endif
    GridReader(Grid<DataType> const & _subgrid, DataType * _buf): subgrid(_subgrid), buf(_buf), buf_size(subgrid.getDataSize()), local(0){}
    GridReader(ContentPackage<DataType, Grid<DataType> > & md): subgrid(md.getMetadata()), buf(md.getData()), buf_size(subgrid.getDataSize()), local(0){}

    DataType * getData(){
      if (!buf) {
	RVLException e;
	e<<"Error: Reader::getData\n";
	e<<"data not initialized yet\n";
	throw e;
      }      
      return buf;
    }
    
    DataType const * getData() const {  
      if (!buf) {
	RVLException e;
	e<<"Error: Reader::getData\n";
	e<<"data not initialized yet\n";
	throw e;
      } 
      return buf;
    }

    // Grid<DataType> & getsubGrid(){
    //  return subgrid;
    //}
    
    Grid<DataType> const & getsubGrid() const {  
      return subgrid;
    }
    
    void reset(){
      local = 0;
    }
    
    // using RVL::LocalReduction<DataType>::operator();

    using RVL::UnaryLocalConstEval<DataType>::operator();
    void operator()(LocalDataContainer<DataType> const & x){
      try{
	ContentPackage<DataType, Grid<DataType> > const & gx = 
	  dynamic_cast<ContentPackage<DataType, Grid<DataType> > const &>(x);
	int gxsize = gx.getSize();
	int lcount = 0;
	/* make comparison and read information from subgrid; currently, assume 
	   subgrid belongs to grid; later, need to release this assumption via some 
	   sampling strategy
	 */
	if(buf){
	  int rem_size = buf_size - local;
	  if(rem_size) {
	    Grid<DataType> const & gxgrid = gx.getMetadata();
	    bool oinsub = iscoloinsubgrid(gxgrid);
	    if(oinsub){
	      DataType *tmp = buf + local;
	      DataType const * sval = gx.getData();
	      for(int i = 0; i< gxsize; i++) {
		if(isinsubaxis(gxgrid,i,1)){
		  tmp[i] = sval[i];
		  lcount++;
		}
		if (rem_size == lcount) break;
	      }
	      local += lcount;
	      lcount = 0;
	    }
	  }
	}
	else {
	  RVLException e;
	  e<<"Error: GridReader::operator()\n";
	  e<<"failed to read data into buffer: not enough memory allocated\n";
	}
	
      }
      catch (RVLException & e) {
	e<<" called from GridReader::operator()\n";
	e<<" called with incompatible data\n";
	throw e;
      }
    }  
    
    string getName() const { string tmp = "GridReader"; return tmp; }

    ~GridReader(){
      buf=NULL;
    }

    bool isinsubaxis(Grid<DataType> const &ggrid,int gindex, int gaxisid) const{
      int eid = ggrid.getExdAxis().id;
      if (gaxisid == eid){
	int sen = subgrid.getExdAxis().n;
	DataType seo = subgrid.getExdAxis().o;
	DataType sed = subgrid.getExdAxis().d;
	DataType geo = ggrid.getExdAxis().o;
	DataType ged = ggrid.getExdAxis().d;
	int subind = static_cast<int>(round(fabs(geo + gindex * ged - seo)/sed));
	if (subind < sen) return true;
	else return false;
      }
      else {
	int sn = subgrid.getAxis(gaxisid-1).n;
	DataType so = subgrid.getAxis(gaxisid-1).o;
	DataType sd = subgrid.getAxis(gaxisid-1).d;
	DataType go = ggrid.getAxis(gaxisid-1).o;
	DataType gd = ggrid.getAxis(gaxisid-1).d;
	int subind = static_cast<int>(round(fabs(go + gindex * gd - so)/sd));
	if (subind < sn) return true;
	else return false;
      }
    }
    
    /* 
       find out if the origin associated with 1D column in subgrid
    */
    bool iscoloinsubgrid(Grid<DataType> const &ggrid) const{
      if (!(isinsubaxis(ggrid,0,ggrid.getExdAxis().id))) return false;
      for(int i = RARR_MAX_NDIM; i > 1; i--)
	if(!(isinsubaxis(ggrid,0,i))) return false;
      return true;
    }

    
  };
  

  /** Writer
   * write content on a subgrid from memory to grid          
  */
  template<typename DataType> 
  class GridWriter: public UnaryLocalFunctionObject<DataType> {
    
  private:
    Grid<DataType> subgrid;
    DataType const * buf;
    int buf_size;
    int append_flag; 
    mutable int local;
    
  public:
    GridWriter(Grid<DataType> const & _subgrid, DataType const * _buf, int _apdflag = 0): subgrid(_subgrid), buf(_buf), buf_size(subgrid.getDataSize()), append_flag(_apdflag), local(0){ }
    GridWriter(ContentPackage<DataType, Grid<DataType> > const & md, int _apdflag = 0): subgrid(md.getMetadata()), buf(md.getData()), buf_size(subgrid.getDataSize()), append_flag(_apdflag), local(0){ }
    // GridWriter(GridReader<DataType> const & r): subgrid(r.getsubGrid()), buf(r.getData()), local(0){ }
    
    void set_apdflag(int _apdflag){
      append_flag = _apdflag;
    }
    
    void reset(){
      local = 0;
    }
    
    using RVL::LocalEvaluation<DataType>::operator();

    void operator()(LocalDataContainer<DataType> & x){
      try{
	ContentPackage<DataType, Grid<DataType> > & gx = 
	  dynamic_cast<ContentPackage<DataType, Grid<DataType> > &>(x);
	int gxsize = gx.getSize();
	int lcount = 0;
	/* make comparison and write information from subgrid to grid;
	   currently, assume subgrid belongs to grid; later, need to
	   release this assumption via some sampling strategy
	 */
	if(buf){
	  int rem_size = buf_size - local;
	  if (rem_size) {
	    Grid<DataType> const & gxgrid = gx.getMetadata();
	    bool oinsub = iscoloinsubgrid(gxgrid);
	    if(oinsub){
	      DataType const *tmp = buf + local;
	      DataType * pval = gx.getData();
	      if (!append_flag) {
		for(int i = 0; i< gxsize; i++) {
		  if(isinsubaxis(gxgrid,i,1)){
		    pval[i] = tmp[i];
		    lcount++;
		  }
		  if (rem_size == lcount) break;
		}
	      }
	      else {
		for(int i = 0; i< gxsize; i++) {
		  if(isinsubaxis(gxgrid,i,1)){
		    pval[i] += tmp[i];
		    lcount++;
		  }
		  if (rem_size == lcount) break;
		}
	      }
	      local += lcount;
	      lcount = 0;
	    }
	  }
	}
	else {
	  RVLException e;
	  e<<"Error: GridWriter::operator()\n";
	  e<<"failed to write data into buffer: data length not compatible\n";
	}
	
      }
      catch (RVLException & e) {
	e<<" called from Writer::operator()\n";
	e<<"called with incompatible data\n";
	throw e;
      }
    }    
    
    string getName() const { string tmp = "Writer"; return tmp; }
    
    ~GridWriter(){
      buf = NULL;
    }
    
    bool isinsubaxis(Grid<DataType> const &ggrid,int gindex, int gaxisid) const{
      int eid = ggrid.getExdAxis().id;
      if (gaxisid == eid){
	int sen = subgrid.getExdAxis().n;
	DataType seo = subgrid.getExdAxis().o;
	DataType sed = subgrid.getExdAxis().d;
	DataType geo = ggrid.getExdAxis().o;
	DataType ged = ggrid.getExdAxis().d;
	int subind = static_cast<int>(round(fabs(geo + gindex * ged - seo)/sed));
	if (subind < sen) return true;
	else return false;
      }
      else {
	int sn = subgrid.getAxis(gaxisid-1).n;
	DataType so = subgrid.getAxis(gaxisid-1).o;
	DataType sd = subgrid.getAxis(gaxisid-1).d;
	DataType go = ggrid.getAxis(gaxisid-1).o;
	DataType gd = ggrid.getAxis(gaxisid-1).d;
	int subind = static_cast<int>(round(fabs(go + gindex * gd - so)/sd));
	if (subind < sn) return true;
	else return false;
      }
    }
    
    /* 
       find out if the origin associated with 1D column in subgrid
    */
    bool iscoloinsubgrid(Grid<DataType> const &ggrid) const{
      if (!(isinsubaxis(ggrid,0,ggrid.getExdAxis().id))) return false;
      for(int i = RARR_MAX_NDIM; i > 1; i--)
	if(!(isinsubaxis(ggrid,0,i))) return false;
      return true;
    }
    

  };
  
  

}

#endif
