#ifndef __TSOPT_EXTOP__
#define __TSOPT_EXTOP__

#include "parserpp.hh"
#include "linop_base.hh"

#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#else
#include "gridpp.hh"
#endif

namespace TSOpt {

using TSOpt::GridSpace;
using namespace RVL;

template<typename Scalar>
class ExtOp: public LinearOp<Scalar> {
private:
    
  GridSpace const & dom;
  GridSpace const & rng;
    grid const & gdom;
    grid const & grng;
    string datatype;
    
  FILE * str;

  // default construction disabled
  ExtOp();
    
protected:
    
    void apply(const Vector<Scalar> & x,
               Vector<Scalar> & y) const {
        try {
            // create empty param lists
            PARARRAY * xpar = ps_new();
            PARARRAY * ypar = ps_new();
            AssignParams xap(*xpar,str);
            AssignParams yap(*ypar,str);
            
            // transfer filenames from input and output vectors to param lists
            x.eval(xap);  y.eval(yap);
            
            // extract filenames from param lists
            string xname, yname;
            if (!parse<string>(*xpar,datatype,xname)) {
                RVLException e;
                e<<"Error: ExtOp::apply - failed to extract input filenames\n";
                throw e;
            }
            if (!parse<string>(*ypar,datatype,yname)) {
                RVLException e;
                e<<"Error: ExtOp::apply - failed to extract output filenames\n";
                throw e;
            }

            // Declare temperary working space
            string cmd, axisno,smplno;
            string tmpf1, tmpf2;
            stringstream con1,con2;
            tmpf1 = xname;  tmpf2 = "spraytmp.rsf";
            
            // Do spray if gdim > dim recusively
            if(grng.gdim > grng.dim){
                con1<<grng.dim+1;            axisno = con1.str();
                con2<<grng.axes[grng.dim].n; smplno = con2.str();

                cmd = "$RSFROOT/bin/sfspray < " + tmpf1 + " > " + tmpf2 + " axis=" + axisno + " n=" + smplno;
//                cout<<"ExtOp: executing command = "<<cmd<<endl;
                if(system(cmd.c_str())){
                    RVLException e;
                    e<<"Error: ExtOp::apply - failed to excute sfspray\n";
                    throw e;
                }
                
                tmpf1=tmpf2; tmpf2="spraytmp1.rsf";
                
                for(int i=grng.dim+2;i<=grng.gdim;++i){
                    con1.str(""); con1<<i;                axisno = con1.str();
                    con2.str(""); con2<<grng.axes[i-1].n; smplno = con2.str();
                
                    cmd = "$RSFROOT/bin/sfspray < " + tmpf1 + " > " + tmpf2 + " axis=" + axisno + " n="+ smplno;
                    tmpf1.swap(tmpf2);
                    if(system(cmd.c_str())){
                        RVLException e;
                        e<<"Error: ExtOp::apply - failed to excute sfspray\n";
                        throw e;
                    }
                }
                con1.str(""); con1<<grng.dim;  axisno = con1.str();
                con2.str(""); con2<<grng.gdim; smplno = con2.str();
                
                // write in dim and gdim info
                cmd = "$RSFROOT/bin/sfput < " + tmpf1 + " > " + tmpf2 +
                    " dim=" + axisno + "gdim=" + smplno;
                if(system(cmd.c_str())){
                    RVLException e;
                    e<<"Error: ExtOp::apply - failed to excute sfput\n";
                    throw e;
                }
            }
            // copy grid data to y
            Vector<Scalar> z(y.getSpace());
            AssignFilename af(tmpf2);
            z.eval(af);
            y.copy(z);
            
            // clean up
            system("/bin/rm spraytmp* $DATAPATH/spraytmp.rsf@");
            ps_delete(&xpar); ps_delete(&ypar);
        }
        catch (RVLException & e) {
            e<<"\ncalled in ExtOp::apply\n";
            throw e;
        }
        
    }

    
  void applyAdj(const Vector<Scalar> & x,
	     Vector<Scalar> & y) const {
    try {

      // create empty param lists
      PARARRAY * xpar = ps_new();
      AssignParams xap(*xpar,str);

      // transfer filenames from input and output vectors to param lists
      x.eval(xap);
        
      // extract filenames from param lists
      string xname;
      string yname;
      if (!parse<string>(*xpar,datatype,xname)) {
          RVLException e;
          e<<"Error: ExtOp::applyAdj - failed to extract filenames\n";
          throw e;
      }
        
      // build command
        
        string cmd;
        string axisno, tmpf1, tmpf2;
        stringstream convert;
        tmpf1 = xname;
        tmpf2 = "stacktmp.rsf";
        // Do stack if gdim > dim
        if ( grng.gdim > grng.dim){
            convert<<grng.dim+1;
            axisno = convert.str();
            // execute - should monitor for failure

            cmd = "$RSFROOT/bin/sfstack < " + tmpf1 + " > " + tmpf2 +
                " axis=" + axisno + " norm=n";
//            cout<<"ExtOp: executing command = "<<cmd<<endl;
            if(system(cmd.c_str())){
                RVLException e;
                e<<"Error: ExtOp::apply - failed to excute sfstack\n";
                throw e;
            }
            tmpf1=tmpf2;
            tmpf2="stacktmp1.rsf";
            for(int i=grng.dim+2;i<=grng.gdim;++i){
                cmd = "sfstack < " + tmpf1 + " > " + tmpf2 +
                    " axis=" + axisno + " norm=n"; 
                if(system(cmd.c_str())){
                    RVLException e;
                    e<<"Error: ExtOp::apply - failed to excute sfstack\n";
                    throw e;
                }
                tmpf1.swap(tmpf2);
            }
        }
      // copy grid data to y
      Vector<Scalar> z(y.getSpace());
      AssignFilename af(tmpf1);
      z.eval(af);
      y.copy(z);

      // clean up
      system("/bin/rm stacktmp* $DATAPATH/stacktmp.rsf@");
      ps_delete(&xpar);
    }
    catch (RVLException & e) {
      e<<"\ncalled in ExtOp::applyAdj\n";
      throw e;
    }
  }
    
public:

  ExtOp(ExtOp const & A): 
  dom(A.dom), rng(A.rng), str(A.str),
  gdom(A.gdom), grng(A.grng),
  datatype(A.datatype){
  }
      

  ExtOp(GridSpace const & _dom,
      GridSpace const & _rng,
        string dtype,
      FILE * _str):
    dom(_dom), rng(_rng), str(_str),
    gdom(dom.getGrid()),grng(rng.getGrid()),
                             datatype(dtype) {

            }

  ~ExtOp() {}

  // this class is considered terminal, with no overrides foreseen,
  // so clone method is not virtual
  LinearOp<Scalar> * clone() const { return new ExtOp(*this); }

  // access to domain, range
  const Space<float> & getDomain() const { return dom; }
  const Space<float> & getRange() const { return rng; }

  ostream & write(ostream & str) const {
    str<<"test operator for system call op design\n";
    return str;
  }

};
}
#endif