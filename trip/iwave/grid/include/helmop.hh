#ifndef __TSOPT_HELMOP__
#define __TSOPT_HELMOP__

#include "linop_base.hh"
#include "space.hh"
#include "rnmat.hh"

#ifdef IWAVE_USE_MPI
#include "mpigridpp.hh"
#else
#include "gridpp.hh"
#endif

#include <gridio.h>
#include <f2c.h>
#include <stdio.h>
#include <iostream>
using namespace std;

/* helmholtz power function */
extern "C" void helm_(int,         /* bc */
                      integer *,   /* n1 */
                      integer *,   /* n2 */
                      float *,     /* d1 */
                      float *,     /* d2 */
                      float *,     /* w1 */
                      float *,     /* w2 */
                      float *,     /* p */
                      float *,     /* datum */
                      float *,     /* data in */
                      float *,     /* data out */
                      float *,     /* workspace */
                      integer *,   /* length of workspace */
                      integer *    /* error flag */
);

/* lenwork must be > 6*n1*n2+3*max(n2,2*n1)+21 */


namespace TSOpt {
    

    using namespace RVL;
    using namespace RVLAlg;
    using namespace RVLUmin;
    using RVL::ContentPackage;
    using RVL::LocalDataContainer;
    
#ifdef IWAVE_USE_MPI
    typedef MPIGridSpace myGridSpace;
#else
    typedef GridSpace myGridSpace;
#endif
    
    class HelmFO: public BinaryLocalFunctionObject<ireal> {
        
    private:
        ireal scale1, scale2;
        ireal power, datum;
        int DirichletSides;
        IPNT n_arr;
        RPNT d_arr;
        HelmFO();
        
    public:
        HelmFO(IPNT const & _narr,
               RPNT const & _darr,
               ireal _scale1=1.0f,
               ireal _scale2=1.0f,
               ireal _power=0.0f,
               ireal _datum=0.0f,
               int _DirichletSides=0)
        : scale1(_scale1),scale2(_scale2),power(_power), datum(_datum), DirichletSides(_DirichletSides){
        IASN(n_arr,_narr);
        RASN(d_arr,_darr);
        }
        
        HelmFO(HelmFO const & f)
        : scale1(f.scale1), scale2(f.scale2), power(f.power), datum(f.datum), DirichletSides(f.DirichletSides){
        IASN(n_arr,f.n_arr);
        RASN(d_arr,f.d_arr);
        }
        
        using RVL::LocalEvaluation<ireal>::operator();
        void operator()(LocalDataContainer<ireal> & x,
                        LocalDataContainer<ireal> const & y);
        
        string getName() const { string tmp = "HelmFO"; return tmp; }
        
    };
    void HelmFO::operator()(LocalDataContainer<ireal> & x,
                            LocalDataContainer<ireal> const & y){
        try{
        float *indata=NULL;
        float *outdata=NULL;
        float *work=NULL;
        integer f2c_n1;
        integer f2c_n2;
        integer lenwork;
        ContentPackage<ireal, RARR>  & gx =
        dynamic_cast<ContentPackage <ireal, RARR>  &> (x);
        ContentPackage<ireal, RARR> const & gy =
        dynamic_cast<ContentPackage <ireal, RARR> const &> (y);
        
        // precondition - metadata are same dimn
        RARR  & rax = gx.getMetadata();
        RARR const & ray = gy.getMetadata();
        int dimx; int dimy;
        int lendom;
        ra_ndim(&rax,&dimx);
        ra_ndim(&ray,&dimy);
        //cerr << "\n xdim=" << dimx << endl;
        //cerr << "\n ydim=" << dimy << endl;
        if (dimx != dimy) {
            RVLException e;
            e<<"Error: HelmFO::operator()\n";
            e<<"arguments have different dims:\n";
            e<<"dimx="<<dimx<<" dimy="<<dimy<<"\n";
            throw e;
        }
        
        // compute grid params
        IPNT gsx; IPNT gex;
        IPNT gsy; IPNT gey;
        IPNT s; IPNT e;
        ra_a_gse(&rax,gsx,gex);
        ra_a_gse(&ray,gsy,gey);
//        cerr << "\n===========================\n";
//        cerr << "\n gsx[0]=" << gsx[0] << endl;
//        cerr << "\n gex[0]=" << gex[0] << endl;
//        cerr << "\n gsx[1]=" << gsx[1] << endl;
//        cerr << "\n gex[1]=" << gex[1] << endl;
//        cerr << "\n===========================\n";
//        cerr << "\n gsy[0]=" << gsy[0] << endl;
//        cerr << "\n gey[0]=" << gey[0] << endl;
//        cerr << "\n gsy[1]=" << gsy[1] << endl;
//        cerr << "\n gey[1]=" << gey[1] << endl;
//        cerr << "\n===========================\n";
        // calculate grid overlap
        for (int ii=0;ii<dimx;ii++)  {
            s[ii]=max(gsy[ii],gsx[ii]);
            e[ii]=min(gey[ii],gex[ii]);
        }
        
        f2c_n1 = n_arr[0];
        f2c_n2 = n_arr[1];
        lendom=f2c_n1*f2c_n2;
        float _scale1=scale1;
        float _scale2=scale2;
        float _power=power;
        float _datum=datum;
        integer iter=0;
        
        // initialize workspace
        lenwork = 6*n_arr[1]*n_arr[0]+3*iwave_max(n_arr[1],2*n_arr[0])+21;
//        cerr << "\n lenwork=" << lenwork << endl;
//        cerr << "\n length of data = " << get_datasize_grid(gdom) << endl;
//        cerr << "\n n_arr[0] = " << n_arr[0] << endl;
//        cerr << "\n n_arr[1] = " << n_arr[1] << endl;
        //cerr << "\n physical domain size=" << lendom << endl;
        //cerr << "\n retrieveGlobalRank()=" << retrieveGlobalRank() << endl;        
        if (!(work = (float *)malloc(lenwork*sizeof(float)))) {
            RVLException e;
            e<<"Error: HelmOp::apply - failed to allocate " << lenwork << " floats for work buffer\n";
            throw e;
        }
        // allocate data arrays
        if (!(indata = (float *)malloc(lendom*sizeof(float)))) {
            RVLException e;
            e<<"Error: HelmOp::apply - failed to allocate " << lendom << " floats for input data\n";
            throw e;
        }
        if (!(outdata = (float *)malloc(lendom*sizeof(float)))) {
            RVLException e;
            e<<"Error: HelmOp::apply - failed to allocate " << lendom << " floats for output data\n";
            throw e;
        }
        IPNT i;
        integer idx;
#if RARR_MAX_NDIM > 0
        if (dimx==1) {
            for (i[0]=s[0];i[0]<=e[0];i[0]++) {
                indata[i[0]-s[0]]=ray._s1[i[0]];
            }
        helm_(DirichletSides,&f2c_n1,&f2c_n2,
              &(d_arr[0]),&(d_arr[1]),
              &(_scale1),&(_scale2),
              &_power,&_datum,
              indata,
              outdata,
              work,
              &lenwork,
              &iter);
        fprintf(stderr, "\n indata [100] = %f\n", indata[100]);
        fprintf(stderr, "\n outdata [100] = %f\n", outdata[100]);
            for (i[0]=s[0];i[0]<=e[0];i[0]++) {
                rax._s1[i[0]]=outdata[i[0]-s[0]];
            }
        }
#endif
#if RARR_MAX_NDIM > 1
        if (dimx==2) {
            for (i[1]=s[1];i[1]<=e[1];i[1]++) {
                for (i[0]=s[0];i[0]<=e[0];i[0]++) {
                    idx = (i[1]-s[1])*n_arr[0] + i[0]-s[0];
                    indata[idx]=ray._s2[i[1]][i[0]];
                }
            }
        helm_(DirichletSides,&f2c_n1,&f2c_n2,
              &(d_arr[0]),&(d_arr[1]),
              &(_scale1),&(_scale2),
              &_power,&_datum,
              indata,
              outdata,
              work,
              &lenwork,
              &iter);
        fprintf(stderr, "\n indata [100] = %f\n", indata[100]);
        fprintf(stderr, "\n outdata [100] = %f\n", outdata[100]);
        // copy data back
            for (i[1]=s[1];i[1]<=e[1];i[1]++) {
                for (i[0]=s[0];i[0]<=e[0];i[0]++) {
                    idx = (i[1]-s[1])*n_arr[0] + i[0]-s[0];
                    rax._s2[i[1]][i[0]]=outdata[idx];
                }
            }
        }
#endif
#if RARR_MAX_NDIM > 2
        if (dimx==3) {
            //cerr << "\n dim3=" << e[2] << endl;
            for (i[2]=s[2];i[2]<=e[2];i[2]++) {
            for (i[1]=s[1];i[1]<=e[1];i[1]++) {
                for (i[0]=s[0];i[0]<=e[0];i[0]++) {
                    idx = (i[1]-s[1])*n_arr[0] + i[0]-s[0];
                    indata[idx]=ray._s3[i[2]][i[1]][i[0]];
                }
            }
            helm_(DirichletSides,&f2c_n1,&f2c_n2,
              &(d_arr[0]),&(d_arr[1]),
              &(_scale1),&(_scale2),
              &_power,&_datum,
              indata,
              outdata,
              work,
              &lenwork,
              &iter);
            // copy data back
            for (i[1]=s[1];i[1]<=e[1];i[1]++) {
                for (i[0]=s[0];i[0]<=e[0];i[0]++) {
                    idx = (i[1]-s[1])*n_arr[0] + i[0]-s[0];
                    rax._s3[i[2]][i[1]][i[0]]=outdata[idx];
                }
            }    
            }
        fprintf(stderr, "\n indata [100] = %f\n", indata[10]);
        fprintf(stderr, "\n outdata [100] = %f\n", outdata[10]);
        }
#endif
        if (dimx<1 || dimx>3) {
            RVLException e;
            e<<"Error: HelmFO::operator()\n";
            e<<"dim = "<<dimx<<" outside of admissible set {1, 2}\n";
            throw e;
        }
        }
    catch (bad_cast) {
      RVLException e;
      e<<"\nError: HelmFO::operator()\n";
      e<<"at least one arg is not ContentPackage<ireal,RARR>\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from HelmFO::operator()\n";
      throw e;
    }
        
    }
    
    
    
    template<typename Scalar>
    class HelmOp: public LinearOp<Scalar> {
    private:
        
        //GridSpace const & dom;
        Space<Scalar> const & dom;
        //grid const & gdom;
        IPNT n_arr;
        RPNT d_arr; 
        RPNT w_arr;
        Scalar power, datum;
        int DirichletSides;
        
        FILE * str;
        
        // default construction disabled
        HelmOp();
        
    protected:
        
        void apply(const Vector<Scalar> & x,
                   Vector<Scalar> & y) const {
            try {
                HelmFO fo(n_arr,d_arr,w_arr[0],w_arr[1],power,datum,DirichletSides);
                MPISerialFunctionObject<Scalar> mpifo(fo);
                y.eval(mpifo,x);
            }
            catch (RVLException & e) {
                e<<"\ncalled in HelmOp::apply\n";
                throw e;
            }
            
        }
        
        
        void applyAdj(const Vector<Scalar> & x,
                      Vector<Scalar> & y) const {
            try {
                RVLException e;
                e<<"Error: HelmOp::applyAdj\n";
                e<<"  method is not implemented!\n";
                throw e;
                //Helm^T;
            }
            catch (RVLException & e) {
                e<<"\ncalled in HelmOp::applyAdj\n";
                throw e;
            }
        }
        
    public:
        
        HelmOp(HelmOp const & A):
        dom(A.dom), str(A.str),
        power(A.power), datum(A.datum), DirichletSides(A.DirichletSides){
            RASN(w_arr,A.w_arr);
        }
        
        HelmOp(Space<Scalar> const & _dom,
               RPNT _w_arr,
               Scalar _power=0.0f,
               Scalar _datum=0.0f,
               int _DirichletSides=0):
        dom(_dom),
        power(_power), datum(_datum), DirichletSides(_DirichletSides){
        try{
        RASN(w_arr,_w_arr);
        //cerr << "\n before dynamic_cast ProductSpace form input\n ";
        ProductSpace<Scalar> const * pdom = dynamic_cast<ProductSpace<Scalar> const *>(&dom);
        myGridSpace const & gsp = dynamic_cast<myGridSpace const &> ((*pdom)[0]);
        if (retrieveGlobalRank() == 0){
           grid const & g = gsp.getGrid();
           get_d(d_arr,g);
           get_n(n_arr,g);
        } 
            // check dimension
/*            if (gdom.dim != 2) {
                RVLException e;
                e<<"Error: HelmOp::HelmOp - input grid is not 2D\n";
                throw e;
           }
*/      }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: HelmOp constructor\n";
      e<<"  either domain or range is neither product nor a GridSpace,\n";
      e<<"  or some component is not a GridSpace\n";
      throw e;
    }
    catch (RVLException & e) {
      e<<"\ncalled from HelmOp constructor\n";
      throw e;
    }
        }
        
        ~HelmOp() {}
        
        // this class is considered terminal, with no overrides foreseen,
        // so clone method is not virtual
        LinearOp<Scalar> * clone() const { return new HelmOp(*this); }
        
        // access to domain, range
        const Space<float> & getDomain() const { return dom; }
        const Space<float> & getRange() const { return dom; }
        
        ostream & write(ostream & str) const {
            str<<"HelmOp\n";
            return str;
        }
        
    };
}
#endif
