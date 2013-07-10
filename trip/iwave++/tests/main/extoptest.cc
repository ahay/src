#include "parserpp.hh"
#include "gridpp.hh"
#include "linop_base.hh"
#include "op.hh"
#include "extop.hh"
#include "rnmat.hh"
#include "adjtest.hh"

using TSOpt::GridSpace;
using TSOpt::ExtOp;
using namespace RVL;

//template<typename Scalar>
//void ExtOpTest(){
//    
//    
//}
void BuildData(){
    FILE * fp;
    float * x;
    float * y;
    float * z;
    int i;
    
    // open, write grid file
    fp = fopen("test0.rsf","w");
    fprintf(fp,"n1=10 d1=2.0 o1=-3.0\n");
    fprintf(fp,"n2=1 d2=1.0 o2=0.0\n");
    fprintf(fp,"data_format=native_float\n");
    fprintf(fp,"data_type=NONE\n");
    fprintf(fp,"in=./test0.rsf@\n");
    fflush(fp);
    fclose(fp);
    
    x = (float*)malloc(10*sizeof(float));
    for (i=0;i<10;i++) x[i]=1;
    
    fp = fopen("test0.rsf@","w");
    fwrite(x,sizeof(float),10,fp);
    fflush(fp);
    fclose(fp);
    
    // open, write grid file
    fp = fopen("test1.rsf","w");
    fprintf(fp,"n1=10 d1=2.0 o1=-3.0\n");
    fprintf(fp,"n2=1 d2=1.0 o2=0.0\n");
    fprintf(fp,"data_format=native_float\n");
    fprintf(fp,"data_type=NONE\n");
    fprintf(fp,"in=./test1.rsf@\n");
    fflush(fp);
    fclose(fp);
    
    for (i=0;i<10;i++) x[i]=2.0;
    
    fp = fopen("test1.rsf@","w");
    fwrite(x,sizeof(float),10,fp);
    fflush(fp);
    fclose(fp);
    
    y = (float*)malloc(240*sizeof(float));
    for (i=0;i<240;i++) y[i]=3.0;
    
    // open, write grid file
    fp = fopen("test2.rsf","w");
    fprintf(fp,"n1=10 d1=2.0 o1=-3.0\n");
    fprintf(fp,"n2=12 d2=1.0 o2=0.0\n");
    fprintf(fp,"n3=2 d3=1.0 o3=0.0\n");
    fprintf(fp,"data_format=native_float\n");
    fprintf(fp,"data_type=NONE\n");
    fprintf(fp,"in=./test2.rsf@\n");
    fprintf(fp,"\ndim=1 gdim=3\n");
    fflush(fp);
    fclose(fp);
    
    fp = fopen("test2.rsf@","w");
    fwrite(y,sizeof(float),240,fp);
    fflush(fp);
    fclose(fp);
    
    
    z = (float*)malloc(240*sizeof(float));
    for (i=0;i<240;i++) z[i]=4.0;
    
    // open, write grid file
    fp = fopen("test3.rsf","w");
    fprintf(fp,"n1=10 d1=2.0 o1=-3.0\n");
    fprintf(fp,"n2=12 d2=1.0 o2=0.0\n");
    fprintf(fp,"n3=2 d3=1.0 o3=0.0\n");
    fprintf(fp,"data_format=native_float\n");
    fprintf(fp,"data_type=NONE\n");
    fprintf(fp,"in=./test3.rsf@\n");
    fprintf(fp,"\ndim=1 gdim=3\n");
    fflush(fp);
    fclose(fp);
    
    fp = fopen("test3.rsf@","w");
    fwrite(z,sizeof(float),240,fp);
    fflush(fp);
    fclose(fp);
}

//template < typedef Scalar >
//bool TestExtOp(){
//
//
//    
//}
int main(int argc, char ** argv) {
    
    typedef float Scalar;
    string xname = "test0.rsf";
    string xname2 = "test1.rsf";
    string yname = "test3.rsf";
    string yname2 = "test2.rsf";
    string datatype = "NONE";
    
    try{
        
        BuildData();
        // create domain and range space
        GridSpace dom(xname,datatype);
        GridSpace dom2(xname2,datatype);
        GridSpace rng(yname,datatype);
        GridSpace rng2(yname2,datatype);
        
        
        // create operator
        cout<<"==============================================\n";
        cout<<"Build ExtOp, do adjoint test using AdjointTest\n";
        cout<<"==============================================\n";

        // create input, output vectors
        Vector<Scalar> x(dom);
        Vector<Scalar> x2(dom2);
        Vector<Scalar> y(rng);
        Vector<Scalar> y2(rng2);
        
        // couple vectors to files
        AssignFilename xaf(xname);
        x.eval(xaf);
        AssignFilename yaf(yname);
        y.eval(yaf);
        AssignFilename xaf2(xname2);
        x2.eval(xaf2);
        AssignFilename yaf2(yname2);
        y2.eval(yaf2);
        
        ExtOp <Scalar>op(dom,rng,datatype,stderr);
        cout << "Operator 1: \n";
        RVLRandomize<Scalar> rnd;
        if(!AdjointTest<Scalar>(op,rnd,cout)){
            cout << "*** adjoint test failed ***\n";
        }
        
        ExtOp <Scalar>op2(dom2,rng2,datatype,stderr);
        cout << "Operator 2: \n";
        RVLRandomize<Scalar> rnd2;
        if(!AdjointTest<Scalar>(op2,rnd2,cout)){
            cout << "*** adjoint test failed ***\n";
        }
        
        cout<<"==============================================\n";
        cout<<"Do adjoint test Manually\n";
        cout<<"==============================================\n";

        Vector<Scalar> y3(rng); y3.copy(y);
        Vector<Scalar> x3(dom); x3.copy(x);
        op.applyOp(x,y);
        cout << "Operator 1: \n";
        Scalar a = y.inner(y3);
        op.applyAdjOp(y3,x);
        Scalar b = x.inner(x3);
        if (fabs(a-b) < x3.norm()*y3.norm()*numeric_limits<Scalar>::epsilon())
            cout << " Passed adjoint test\n";
        else cout << " Adjoint test failed\n";
        cout << "    <Ax,y> = "<< a << endl;
        cout << "    <x,Ay> = "<< b << endl;

        
        Vector<Scalar> y4(rng2); y4.copy(y2);
        Vector<Scalar> x4(dom2); x4.copy(x2);
        
        cout <<"Operator 2: \n";
        op2.applyOp(x2,y2);
        a = y2.inner(y4);
        op2.applyAdjOp(y4,x2);
        b = x2.inner(x4);
        if (fabs(a-b) < x4.norm()*y4.norm()*numeric_limits<Scalar>::epsilon())
            cout << " Passed adjoint test\n";
        else cout << " Adjoint test failed\n";
        cout << "    <Ax,y> = "<< a << endl;
        cout << "    <x,Ay> = "<< b << endl;
        
    }
    catch(RVLException &e){
        e.write(cerr);
        exit(1);
    }
  // strings for names, data type
  }

