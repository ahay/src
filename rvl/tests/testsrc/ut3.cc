/*************************************************************************

Copyright Rice University, 2004.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/

#include "std_cpp_includes.hh"
#include "rn.hh"
#include "productdata.hh"
#include "functions.hh"

using namespace RVL;

/** RVL Unit Test 3
    Purpose: direct test of StdProductDC behaviour (as opposed 
    to indirect tests via ProductSpace in UT4). Also tests
    value extraction (method getValue) and retention in 
    UnaryFOScalarRedn class.
*/
    
int main() {

  try {
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *            BEGIN RVL UNIT TEST 3              *"<<endl;
    cout<<" * Purpose: direct test of StdProductDC          *"<<endl;
    cout<<" * behaviour (as opposed to indirect tests       *"<<endl;
    cout<<" * via ProductSpace in UT4). Also tests          *"<<endl;
    cout<<" * value extraction (method getValue) and        *"<<endl;
    cout<<" * retention in UnaryFOScalarRedn class.         *"<<endl;
    cout<<" *************************************************/"<<endl; 
    cout<<endl;
    cout<<"1. Construct RnArray<float>s of lengths 5 and 10 \n";  
    RnArray<float> dc1(5);
    RnArray<float> dc2(10);
    cout<<endl;
    cout<<"2. Construct RVLAssignConst, val = 2\n";
    float c=2;
    RVLAssignConst<float> f(c);
    cout<<endl;
    cout<<"3. Construct StdProductDataContainer p, push RnArrays onto \n";
    cout<<"   so that it has size 2\n";
    // revision 1209: conform to new StdPDC design
    StdProductDataContainer p;
    RnDataContainerFactory<float> fac5(5);
    RnDataContainerFactory<float> fac10(10);
    p.push(fac5);
    p.push(fac10);
    cout<<endl;
    cout<<"4. Assign val = 2 to all components of the StdProductDC\n";
    cout<<"   by evaluation of the FO constructed in step 2\n";
    vector<DataContainer const *> jnk(0);
    p.eval(f,jnk);
    cout<<endl;
    cout<<"5. Construct another RVLAssignConst<float>, val = 3\n";
    c=3;
    RVLAssignConst<float> g(c);
    cout<<endl;
    cout<<"6. Build two more RnArrays of lengths 5, 10 resp., \n";
    cout<<"   another StdProductDC q, and push the RnArrays onto\n";
    cout<<"   q, so that it also has length 2 and components compatible\n";
    cout<<"   with those of p\n";
    StdProductDataContainer q;
    q.push(fac5);
    q.push(fac10);
    cout<<endl;
    cout<<"7. Evaluate the FO constructed in step 5 on q, so that all\n";
    cout<<"   components have value 3\n";
    q.eval(g,jnk);
    cout<<endl;
    cout<<"8. Construct RVLL2innerProd FO\n";
    RVLL2innerProd< float> h;
    cout<<endl;
    cout<<"9. Evaluate RVLL2innerProd on p, q\n";
    vector<DataContainer const *> tmp(1);
    tmp[0]=&q;
    p.eval(h,tmp);
    cout<<endl;
    float ip = h.getValue();
    cout<<"10. Extract value of inner product, using getValue method of\n";
    cout<<"    UnaryFOScalarRedn class - value = "<<ip
	<<", should = 90\n"; 
    cout<<endl;
    if (fabs(ip-90.0) > 100.0*numeric_limits<float>::epsilon()) {
      RVLException e;
      e<<"Error: ut3\n";
      e<<"inner product should = 90.0, but = "<<ip<<"\n";
      throw e;
    }
    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 3              *"<<endl;
    cout<<" *************************************************/"<<endl; 
  }
  catch (RVLException & e) {
    e.write(cout);
    exit(1);

  }
}

