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
#include "rnspace.hh"
#include "functions.hh"
#include "productspace.hh"

using namespace RVL;

/** RVL Unit Test 1
    Purpose: exercise StdProductSpace, Components
*/

int main() {

  try {
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *            BEGIN RVL UNIT TEST 4              *"<<endl;
    cout<<" * Purpose: exercise StdProductSpace, Components *"<<endl;
    cout<<" *************************************************/"<<endl; 
    cout<<endl;
    cout<<"1. construct 5, 10-diml RnSpace<float>s\n";
    RnSpace<float> sp1(5);
    RnSpace<float> sp2(10);
    cout<<endl;
    cout<<"2. build StdProductSpace with these two spaces as components\n";
    StdProductSpace<float> psp(sp1,sp2);
    cout<<endl;
    cout<<"3. construct Vector in the StdProductSpace\n";
    Vector<float> pv(psp);
    cout<<endl;
    cout<<"4. initialize with val = 2.0 via RVLAssignConst FO eval\n";
    float c=2;
    RVLAssignConst<float> f(c);
    pv.eval(f);
    cout<<endl;
    cout<<"5. construct RVLAssignConst with val = 3.0\n";
    c=3;
    RVLAssignConst<float> g(c);
    {
      cout<<endl;
      cout<<"6. within a new scope, construct Components from Vector in \n";
      cout<<"   product space (step 3)\n";
      Components<float> cv(pv);
      cout<<endl;
      cout<<"7. write out components:\n";
      cv.write(cout);
      cout<<endl;
      cout<<"8. evaluate FO constructed in 5 on 2nd component, setting\n";
      cout<<"   its coordinates all to 3\n";
      cv[1].eval(g);
      cout<<endl;
      cout<<"9. exit scope, which destroys Components object\n";
    }
    cout<<endl;
    cout<<"10. write product vector - first component should still be 2's,\n";
    cout<<"    second component should now be 3's\n";
    pv.write(cout);
    cout<<endl;

    cout<<"11. construct another vector in same product space\n";
    Vector<float> qv(psp);
    cout<<endl;
    cout<<"12. initialize all coords to 3\n";
    qv.eval(g);
    cout<<endl;
    cout<<"13. write it out\n";
    qv.write(cout);
    cout<<endl;
    cout<<"14. compute L2 inner product = "<<pv.inner(qv)
	<<", value should = 120\n";
    cout<<"14a. compute L2 inner product of pv with itself - should be 110.";
    cout<<endl;
    cout<<"15. build another product space, all components having dim 5,\n";
    cout<<"    by default-constructing a StdProductSpace and pushing the \n";
    cout<<"    5D RnSpace constructed in step 1 onto it 3 times\n";
    StdProductSpace<float> rsp(sp1,sp1,sp1);
    cout<<endl;
    cout<<"16. construct a Vector in this 3-component space, construct its\n";
    cout<<"    Components object, initialize the first component to 1, second\n";
    cout<<"    to 2, third to 3, using appropriate FOs\n";
    Vector<float> rv(rsp);
    Components<float> crv(rv);
    c=1.0;
    RVLAssignConst<float> h(c);    
    crv[0].eval(h);
    crv[1].eval(f);
    crv[2].eval(g);
    cout<<endl;
    cout<<"18. write out the Vector\n";
    rv.write(cout);
    cout<<endl;
    cout<<"18.5 write out the Components\n";
    crv.write(cout);
    cout<<endl;
    cout<<"19. construct a vector in the first component space (5D RnSpace)\n";
    Vector<float> s(sp1);
    cout<<endl;
    cout<<"20. overwrite the sum of the three components of the Vector\n";
    cout<<"    constructed in steps 16 and 17 on the Vector construced in \n";
    cout<<"    step 19, using two applications of RVL::Vector::linComb\n";
    s.copy(crv[0]);
    s.linComb(ScalarFieldTraits<float>::One(),crv[1]);
    s.linComb(ScalarFieldTraits<float>::One(),crv[2]);
    cout<<endl;
    cout<<"21. write out the sum\n";
    s.write(cout);
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 4              *"<<endl;
    cout<<" *************************************************/"<<endl; 
  }
  catch (RVLException & e) {
    e.write(cout);
    exit(1);
  }
}

  
