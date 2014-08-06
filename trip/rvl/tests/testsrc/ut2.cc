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

#include "rnspace.hh"

using namespace RVL;

/** RVL Unit Test 2
    Purpose: exercise RnSpace comparison and inner product.
    constructors, also RVLAssignConst.
*/
int main() {

  try {
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *            BEGIN RVL UNIT TEST 2              *"<<endl;
    cout<<" * Purpose: exercise RnSpace comparison and      *"<<endl;
    cout<<" * Vector inner product methods, RVLAssignConst  *"<<endl;
    cout<<" *************************************************/"<<endl; 
    cout<<endl;
    cout<<"1. construct 5-diml RnSpace<float>\n";    
    int n=5;
    RnSpace<float> sp(n);
    cout<<endl;
    cout<<"2. construct RVL::Vector x in this space\n";
    Vector<float> x(sp);
    cout<<endl;
    cout<<"3. construct RVLAssignConst<float> FO, val = 1.0\n";
    float c = 1;
    RVLAssignConst<float> f(c);
    cout<<endl;
    cout<<"4. evaluate, which makes x the const vector, val = 1\n";
    x.eval(f);
    cout<<endl;
    cout<<"5. write out x\n";
    x.write(cout);
    cout<<endl;
    cout<<"6. compute norm squared of x = "<<x.inner(x)<<endl;
    cout<<endl;
    cout<<"7. construct a second vector y in space\n";
    Vector<float> y(x);
    cout<<endl;
    cout<<"8. construct RVLAssignConst FO, val = 2.0\n";
    c=2;
    RVLAssignConst<float> g(c);
    cout<<endl;
    cout<<"9. evaluate, which makes y the const.vector, val = 2\n";
    y.eval(g);
    cout<<endl;
    cout<<"10. write out y\n";
    y.write(cout);
    cout<<endl;
    cout<<"11. compute inner product of x, y = "<<x.inner(y)<<endl;
    cout<<"    should = 10\n";
    cout<<endl;
    cout<<"12. overwrite x with 1*x + 2*y, which should be const\n";
    cout<<"    vector, val = 5\n";
    float b = 2;
    x.linComb(b,y);
    cout<<endl;
    cout<<"13. write out altered x\n";
    x.write(cout);
    cout<<endl;
    cout<<"14. Create another Rn, dim = 4.\n";
    int m=4;
    RnSpace<float> osp(m);
    cout<<endl;
    cout<<"15. Create a vector z in the second Rn\n";
    Vector<float> z(osp);
    cout<<endl;
    cout<<"16. evaluate RVLAssignConst, val = 2.0, on z\n";
    z.eval(g);
    cout<<endl;
    cout<<"17. Attempt to evaluate inner product of x and z\n";
    cout<<"    This should FAIL (throw an exception) because the\n";
    cout<<"    spaces have different dimensions.\n";
    cout<<x.inner(z)<<endl;
    exit(1);
  }
  catch (RVLException & e) {
    e.write(cout);
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 2              *"<<endl;
    cout<<" *************************************************/"<<endl; 
    exit(0);
  }
}







