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

#include "rn.hh"
#include "functions.hh"

using namespace std;
using namespace RVL;

/** RVL Unit Test 1
    Purpose: exercise RnArray and all levels of the FO interface, 
    as well as the DataContainer hierarchy.
*/
    
int main() {

  try {

    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *            BEGIN RVL UNIT TEST 0              *"<<endl;
    cout<<" * Purpose: exercise DC and FO hierarchies using *"<<endl;
    cout<<" * RnArray concrete LDC class and FOs defined in *"<<endl;
    cout<<" * functions.H                                   *"<<endl;
    cout<<" *************************************************/"<<endl; 
    cout<<endl;
    cout<<"1. construct 10-diml RnArray<float>\n";
    RnArray<float> x(10);
    cout<<endl;
    cout<<"2. Assign constant val 1.0 to all components\n";
    RVLAssignConst<float> f(ScalarFieldTraits<float>::One());
    f(x);
    cout<<endl<<"result:"<<endl<<endl;
    x.write(cout);
    cout<<endl;
    cout<<"3. Assign const val 2.0 to all components, using UnaryFO\n";
    cout<<"   operator() interface\n";
    RVLAssignConst<float> g(2.0);
    g(x);
    cout<<endl<<"result:"<<endl<<endl;
    x.write(cout);
    cout<<endl;    
    cout<<"4. Compute inner product of vector with itself; result\n";
    cout<<"   should be = 40.0\n";
    RVLL2innerProd<float> ip;
    ip(x,x);
    cout<<endl<<"<x,x> = "<<ip.getValue()<<endl<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 0              *"<<endl;
    cout<<" *************************************************/"<<endl; 
  }
  catch (RVLException & e) {
    e.write(cout);
    exit(1);
  }
}






