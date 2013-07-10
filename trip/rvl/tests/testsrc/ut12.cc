/*************************************************************************

Copyright Rice University, 2011.
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

#include "rnmat.hh"
#include "adjtest.hh"

using namespace RVL;

/** RVL Unit Test 12
    Purpose: test detection of input/domain mismatch
*/
    
int main() {

  try {
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *            BEGIN RVL UNIT TEST 12             *"<<endl;
    cout<<" * Purpose: test sanity test in LinOp            *"<<endl;
    cout<<" *************************************************/"<<endl; 
    cout<<endl;
    cout<<"1. Construct domain, range spaces"<<endl;
    int rows=261;
    int cols=121;
    RnSpace<float> dom(cols);
    RnSpace<float> rng(rows);
    cout<<endl;
    cout<<"2. Build GenMat\n";
    GenMat<float> A(dom,rng);
    cout<<endl;
    cout<<"3. Construct input vector in wrong space\n";
    Vector<float> v(rng);
    Vector<float> w(rng);
    A.applyOp(v,w);

    cout<<" something wrong - did not detect\n";
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 12             *"<<endl;
    cout<<" *************************************************/"<<endl; 
  }
  catch (RVLException & e) {
    ofstream str("testsrc/ut12/ut12.out");
    cout<<endl;
    cout<<"successful conclusion - error trapped - see output in ut12/ut12.out\n";
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 12             *"<<endl;
    cout<<" *************************************************/"<<endl;     
    e.write(str);
    exit(0);

  }
}

