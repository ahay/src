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

/** RVL Unit Test 13
    Purpose: test LinearOpFO construction
*/
    
int main() {

  try {
    PlantSeeds(getpid());
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *            BEGIN RVL UNIT TEST 13             *"<<endl;
    cout<<" * Purpose: compare LinearOpFO with explicit     *"<<endl;
    cout<<" * construction                                  *"<<endl;
    cout<<" *************************************************/"<<endl; 
    cout<<endl;
    cout<<"1. Construct domain, range spaces"<<endl;
    int rows=261;
    int cols=121;
    RnSpace<float> dom(cols);
    RnSpace<float> rng(rows);
    cout<<endl;
    cout<<"2. Build matvec\n";
    matvec<float> m(rng.getSize(),dom.getSize());
    for (int i=0;i<rows;i++) {
      for (int j=0;j<cols;j++) {
	float x = float(rand())/float(RAND_MAX) -0.5;
	m.getElement(i,j)=x;
      }
    }
    cout<<endl;
    cout<<"3. Build GenMat\n";
    GenMat<float> A(dom,rng);
    for (int i=0;i<rows;i++) {
      for (int j=0;j<cols;j++) {
	A.setElement(i,j,m.getElement(i,j));
      }    
    }
    cout<<endl;
    cout<<"4. Build LinearOpFO\n";
    fmatvec<float> fwd(m);
    amatvec<float> adj(m);
    LinearOpFO<float> B(dom,rng,fwd,adj);
    cout<<endl;
    cout<<"5. Compare actions - norm of diff should = 0 exactly\n";
    Vector<float> inp1(dom);
    Vector<float> inp2(dom);
    Vector<float> outp1(rng);
    Vector<float> outp2(rng);

    RVLRandomize<float> rnd;
    inp1.eval(rnd);
    A.applyOp(inp1,outp1);
    B.applyOp(inp1,outp2);
    outp1.linComb(-1.0f,outp2);
    cout<<"norm of difference in fwd dir = "<<outp1.norm()<<endl;
    outp1.eval(rnd);
    A.applyAdjOp(outp1,inp1);
    B.applyAdjOp(outp1,inp2);
    inp1.linComb(-1.0f,inp2);
    cout<<"norm of difference in fwd dir = "<<inp1.norm()<<endl;

    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 13             *"<<endl;
    cout<<" *************************************************/"<<endl; 
  }
  catch (RVLException & e) {
    e.write(cout);
    exit(1);

  }
}

