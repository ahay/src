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

int main() {

  try {
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *            BEGIN RVL UNIT TEST 17             *"<<endl;
    cout<<" * Purpose: compare one-, two-line versions of   *"<<endl;
    cout<<" *             y = Ax + by                       *"<<endl;
    cout<<" * (complex double version)                      *"<<endl;
    cout<<" *************************************************/"<<endl; 
    cout<<endl;
    cout<<"1. Construct domain, range spaces"<<endl;
    int rows=121;
    int cols=121;
    RnSpace<std::complex<double> > dom(cols);
    RnSpace<std::complex<double> > rng(rows);
    cout<<endl;
    cout<<"2. Build GenMat (square matrix), initialize with rands\n";
    GenMat<std::complex<double> > A(dom,rng);
    RVLRandomize<std::complex<double> > rnd;
    A.eval(rnd);

    cout<<endl;
    cout<<"3. Construct input vector in domain, output in range\n";
    Vector<std::complex<double> > v(dom);
    Vector<std::complex<double> > w1(rng);
    Vector<std::complex<double> > w2(rng);
    Vector<std::complex<double> > tmp(rng);
    v.eval(rnd);
    w1.eval(rnd);
    w2.copy(w1);

    cout<<"4. y = Ax+by in two lines\n";
    std::complex<double>  a=2.0f;
    std::complex<double>  b=3.0f;

    A.applyOp(v,tmp);
    w1.linComb(a,tmp,b);

    cout<<"5. y = Ax+by in one line\n";
    A.applyOp(a,v,b,w2);

    Vector<std::complex<double> > dw(rng);
    dw.copy(w1);
    dw.linComb(-1.0,w2);
    double  t = dw.norm();

    ofstream str("ut17.out");
    str<<"norm of difference = "<<t<<endl;
    str.close(); 

    bool res = (abs(t) < 100.0 * numeric_limits<double>::epsilon());
 
    cout<<endl;
    cout<<"test result = "<<res<<endl;
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 17             *"<<endl;
    cout<<" *************************************************/"<<endl; 

    if (res) exit(0);
    exit(1);

  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}

