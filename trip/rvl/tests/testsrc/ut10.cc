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
#include "fcnfo.hh"
#include "vdtokb.hh"

using namespace RVL;

/** RVL Unit Test 10
    Purpose: direct test of FcnFO
*/

    
int main() {

  try {
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *            BEGIN RVL UNIT TEST 10             *"<<endl;
    cout<<" * Purpose: test conversion of C functions       *"<<endl;
    cout<<" * to FOs                                        *"<<endl;
    cout<<" *************************************************/"<<endl; 
    cout<<endl;
    cout<<"1. Construct RnArray<float>s of length 1 for v, d \n";  
    RnArray<float> v(1);
    RnArray<float> d(1);
    RnArray<float> kappa(1);
    RnArray<float> bouy(1);
    cout<<endl;
    cout<<"2. Assign v, d\n";
    v.getData()[0]=1.5;
    d.getData()[0]=1000.0;
    kappa.getData()[0]=0.0;
    bouy.getData()[0]=0.0;

    ScalarFO2<float,k> kfo;
    ScalarFO3<float,dkdv> dkdvfo;
    ScalarFO3<float,dkdd> dkddfo;
    ScalarFO2<float,b> bfo;
    ScalarFO3<float,dbdv> dbdvfo;
    ScalarFO3<float,dbdd> dbddfo;
    
    // set up BlockFO vectors
    std::vector<FunctionObject *> f;
    f.push_back(&kfo);
    f.push_back(&bfo);

    std::vector<FunctionObject *> dff0;
    dff0.push_back(&dkdvfo);
    dff0.push_back(&dkddfo);
    std::vector<FunctionObject *> dff1;
    dff1.push_back(&dbdvfo);
    dff1.push_back(&dbddfo);
    std::vector<std::vector<FunctionObject *> > dff(2);
    dff.push_back(dff0);
    dff.push_back(dff1);

    std::vector<FunctionObject *> dfa0;
    dfa0.push_back(&dkdvfo);
    dfa0.push_back(&dbdvfo);
    std::vector<FunctionObject *> dfa1;
    dfa1.push_back(&dkddfo);
    dfa1.push_back(&dbddfo);
    std::vector<std::vector<FunctionObject *> > dfa;
    dfa.push_back(dfa0);
    dfa.push_back(dfa1);
    
    RnArray<float> dv(1);
    dv.getData()[0]=1;
    RnArray<float> dd(1);
    dd.getData()[0]=-100.0;
    RnArray<float> dk(1);
    dk.getData()[0]=0.0;
    RnArray<float> db(1);
    db.getData()[0]=0.0;
    
    vector<DataContainer const *> uvec(2);
    uvec[0]=&v;
    uvec[1]=&d;

    vector<LocalDataContainer<float> *> cvec(2);
    cvec[0]=&kappa;
    cvec[1]=&bouy;

    cvec[0]->eval(*(f[0]),uvec);
    cvec[1]->eval(*(f[1]),uvec);

    cout<<"kappa="<<kappa.getData()[0]<<endl;
    cout<<"bouy ="<<bouy.getData()[0]<<endl;

    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 10             *"<<endl;
    cout<<" *************************************************/"<<endl; 
  }
  catch (RVLException & e) {
    e.write(cout);
    exit(1);

  }
}

