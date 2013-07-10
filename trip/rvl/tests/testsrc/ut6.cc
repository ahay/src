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

int main() {
  try {
    //    srand(getpid());
    //    srand(19490615);
    PlantSeeds(19490615);

    cout<<"/*************************************************"<<endl;
    cout<<" *            BEGIN RVL UNIT TEST 6              *"<<endl;
    cout<<" * Purpose: test ASCII i/o FOs for RnArray       *"<<endl;
    cout<<" *************************************************/"<<endl;

    int dim=10;
    cout<<endl<<"1. construct RnSpace<double>, Vector in it, random initialization\n";
    RnSpace<double> sp(dim);
    Vector<double> vec1(sp);
    RVLRandomize<double> f;
    vec1.eval(f);
    cout<<endl<<"2. construct ASCIIWriter FO for file ut6.dat\n"; 
    ASCIIWriter<double> svr("./ut6.dat");
    cout<<endl<<"3. save Vector to ut6/ut6.dat by evaluating FO\n";
    vec1.eval(svr);
    cout<<endl<<"4. construct a second vector in the space\n";
    Vector<double> vec2(sp);
    cout<<endl<<"5. construct ASCIIReader FO for file ut6.dat\n";
    ASCIIReader<double> ldr("./ut6.dat");
    cout<<endl<<"6. load ut6.aux into vector constructed in step 4\n"; 
    vec2.eval(ldr);
    cout<<endl<<"7. write both vectors to stdout: they should be same\n";
    cout<<"   also same as contents of ut6.dat, which you should check\n";
    cout<<endl<<"   first vector:\n"<<endl;
    vec1.write(cout);
    cout<<endl<<"   second vector:\n"<<endl;
    vec2.write(cout);
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 6              *"<<endl;
    cout<<" *************************************************/"<<endl; 
  }
  catch (RVLException & e) {
    e.write(cout);
    exit(1);
  }
}

  
