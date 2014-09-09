// identity.H
// created by ADP 12/23/03

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


#ifndef __RVL_FUNCTION_SETTOCOLUMNOFIDENTITY_
#define __RVL_FUNCTION_SETTOCOLUMNOFIDENTITY_

#include "data.hh"

namespace RVL {

  /** A function object used in the CompassSearch class of gps.H. 
      It sets the ith element to 1, and all the rest to zero.  

      6/17/04 Modified to work correctly when called multiple times.
              Will set the ith element seen to 1, and all others to zero.
	      Thus, still valid when i > n.
*/
template<class Scalar>
class SetToColumnOfIdentity : public UnaryFunctionObject<Scalar> {
private:
  long i;
  long seensofar;
public:
  /** the integer indicates the column of identity to which we set vectors
   */
  SetToColumnOfIdentity(long _i): i(_i), seensofar(0) {}

  /** never reads data
  */
   virtual bool readsData(int i=0) { return false; }

  /** always writes data      
  */
  virtual bool writesData(int i=0) { return true; }

  /** Name method */
  virtual string name() { return "SetToColumnOfIdentity";}

  /** report to exception */
  virtual void write(RVLException & e) {
    e << "SetToColumnOfIdentity with column = " << i << '\n';
  }

  /** report to ostream */
  virtual ostream & write(ostream & str) {
    str << "SetToColumnOfIdentity with column = " << i << '\n';
    return str;
  }

  /** Evaluation method */
  virtual void operator () (LocalDataContainer<Scalar> & A) {
    int n = A.getSize();
    Scalar * d = A.getData();
    for(int j = 0; j < n; j++ )
      d[j] = 0;

    if((i > seensofar)&&(i < n+seensofar ))
      d[i-seensofar] = 1;
  }
};

}
#endif
