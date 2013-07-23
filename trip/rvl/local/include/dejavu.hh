/*************************************************************************

Copyright Rice University, 2004, 2005, 2006
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

#ifndef __RVL_DEJAVU
#define __RVL_DEJAVU

#include "std_cpp_includes.hh"
#include "except.hh"

namespace RVL {

  /** Helper function - checks whether vec[i] is identical to vec[j]
      for some j<i. If so, on return i changes to j, otherwise is
      unchanged. Typical use:
      
      std::vector<T *> vec;
      ...
      size_t j=i;
      deja-vu<T>(&j,vec);
      if (j==i) { do something for unaliased T }
      else { do something else for aliased T }
      
      Used to avoid aliasing in PackageContainer::eval methods. */
  
  template<class T>
  void dejavu(size_t * i, std::vector<T *> vec) {
    if (*i>vec.size()-1) { 
      RVLException e;
      e<<"Error: function dejavu\n";
      e<<"input int "<<*i<<" out of range [0, "<<vec.size()-1<<"].\n";
      throw e;
    }
    // save init value of *i as loop limit
    size_t jsave = *i;
    // in loop, reset *i if vec value is aliased to previous, AND
    // if it's still = saved value - the second condn prevents repeated
    // assignments
    for (size_t j=0;j<jsave;j++) { 
      if ((vec[j] == vec[*i]) &&
	  (*i == jsave)) { *i=j; }
    }
  }    

}

#endif
