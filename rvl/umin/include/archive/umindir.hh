
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

#ifndef __ALG_UMINDIR_H_
#define __ALG_UMINDIR_H_

#include "space.hh"
#include "functional.hh"
#include "lnsrch.hh"

namespace RVLUmin {

  using RVL::Vector;
  using RVL::FunctionalEvaluation;
  using RVLAlg::Terminator;

  /** Abstract interface for computation of search 
      directions, in support of line search methods.
      Since computation can succeed or fail, must have 
      character of terminator.
  */

  template<typename Scalar>
  class UMinDir: public Terminator {
  public:
    UMinDir() {}
    UMinDir(UMinDir<Scalar> const &) {}
    virtual ~UMinDir() {}

    // note that Terminator::query must be defined by instantiable subclass 
    /** Returns search direction in mutable first argument */
    virtual void calcDir(Vector<Scalar> & dir,
			 FunctionalEvaluation<Scalar> & fx) = 0;
    /** Use data generated during line search to update any internals */
    virtual void updateDir(LineSearchAlg<Scalar> const & ls) = 0;
    /** Reinitialize direction computation */
    virtual void resetDir() = 0;
    /** verbose output stream */
    
    virtual ostream & write(ostream & str) const = 0;
  };
}

#endif
