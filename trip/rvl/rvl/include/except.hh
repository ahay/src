/*************************************************************************

Copyright Rice University, 2004, 2005, 2006.
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

#ifndef __RVL_EXCEPT
#define __RVL_EXCEPT

#define BUFLEN 100

#include "std_cpp_includes.hh"

namespace RVL {

  /** An implementation of the std::exception interface, with additional
      methods so it can be used more like a output stream.
  */

  class RVLException: public std::exception {
  private:
    string msg;
  public:
    RVLException(): msg("") {}
    RVLException(const RVLException & s): msg("") { msg +=s.msg; }
    virtual ~RVLException() throw() {}

    const char* what() const throw() { return msg.c_str(); }

    RVLException & operator<< ( string str ) { 
      msg += str;
      return *this;
    }
    RVLException & operator<< ( const char* str ) { 
      msg += str;
      return *this;
    }
    RVLException & operator<< ( int i ) {
      char buf[ BUFLEN ];
      sprintf( buf, "%d", i );
      msg += buf;
      return *this;
    }
    RVLException & operator<< ( unsigned int i ) {
      char buf[ BUFLEN ];
      sprintf( buf, "%u", i );
      msg += buf;
      return *this;
    }
    RVLException & operator<< ( long i ) {
      char buf[ BUFLEN ];
      sprintf( buf, "%ld", i );
      msg += buf;
      return *this;
    }
    RVLException & operator<< ( unsigned long i ) {
      char buf[ BUFLEN ];
      sprintf( buf, "%lu", i );
      msg += buf;
      return *this;
    }
    RVLException & operator<< ( short i ) {
      char buf[ BUFLEN ];
      sprintf( buf, "%d", i );
      msg += buf;
      return *this;
    }
    RVLException & operator<< ( unsigned short i ) {
      char buf[ BUFLEN ];
      sprintf( buf, "%d", i );
      msg += buf;
      return *this;
    }
    /*
    RVLException & operator<< ( off_t i ) {
      char buf[ BUFLEN ];
      sprintf( buf, "%zd", i );
      msg += buf;
      return *this;
    }
    */
    RVLException & operator<< ( double d ) {
      char buf[ BUFLEN ];
      sprintf( buf, "%g", d );
      msg += buf;
      return *this;
    }
    RVLException & operator<< ( float d ) {
      char buf[ BUFLEN ];
      sprintf( buf, "%g", d );
      msg += buf;
      return *this;
    }
    template<class T>
    RVLException & operator<< ( complex<T> d ) {
      char buf[ BUFLEN ];
      sprintf( buf, "(%g,%g)", d.real(),d.imag() );
      msg += buf;
      return *this;
    }

    RVLException & operator<< ( char c ) {
      char buf[ BUFLEN ];
      buf[ 0 ] = c;
      buf[ 1 ] = '\0';
      msg += buf;
      return *this;
    }

    ostream & write(ostream & str) const {
      str<<msg<<endl;
      return str;
    }
  };

}

#endif
