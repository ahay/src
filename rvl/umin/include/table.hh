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

#ifndef TABLE_H
#define TABLE_H

#include "std_cpp_includes.hh"
#include "except.hh"
#include "utility.hh"

namespace RVL {

  // ASCII Parameter Table Object.
  //    This class implements the HCL associative array (hash table)
  //    interface defined in the {\bf AbstractTable} base class.

  //    The principal Table constructor reads the Table's internal data
  //    from a "config" file. Each line of the config file is either of
  //    the form "Key = Value" or "Object::Key = Value".  The first form
  //    initializes the parameter associated to "Key" in the table to the
  //    value "Value". The second form facilitates the use of one config
  //    file for several different objects. If the "Object" matches one of
  //    the object names provided by the user (as argument to the
  //    constructor, see below), then the parameter is read/written as
  //    before. Otherwise, the line is ignored.  The type of "Value" can
  //    be integer, float, double, or string. String values are optionally
  //    enclosed by double-quotes.


  class Table {

  private:

    typedef map<string,string> str_map;
    typedef str_map::iterator iterator;
    typedef str_map::const_iterator const_iterator;
    typedef str_map::value_type value_type;
  
    string prefixes;
    str_map table;

    static bool isPunct(char c);
    static string stripWord(const string word);
  
  public:

    // usual constructor
    Table(string fname="",string prefixes="");

    // Copy constructor 
    Table(const Table & T);

    // Destructor 
    virtual ~Table();

    // Merge entries from a file into the parameter table 
    int mergeFile(string fname);

    // Returns size of the table
    int getSize() const { return table.size(); }

    // Returns ith key (do we want this one?)
    string getKey(unsigned i);

    // Returns ith value (do we want this one?)
    string strValue(string key);

    // Methods to access values in the table 
    int getValue(string key,int & value) const;
    int getValue(string key,double & value) const;
    int getValue(string key,float & value) const;
    int getValue(string key,string & value,int length=0) const;
  
    // Methods to insert new pairs (key-value) in the table 
    void putValue(string key,int value);
    void putValue(string key,float value);
    void putValue(string key,double value);
    void putValue(string key,const string value);

    // idem (?do we really need those?)
    int getArrayValue(string key,int ind,int & value) const;
    int getArrayValue(string key,int ind,float & value) const;
    int getArrayValue(string key,int ind,double & value) const;
    int getArrayValue(string key,int ind,string & value,int length=0) const;

    // idem (?do we really need those?)
    void putArrayValue(string key,int ind,int value);
    void putArrayValue(string key,int ind,float value);
    void putArrayValue(string key,int ind,double value);
    void putArrayValue(string key,int ind,const string value);
  
    ostream & write(ostream & os);
  
  };

  ostream & operator << (ostream & os,Table & t);

  // convenient helper functions added 12.05 WWS
  // deprecated 04.07 WWS
  string getStringFromTable(Table const & par, string key);
  int getIntFromTable(Table const & par, string key);
  float getFloatFromTable(Table const & par, string key);
  double getDoubleFromTable(Table const & par, string key);

  template<typename Scalar>
  Scalar getValueFromTable(Table const & par, string key);

  template<>
  int getValueFromTable<int>(Table const & par, string key);

  template<>
  double getValueFromTable<double>(Table const & par, string key);

  template<>
  float getValueFromTable<float>(Table const & par, string key);

  template<>
  string getValueFromTable<string>(Table const & par, string key);
  
  template<>
  bool getValueFromTable<bool>(Table const & par, string key);
  
}

#endif


