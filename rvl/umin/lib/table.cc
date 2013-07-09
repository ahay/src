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

#include "table.hh"

namespace RVL {

  Table::Table(string fname,string prefix_str)
    : prefixes(prefix_str) {
    try {
      if (fname!="") {
	if (mergeFile(fname)) {
	  RVLException e; 
	  e<<"error: Table constructor from Table::Merge\n";
	  e<<"failed to merge data from file \""<<fname<< "\"";
	  throw e;
	}
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from Table::Table(string,string)\n";
      throw e;
    }
  }

  Table::Table(const Table & T) {
    prefixes = T.prefixes;
    table = T.table;
  }

  Table::~Table() { table.clear(); }

  string Table::getKey(unsigned i) {
    if (i>table.size()) {
      cerr<<"warning: Table::getKey: index greater than table size";
      return "";
    }
    else {
      const_iterator it = table.begin();
      for (unsigned iter=0;iter<=i;++iter) ++it;
      return it->first;
    }
  }

  string Table::strValue(string key)
  {
    const_iterator it(table.find(key)),end_it(table.end());
    if (it==end_it) {
      cerr<<"warning: Table::getValue: key \""
	  <<key<<"\" does not exist"<<endl;
      return "";
    }
    else
      return it->second;
  }

  int Table::getValue(string key,
		      int & value) const {
    const_iterator it(table.find(key)),end_it(table.end());
    if (it==end_it)
      return 1;
    istringstream inputString(it->second); 
    inputString>>value;
    return 0;
  }

  int Table::getValue(string key,
		      double & value) const {
    const_iterator it(table.find(key)),end_it(table.end());
    if (it==end_it)
      return 1;
    istringstream inputString(it->second); 
    inputString>>value;
    return 0;
  }

  int Table::getValue(string key,
		      float & value) const { 
    const_iterator it(table.find(key)),end_it(table.end());
    if (it==end_it)
      return 1;
    istringstream inputString(it->second); 
    inputString>>value;
    return 0;
  }

  int Table::getValue(string key,
		      string & value,
		      int length) const {
    const_iterator it(table.find(key)),end_it(table.end());
    if (it==end_it) 
      return 1;
    if (length>0) {
      int len=(it->second).length();
      if (length<len) {
	cout<<"warning: Table::getValue: truncating string for key \""
	    <<key <<'\"'<<endl;
	value.assign(it->second,0,length);
      }
      else
	value.assign(it->second);
      return 0;
    }
    value.assign(it->second);
    return 0;
  }

  void Table::putValue(string key,
		       int value) {
    ostringstream outputString; 
    outputString<<value;
    // if (outputString.fail()) throw exception?
    //table.insert(value_type(key,outputString.str()));
    table[key]=outputString.str();
  } 

  void Table::putValue(string key,
		       double value) {
    ostringstream outputString; 
    outputString<<value;
    //table.insert(value_type(key,outputString.str()));
    table[key]=outputString.str();
  }

  void Table::putValue(string key,
		       float value) {
    ostringstream outputString; 
    outputString<<value;
    //table.insert(value_type(key,outputString.str()));
    table[key]=outputString.str();
  }

  void Table::putValue(string key,
		       const string value) {
    ostringstream outputString(value);
    //table.insert(value_type(key,outputString.str()));
    table[key]=outputString.str();
  } 

  int Table::getArrayValue(string key,
			   int ind,
			   int & value) const {
    ostringstream outputString;
    outputString<<ind;
    string newkey=key+outputString.str();
    return getValue(newkey,value);
  } 

  int Table::getArrayValue(string key,
			   int ind,
			   float & value) const {
    ostringstream outputString;
    outputString<<ind;
    string newkey=key+outputString.str();
    return getValue(newkey,value);  
  } 

  int Table::getArrayValue(string key,
			   int ind,
			   double & value) const { 
    ostringstream outputString;
    outputString<<ind;
    string newkey=key+outputString.str();
    return getValue(newkey,value); 
  } 

  int Table::getArrayValue(string key,
			   int ind,
			   string & value,
			   int length) const { 
    ostringstream outputString;
    outputString<<ind;
    string newkey=key+outputString.str();
    return getValue(newkey,value,length);
  }


  void Table::putArrayValue(string key,
			    int ind,
			    int value) {
    ostringstream outputString;
    outputString<<ind;
    string newkey=key+outputString.str();
    putValue(newkey,value);
  }

  void Table::putArrayValue(string key,
			    int ind,
			    float value) {
    ostringstream outputString;
    outputString<<ind;
    string newkey=key+outputString.str();
    putValue(newkey,value);
  }

  void Table::putArrayValue(string key,
			    int ind,
			    double value) {
    ostringstream outputString;
    outputString<<ind;
    string newkey=key+outputString.str();
    putValue(newkey,value);
  }

  void Table::putArrayValue(string key,
			    int ind,
			    const string value) {
    ostringstream outputString;
    outputString<<ind;
    string newkey=key+outputString.str();  
    putValue(newkey,value);
  }

  int Table::mergeFile(string fname) {

    try {
      ifstream inFile(fname.c_str(),ios::in);
      string line,subline;
      unsigned line_cnt = 0;

      if (!inFile)
	return 1;

      while (getline(inFile,line) && !inFile.bad()) {

	++line_cnt;
	if (line.length()<string::npos && line.find("=")!=string::npos) {  
	  // Ignore comments 
	  string::size_type x = line.find("//");
	  if (x<line.length()) line.erase(x);
	  // Remove all spaces between words
	  string::size_type y =line.find(' ');
	  while (y < line.length()) {
	    line.erase(y,1);
	    y=line.find(' ',y+1);
	  }
	  // Check whether the table is constructed for a particular class
	  if (prefixes.empty()) {    
	    // Read in only "general-purpose" parameters
	    if (line.find("::")>line.length()) {
	      unsigned cur = line.find_first_of("=");
	      subline = line.substr(0,cur);
	      line = line.substr(cur+1,line.length());
	      table.insert(value_type(subline,stripWord(line)));
	    }
	  }
	  else {
	    // Read in only those parameters pertaining to the object
	    // identified by the prefix
	    if (line.find("::")!=string::npos) {
	      unsigned pos = line.find_first_of(":");
	      subline = line.substr(0,pos);
	      if (prefixes.find(subline)!=string::npos) {
		// prefixes.resize(prefixes.length()+subline.length());
		// prefixes.append(subline); prefixes.append(" ");
		unsigned cur = line.find_first_of("=");
		subline = line.substr(pos+2, cur-pos-2);
		line = line.substr(cur+1, line.length());
		table.insert(value_type(subline,stripWord(line)));
	      }
	    }
	  }
	}
	else if (!line.empty()) {
	  cerr<<"warning: Table::merge skipped line #"<<line_cnt
	      <<" from input file \""<<fname << "\""<<endl;
	}
      }
      return 0;
    }
    catch (RVLException & e) {
      e<<"\ncalled from Table::mergeFile\n";
      throw e;
    }
    catch (...) {
      RVLException e;
      e<<"\ncalled from Table::mergeFile\n";
      throw e;
    }
  }
  
  ostream & Table::write(ostream & os) {

    iterator it=table.begin(),end_it=table.end();
    ostringstream outputString(prefixes);
  
    os<<"Table ( "<<outputString.str();
    os<<" )"<<endl<<"{"<<endl;
    for (;it!=end_it;++it) {
      os<<"  "<<it->first<<" = "<<it->second<<endl;
    }
    os<<"}"<<endl;
    return os;
  }

  ostream & operator<<(ostream & os,Table & T) {
    return T.write(os);
  }

  string Table::stripWord(const string word) {
    // get first index into word which is not a punct.
    int first = 0;
    while ((isPunct(word[first])) && (first<int(word.length())-1))
      first++;
    // get first index from the right into word without a punct.
    int last=word.length()-1; 
    while ((isPunct(word[last])) && (last>0))
      last--;
  
    if (first>last)
      return( "" );
    else
      return (word.substr(first,last-first+1));
  }

  bool Table::isPunct(char c) {
    char punct[] = {' ',',',';',':','!','(',')','"','\''};
    const int num_of_punct = sizeof(punct)/sizeof(char);
    for (int i=0;i < num_of_punct;i++)
      if (c==punct[i])
	return (true);
    return false;
  }

  template<typename Scalar>
  Scalar getValueFromTable(Table const & tbl, string key) {
    RVLException e;
    e<<"Error: getValueFromTable\n";
    e<<"unknown type";
    throw e;
  }

  template<>
  int getValueFromTable<int>(Table const & tbl, string key) {
    int t;
    if (tbl.getValue(key,t)) {
      RVLException e;
      e<<"Error: getValueFromTable (int) key = "<<key<<"\n";
      throw e;
    }
    return t;
  }

  template<>
  double getValueFromTable<double>(Table const & tbl, string key) {
    double t;
    if (tbl.getValue(key,t)) {
      RVLException e;
      e<<"Error: getValueFromTable (double) key = "<<key<<"\n";
      throw e;
    }
    return t;
  }

  template<>
  float getValueFromTable<float>(Table const & tbl, string key) {
    float t;
    if (tbl.getValue(key,t)) {
      RVLException e;
      e<<"Error: getValueFromTable (float) key = "<<key<<"\n";
      throw e;
    }
    return t;
  }
  
  template<>
  string getValueFromTable<string>(Table const & tbl, string key) {
    string t;
    if (tbl.getValue(key,t)) {
      RVLException e;
      e<<"Error: getValueFromTable (string) key = "<<key<<"\n";
      throw e;
    }
    return t;
  }
  
  template<>
  bool getValueFromTable<bool>(Table const & tbl, string key) {
    int t;
    if (tbl.getValue(key,t)) {
      RVLException e;
      e<<"Error: getValueFromTable (int) key = "<<key<<"\n";
      throw e;
    }
    if (t!=0) return true;
    return false;
  }


  // the following functions are deprecated - 04.07 

  string getStringFromTable(Table const & tbl, string key) {
    string val;
    if (tbl.getValue(key,val)) {
      RVLException e;
      e<<"Error: getStringFromTable helper function\n";
      e<<"failed to extract value for key = "<<key<<"\n";
      throw e;
    }
    return val;
  }

  int getIntFromTable(Table const & tbl, string key) {
    int val;
    if (tbl.getValue(key,val)) {
      RVLException e;
      e<<"Error: getIntFromTable helper function\n";
      e<<"failed to extract value for key = "<<key<<"\n";
      throw e;
    }
    return val;
  }

  float getFloatFromTable(Table const & tbl, string key) {
    float val;
    if (tbl.getValue(key,val)) {
      RVLException e;
      e<<"Error: getFloatFromTable helper function\n";
      e<<"failed to extract value for key = "<<key<<"\n";
      throw e;
    }
    return val;
  }

  double getDoubleFromTable(Table const & tbl, string key) {
    double val;
    if (tbl.getValue(key,val)) {
      RVLException e;
      e<<"Error: getDoubleFromTable helper function\n";
      e<<"failed to extract value for key = "<<key<<"\n";
      throw e;
    }
    return val;
  }


}
