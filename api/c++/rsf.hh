// File: rsf.hh
// --------------
// RSF C++ interface
///////////////////////
#ifndef _rsf_hh
#define _rsf_hh

#include <valarray>

// to simplify main program includes
extern "C" {
#include <rsf.h>
}

class iRSF {
public:
    // constructor
    iRSF (const char* file = "in");
    // destructor
    ~iRSF();
    // file size
    int size (int dim=0);
    // reading data
    const iRSF& operator>> (std::valarray <float> &vect) const;
    // reading parameters
    void get (const char* name,   int& value,   int defolt) const;
    void get (const char* name,   int& value) const;
    void get (const char* name, float& value, float defolt) const;
    void get (const char* name, float& value) const;
    void get (const char* name, double& value, double defolt) const;
    void get (const char* name, double& value) const;
    void get (const char* name,  bool& value,  bool defolt) const;
    void get (const char* name,  bool& value) const;
//    void get (const char* name,  char* value, const char* defolt) const;
//    void get (const char* name,  char* value) const;
    // reading parameter arrays
    void get (const char* name, int size, int*   value, 
	      const int*   defolt) const;
    void get (const char* name, int size, int*   value) const;
    void get (const char* name, int size, float* value, 
 	      const float* defolt) const;
    void get (const char* name, int size, float* value) const;
    void get (const char* name, int size, bool* value, 
 	      const bool* defolt) const;
    void get (const char* name, int size, bool* value) const;
private:
    sf_file file_;
    // copy constructor - undefined to prevent misuse
    iRSF( const iRSF& src); 				
    // assignment operator - undefined to prevent misuse
    iRSF& operator=( const iRSF& src);
};

class oRSF {
public:
    // constructors 
    oRSF (const char* file = "out");
    // destructor 
    ~oRSF();
    // writing data
    const oRSF& operator<< (std::valarray <float> &vect) const;
    // writing parameters
    void put (const char* name, int value) const;
    void put (const char* name, float value) const;
    void put (const char* name, const char* value) const; 
    // writing parameter arrays
    void put (const char* name, int size, const int*   value) const;
//     void put (const char* name, int size, const float* value) const;
private:
    sf_file file_; 
    // copy constructor - undefined to prevent misuse
    oRSF( const oRSF& src); 				
    // assignment operator - undefined to prevent misuse
    oRSF& operator=( const oRSF& src); 
};

#endif

// 	$Id: rsf.hh 969 2005-01-21 03:20:13Z fomels $	
