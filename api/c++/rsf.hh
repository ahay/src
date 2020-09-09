// File: rsf.hh
// --------------
// RSF C++ interface
///////////////////////
#ifndef _rsf_hh
#define _rsf_hh

#include <valarray>
#include <string>

#include <complex>
typedef std::complex<float> Complex8;
typedef std::complex<double> Complex16;

// to simplify main program includes
extern "C" {
#define cpx8 Complex8
#define cpx16 Complex16

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
    // file data type
    int type (void);
    // reading data
    const iRSF& operator>> (float &value) const;
    const iRSF& operator>> (int &value) const;
    const iRSF& operator>> (sf_complex &value) const;
    const iRSF& operator>> (std::valarray <float> &vect) const;
    const iRSF& operator>> (std::valarray <int> &vect) const;
    const iRSF& operator>> (std::valarray <sf_complex> &vect) const;
    // direct access
    void unpipe(off_t size);
    off_t tell(void);
    void seek(off_t offset, int whence);
    // reading parameters
    void get (const char* name,   int& value,   int defolt) const;
    void get (const char* name,   int& value) const;
    void get (const char* name, float& value, float defolt) const;
    void get (const char* name, float& value) const;
    void get (const char* name, double& value, double defolt) const;
    void get (const char* name, double& value) const;
    void get (const char* name,  bool& value,  bool defolt) const;
    void get (const char* name,  bool& value) const;
    void get (const char* name,  std::string& value, const std::string& defolt) const;
    void get (const char* name,  std::string& value) const;
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
    const oRSF& operator<< (float &value) const;
    const oRSF& operator<< (int &value) const;
    const oRSF& operator<< (sf_complex &value) const;
    const oRSF& operator<< (std::valarray <float> &vect) const;
    const oRSF& operator<< (std::valarray <int> &vect) const;
    const oRSF& operator<< (std::valarray <sf_complex> &vect) const;
    off_t tell(void);
    void seek(off_t offset, int whence);
    // set file data type
    void type (sf_datatype type);
    // writing parameters
    void put (const char* name, int value) const;
    void put (const char* name, float value) const;
    void put (const char* name, const char* value) const; 
    // writing parameter arrays
    void put (const char* name, int size, const int*   value) const;
//     void put (const char* name, int size, const float* value) const;
    void flush(void);
    int bufsiz(void);
private:
    sf_file file_; 
    // copy constructor - undefined to prevent misuse
    oRSF( const oRSF& src); 				
    // assignment operator - undefined to prevent misuse
    oRSF& operator=( const oRSF& src); 
};

#endif

// 	$Id: rsf.hh 969 2005-01-21 03:20:13Z fomels $	
