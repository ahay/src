#ifndef CUB_H
#define CUB_H

#include <valarray>
#include <iostream>
#include <complex>
#include <rsf.h>
using namespace std;

class CUB {
public:
    // constructor
    CUB ();
    CUB (const char* tag, const char* typ);
    // destructor
    ~CUB();
		

    void init(const char* tag, const char* typ);
    
    void closedelete();

    void headin();       // input  header
    void headou();       // output header
    void clone(CUB c);   // clone a header from another header
    void report();       // report
    int  size();         // return cube size
    int  esize();        // return esize
    int  num_dimesions();// returns the number of dimensions
    void settype( sf_datatype type); // Set datatype

    sf_axis getax(int);      // get an axis from CUB
    void putax(int,sf_axis); // put an axis in   CUB
    void setup(int);     // setup cube dimensions 

    sf_axis *ax; // different axes
    int   nd; // number of dimensions

    // whence SEEK_SET=0 (start of file), SEEK_CUR=1 (current location), 
    // SEEK_END=2
    void  seek(off_t offset, int whence) const;	

    const CUB& operator>> (std::valarray <complex<float> > &vect) const;
    const CUB& operator<< (std::valarray <complex<float> > &vect) const;

    const CUB& operator>> (std::valarray <sf_complex> &vect) const;
    const CUB& operator<< (std::valarray <sf_complex> &vect) const;

    const CUB& operator>> (std::valarray <float> &vect) const;
    const CUB& operator<< (std::valarray <float> &vect) const;

    const CUB& operator>> (std::valarray <int> &vect) const;
    const CUB& operator<< (std::valarray <int> &vect) const;

    const CUB& operator>> (std::valarray <short> &vect) const;
    const CUB& operator<< (std::valarray <short> &vect) const;

    const CUB& operator>> (std::valarray <char> &vect) const;
    const CUB& operator<< (std::valarray <char> &vect) const;

    void  read(complex<float> &in,int num_values);
    void  write(complex<float> &in,int num_values);

    void  read(sf_complex &in,int num_values);
    void  write(sf_complex &in,int num_values);

    void  read(float &in,int num_values);
    void  write(float &in,int num_values);

    void  read(int &in,int num_values);
    void  write(int &in,int num_values);

    void  read(short &in,int num_values);
    void  write(short &in,int num_values);

    void  read(char &in,int num_values);
    void  write(char &in,int num_values);

private:
    sf_file file_; // file structure
};

#endif
