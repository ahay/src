#include <valarray>
#include <iostream>
#include <rsf.h>
using namespace std;

class CUB {
public:
    // constructor
    CUB (const char* tag, const char* typ);
    // destructor
    ~CUB();

    void headin();       // input  header
    void headou();       // output header
    void clone(CUB c);   // clone a header from another header
    void report();       // report
    int  size();         // return cube size
    int  esize();        // return esize

    sf_axis getax(int);      // get an axis from CUB
    void putax(int,sf_axis); // put an axis in   CUB
    void setup(int,int);     // setup cube dimensions

    int   *n; // n for cube
    float *o; // o for cube
    float *d; // d for cube
    int    e; // esize
    int   nd; // number of dimensions

    const CUB& operator>> (std::valarray <float> &vect) const;
    const CUB& operator<< (std::valarray <float> &vect) const;
private:
    sf_file file_; // file structure
};

