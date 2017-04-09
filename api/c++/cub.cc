#include<complex>
#include "rsf.hh"
#include "cub.hh"
using namespace std;

//------------------------------------------------------------
CUB::CUB( )
{
    file_ = 0;
}
//------------------------------------------------------------
//------------------------------------------------------------
CUB::CUB( const char* tag, const char* typ)
{
    if(tag) {
	if(typ[0] == 'i') {
	    file_ = sf_input (tag);
	} else {
	    file_ = sf_output(tag);
	}
    } else {
        file_ = 0;
    }
}
//------------------------------------------------------------
CUB::~CUB()
{
    sf_fileclose(file_);
}
//------------------------------------------------------------
void CUB::init( const char* tag, const char* typ)
{ 
    if(tag) {
	if(typ[0] == 'i') {
	    file_ = sf_input (tag);
	} else {
	    file_ = sf_output(tag);
	}
    } else {
        file_ = 0;
    }
}
//------------------------------------------------------------
void CUB::closedelete()
{
	sf_fileclosedelete(file_);
}
//------------------------------------------------------------
void CUB::headin()
{
    int nn[SF_MAX_DIM];
    char key[3];

    nd = sf_filedims(file_,nn);
    n=sf_intalloc  (nd);
    o=sf_floatalloc(nd);
    d=sf_floatalloc(nd);

    // loop over dimensions
    // read axes
    for(int id=0;id<nd;id++) {	
	(void) snprintf(key,3,"n%d",id+1);
	if (!sf_histint  (file_,key,&n[id])) n[id] = 1;
	(void) snprintf(key,3,"o%d",id+1);
	if (!sf_histfloat(file_,key,&o[id])) o[id] = 0;
	(void) snprintf(key,3,"d%d",id+1);
	if (!sf_histfloat(file_,key,&d[id])) d[id] = 1;
    }
    // find esize
    e = sf_esize(file_);
}
//------------------------------------------------------------
void CUB::setup(int kd, int esize)
{
    nd=kd;
    n=sf_intalloc  (nd);
    o=sf_floatalloc(nd);
    d=sf_floatalloc(nd);

    e=esize;
}
//------------------------------------------------------------
void CUB::clone(CUB c)
{
    nd = c.nd;
    n=sf_intalloc  (nd);
    o=sf_floatalloc(nd);
    d=sf_floatalloc(nd);

    // loop over dimensions
    // copy axes
    for (int id=0;id<nd;id++) {
	n[id] = c.n[id];
	o[id] = c.o[id];
	d[id] = c.d[id];
    }
    // copy esize
    e = c.e;
}
//------------------------------------------------------------
void CUB::headou()
{
    char key[3];

    // loop over dimensions
    // write axes
    for(int id=0;id<nd;id++) {
	(void) snprintf(key,3,"n%d",id+1);
	sf_putint  (file_,key,n[id]);
	(void) snprintf(key,3,"o%d",id+1);
	sf_putfloat(file_,key,o[id]);
	(void) snprintf(key,3,"d%d",id+1);
	sf_putfloat(file_,key,d[id]);
    }
    // write esize
    sf_putint(file_,"esize",e);
}
//------------------------------------------------------------
void CUB::report()
{
    std::cerr << " REPORT " << endl;
    std::cerr << "nd=" << nd << " esize=" << e << endl;
    for(int id=0;id<nd;id++) {
	cerr << n[id] 
	     << " " 
	     << o[id] 
	     << " " 
	     << d[id] << endl;
    }
}
//------------------------------------------------------------
int CUB::size()
{
    int nn=1;
    for(int id=0;id<nd;id++) {
	nn*=n[id];
    }
    return(nn);
}
int CUB::esize()
{
    return(e);
}
//------------------------------------------------------------
int  CUB::num_dimesions()
{
    return(nd);
}
//------------------------------------------------------------
void CUB::settype(sf_datatype type)
{
  sf_settype(file_,  type); 
}
//------------------------------------------------------------
void CUB::seek(off_t offset, int whence) const
{
    sf_seek(file_,offset,whence);
}
//------------------------------------------------------------
const CUB&
CUB::operator>> (std::valarray<complex<float> > &array) const
{
    sf_complexread((reinterpret_cast<sf_complex*>(&array[0])), array.size(),file_);
    return *this;
}
//------------------------------------------------------------
const CUB&
CUB::operator<< (std::valarray<complex<float> > &array) const
{
    sf_complexwrite((reinterpret_cast<sf_complex*>(&array[0])),array.size(),file_);
    return *this;
}
//------------------------------------------------------------
const CUB&
CUB::operator>> (std::valarray<sf_complex> &array) const
{
    sf_complexread(&(array[0]),array.size(),file_);
    return *this;
}
//------------------------------------------------------------
const CUB&
CUB::operator<< (std::valarray<sf_complex> &array) const
{
    sf_complexwrite(&(array[0]),array.size(),file_);
    return *this;
}
//------------------------------------------------------------
const CUB&
CUB::operator>> (std::valarray<float> &array) const
{
    sf_floatread(&(array[0]),array.size(),file_);
    return *this;
}
//------------------------------------------------------------
const CUB&
CUB::operator<< (std::valarray<float> &array) const
{
    sf_floatwrite(&(array[0]),array.size(),file_);
    return *this;
}
//------------------------------------------------------------
const CUB&
CUB::operator>> (std::valarray<int> &array) const
{
    sf_intread(&(array[0]),array.size(),file_);
    return *this;
}
//------------------------------------------------------------
const CUB&
CUB::operator<< (std::valarray<int> &array) const
{
    sf_intwrite(&(array[0]),array.size(),file_);
    return *this;
}
//------------------------------------------------------------
const CUB&
CUB::operator>> (std::valarray<short> &array) const
{
    sf_shortread(&(array[0]),array.size(),file_);
    return *this;
}
//------------------------------------------------------------
const CUB&
CUB::operator<< (std::valarray<short> &array) const
{
    sf_shortwrite(&(array[0]),array.size(),file_);
    return *this;
}
//------------------------------------------------------------
const CUB&
CUB::operator>> (std::valarray<char> &array) const
{
    sf_charread(&(array[0]),array.size(),file_);
    return *this;
}
//------------------------------------------------------------
const CUB&
CUB::operator<< (std::valarray<char> &array) const
{
    sf_charwrite(&(array[0]),array.size(),file_);
    return *this;
}                                                                                                                      
//------------------------------------------------------------
void CUB::read (complex<float> &in,int num_values)
{
    sf_complexread(reinterpret_cast<sf_complex*>(&in),num_values,file_);
}
//------------------------------------------------------------
void CUB::write (complex<float> &in,int num_values)
{
    sf_complexwrite((reinterpret_cast<sf_complex*>(&in)),num_values,file_);
}

//------------------------------------------------------------
void CUB::read (sf_complex &in,int num_values)
{
    sf_complexread(&in,num_values,file_);
}
//------------------------------------------------------------
void CUB::write (sf_complex &in,int num_values)
{
    sf_complexwrite(&in,num_values,file_);
}
//------------------------------------------------------------
void CUB::read (float &in,int num_values)
{
    sf_floatread(&in,num_values,file_);
}
//------------------------------------------------------------
void CUB::write (float &in,int num_values)
{
    sf_floatwrite(&in,num_values,file_);
}
//------------------------------------------------------------
void CUB::read (int &in,int num_values)
{
    sf_intread(&in,num_values,file_);
}
//------------------------------------------------------------
void CUB::write (int &in,int num_values)
{
    sf_intwrite(&in,num_values,file_);
}
//------------------------------------------------------------
void CUB::read (short &in,int num_values)
{
    sf_shortread(&in,num_values,file_);
}
//------------------------------------------------------------
void CUB::write (short &in,int num_values)
{
    sf_shortwrite(&in,num_values,file_);
}//------------------------------------------------------------
void CUB::read (char &in,int num_values)
{
    sf_charread(&in,num_values,file_);
}
//------------------------------------------------------------
void CUB::write (char &in,int num_values)
{
    sf_charwrite(&in,num_values,file_);
}


//------------------------------------------------------------
sf_axis CUB::getax(int id)
{
    return sf_maxa(n[id],o[id],d[id]);
}

void CUB::putax(int id, sf_axis a)
{
    n[id] = sf_n(a);
    o[id] = sf_o(a);
    d[id] = sf_d(a);
}
//------------------------------------------------------------
