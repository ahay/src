#include "rsf.hh"
#include "cub.hh"
using namespace std;

//------------------------------------------------------------
CUB::CUB( const char* tag, const char* typ)
{
    if(tag) {
	if(typ == "i") {
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
    // read esize
    if (!sf_histint(file_,"esize",&e)) e = 4;
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
void CUB::getax(int id, axa *a)
{
    a->n = n[id];
    a->o = o[id];
    a->d = d[id];
}

void CUB::putax(int id, axa *a)
{
    n[id] = a->n;
    o[id] = a->o;
    d[id] = a->d;
}
//------------------------------------------------------------
