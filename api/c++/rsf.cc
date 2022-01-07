// File: rsf.cc
// --------------
// RSF C++ interface
/////////////////////////////////////////////////

#include <vector>
#include <string>
using std::string;

#include "rsf.hh"

// ***** Public methods ***********
///////////////////////////////////

// iRSF Constructor
// ----------------
// opens and RSF file
iRSF::iRSF( const char* file /* = "in" */)
{
    if (file) { 
	file_ = sf_input(file);
    } else { // NULL means "command line parameters"
	file_ = 0;
    }
}

// oRSF Constructor
// ----------------
// Copies the file tag, writes out the grid
oRSF::oRSF(const char* file /* = "out" */)
{
    if (NULL == file) sf_error("cannot initialize oRSF with NULL"); /* yet */
    file_ = sf_output(file);
}

// iRSF Destructor
// ----------------
iRSF::~iRSF()
{
    if (file_) sf_fileclose(file_);
}

// oRSF Destructor
// ----------------
oRSF::~oRSF()
{
    sf_fileclose(file_);
}

// flush header
// ----------------
void oRSF::flush()
{
    sf_fileflush(file_,NULL);
}


// file size
// ----------------
int 
iRSF::size (int dim /* =0 */)
{
    return sf_leftsize(file_,dim);
}

// file data type
int 
iRSF::type (void)
{
    return (int) sf_gettype(file_);
}

// Reading data
// ------------
const iRSF&
iRSF::operator>> (float &value) const
{
    if (NULL == file_) sf_error("Cannot write data to parameter file");
    sf_floatread(&value,1,file_);
    return *this;
}

const iRSF&
iRSF::operator>> (std::valarray<float> &array) const
{
    if (NULL == file_) sf_error("Cannot write data to parameter file");
    sf_floatread(&(array[0]),array.size(),file_);
    return *this;
}

const iRSF&
iRSF::operator>> (int &value) const
{
    if (NULL == file_) sf_error("Cannot write data to parameter file");
    sf_intread(&value,1,file_);
    return *this;
}

const iRSF&
iRSF::operator>> (std::valarray<int> &array) const
{
    if (NULL == file_) sf_error("Cannot write data to parameter file");
    sf_intread(&(array[0]),array.size(),file_);
    return *this;
}

const iRSF&
iRSF::operator>> (sf_complex &value) const
{
    if (NULL == file_) sf_error("Cannot write data to parameter file");
    sf_complexread(&value,1,file_);
    return *this;
}

const iRSF&
iRSF::operator>> (std::valarray<sf_complex> &array) const
{
    if (NULL == file_) sf_error("Cannot write data to parameter file");
    sf_complexread(&(array[0]),array.size(),file_);
    return *this;
}

void
iRSF::unpipe(off_t size) // prepare for direct access
{
    sf_unpipe(file_,size);
}

off_t 
iRSF::tell(void) // position in a file
{
    return sf_tell(file_);
}

void
iRSF::seek( off_t offset, int whence) // seek to a position in a file
{
    sf_seek(file_,offset, whence);
}

// Writing data
// ------------
const oRSF&
oRSF::operator<< (float &value) const
{
    sf_floatwrite(&value,1,file_);
    return *this;
}

const oRSF&
oRSF::operator<< (std::valarray<float> &array) const
{
    sf_floatwrite(&(array[0]),array.size(),file_);
    return *this;
}

const oRSF&
oRSF::operator<< (int &value) const
{
    sf_intwrite(&value,1,file_);
    return *this;
}

const oRSF&
oRSF::operator<< (std::valarray<int> &array) const
{
    sf_intwrite(&(array[0]),array.size(),file_);
    return *this;
}

const oRSF&
oRSF::operator<< (sf_complex &value) const
{
    sf_complexwrite(&value,1,file_);
    return *this;
}

const oRSF&
oRSF::operator<< (std::valarray<sf_complex> &array) const
{
    sf_complexwrite(&(array[0]),array.size(),file_);
    return *this;
}

off_t 
oRSF::tell(void) // position in a file
{
    return sf_tell(file_);
}

void
oRSF::seek( off_t offset, int whence) // seek to a position in a file
{
    sf_seek(file_,offset, whence);
}

// file data type
int 
oRSF::bufsiz (void)
{
    return (int) sf_bufsiz(file_);
}

// Reading parameters
/////////////////////
void 
iRSF::get (const char* name, int &value, int defolt) const
{
    if (file_)  {
	if (!sf_histint (file_,name,&value)) value = defolt;
    } else {
	if (!sf_getint (name, &value)) value = defolt;
    }
}

void 
iRSF::get (const char* name, int& value) const
{
    if (file_)  {
	if (!sf_histint(file_,name,&value)) 
	    sf_error("missing history value: %s",name);
    } else {
	if (!sf_getint(name,&value)) 
	    sf_error("missing parameter value: %s",name);
    }
}

void 
iRSF::get (const char* name, bool &value, bool defolt) const
{
    bool cvalue[4];

    if (file_)  {
	if (sf_histbool (file_,name,cvalue)) value=cvalue[0]; 
	else value = defolt;
    } else {
	if (sf_getbool (name,cvalue)) value=cvalue[0]; 
	else value = defolt;
    }
}

void 
iRSF::get (const char* name, bool& value) const
{
    bool cvalue[4];	
	
    if (file_)  {
	if (!sf_histbool(file_,name,cvalue)) 
	    sf_error("missing history value: %s",name);
    } else {
	if (!sf_getbool(name,cvalue)) 
	    sf_error("missing parameter value: %s",name);
    }

    value = cvalue[0];
}

void 
iRSF::get (const char* name, float& value, float defolt) const
{
    if (file_)  {
	if (!sf_histfloat (file_,name,&value)) value = defolt;
    } else {
	if (!sf_getfloat (name,&value)) value = defolt;
    }
}

void 
iRSF::get (const char* name, double& value, double defolt) const
{
    if (file_)  {
	if (!sf_histdouble (file_,name,&value)) value = defolt;
    } else {
	if (!sf_getdouble (name,&value)) value = defolt;
    }
}

void 
iRSF::get (const char* name, float& value) const
{
    if (file_)  {
	if (!sf_histfloat(file_,name,&value)) 
	    sf_error("missing history value: %s",name);
    } else {
	if (!sf_getfloat(name,&value)) 
	    sf_error("missing parameter value: %s",name);
    }
}

void 
iRSF::get (const char* name, double& value) const
{
    if (file_)  {
	if (!sf_histdouble(file_,name,&value)) 
	    sf_error("missing history value: %s",name);
    } else {
	if (!sf_getdouble(name,&value)) 
	    sf_error("missing parameter value: %s",name);
    }
}

void
iRSF::get (const char* name, string& value, const string& defolt) const
{
    char* retval;
    if (file_)
        retval = sf_histstring( file_, name );
    else
        retval = sf_getstring( name );
    if( retval ) {
        value = retval;
        free( retval ); // retval was allocated with malloc and we now own it!
    }
    else
        value = defolt;
}

void
iRSF::get (const char* name, string& value) const
{
    char* retval;
    if (file_) {
        retval = sf_histstring( file_, name );
        if( !retval )
                sf_error("missing history value: %s",name);
    } else {
        retval = sf_getstring( name );
        if( !retval )
                sf_error("missing parameter value: %s",name);
    }
    value = retval;
    free( retval ); // retval was allocated with malloc and we now own it!
}

// Reading parameter arrays
///////////////////////////
void 
iRSF::get (const char* name, int size, int*   value, const int*   defolt) const
{
    if (file_) {
	if (!sf_histints(file_,name, value,size)) {
	    for (int i = 0; i < size; i++) {
		value[i] = defolt[i];
	    }
	}
    } else {
	if (!sf_getints(name, value,size)) {
	    for (int i = 0; i < size; i++) {
		value[i] = defolt[i];
	    }
	}
    }
}

void 
iRSF::get (const char* name, int size, int*   value) const
{
    if (file_) {
	if (!sf_histints(file_,name, value,size))
	    sf_error("missing history value: %s",name);
    } else {
	if (!sf_getints(name, value,size)) 
	    sf_error("missing parameter value: %s",name);
    }
}

void 
iRSF::get (const char* name, int size, float* value, const float* defolt) const
{
    if (file_) {
	if (!sf_histfloats(file_,name,value,size)) {
	    for (int i = 0; i < size; i++) {
		value[i] = defolt[i];
	    }
	}
    } else {
	if (!sf_getfloats(name, value,size)) {
	    for (int i = 0; i < size; i++) {
		value[i] = defolt[i];
	    }
	}
    }
}

void 
iRSF::get (const char* name, int size, float* value) const
{
    if (file_) {
	if (!sf_histfloats(file_,name,value,size))
	    sf_error("missing history value: %s",name);
    } else {
	if (!sf_getfloats(name, value,size)) 
	    sf_error("missing parameter value: %s",name);
    }
}

void
iRSF::get (const char* name, int size, bool* value, const bool* defolt) const
{
    //// BUG: in gcc < 3 where sizeof(bool)=4

    if (file_) {
	if (!sf_histbools(file_,name,value,size)) {
	    for (int i = 0; i < size; i++) {
		value[i] = defolt[i];
	    }
	}
    } else {
	if (!sf_getbools(name, value,size)) {
	    for (int i = 0; i < size; i++) {
		value[i] = defolt[i];
	    }
	}
    }
}

void 
iRSF::get (const char* name, int size, bool* value) const
{
    //// BUG: in gcc < 3 where sizeof(bool)=4

    if (file_) {
	if (!sf_histbools(file_,name,value,size))
	    sf_error("missing history value: %s",name);
    } else {
	if (!sf_getbools(name, value,size)) 
	    sf_error("missing parameter value: %s",name);
    }
}

// Writing parameters
/////////////////////

// set file data type
void
oRSF::type (sf_datatype type)
{
    sf_settype(file_,type);
}


void 
oRSF::put (const char* name, int value) const
{
    sf_putint(file_,name,value);
}
 
void 
oRSF::put (const char* name, float value) const
{
    sf_putfloat(file_,name,value);
}

void 
oRSF::put (const char* name, const char* value) const
{
    sf_putstring(file_,name,value);
}

// Writing parameter arrays
///////////////////////////
void 
oRSF::put (const char* name, int size, const int*   value) const
{ 
    sf_putints(file_,name,value,size);
}

// void 
// oRSF::put (const char* name, int size, const float* value) const
// {
//     sf_putfloats(file_,name,value,size);
// }

// 	$Id: rsf.cc 969 2005-01-21 03:20:13Z fomels $	
