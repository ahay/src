// File: rsf.cc
// --------------
// RSF C++ interface
/////////////////////////////////////////////////

#include <vector>

#include "rsf.hh"

#include "_bool.h"

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


// file size
// ----------------
int 
iRSF::size (int dim /* =0 */)
{
    return sf_leftsize(file_,dim);
}

// Reading data
// ------------
const iRSF&
iRSF::operator>> (std::valarray<float> &array) const
{
    if (NULL == file_) sf_error("Cannot write data to parameter file");
    sf_floatread(&(array[0]),array.size(),file_);
    return *this;
}

// Writing data
// ------------
const oRSF&
oRSF::operator<< (std::valarray<float> &array) const
{
    sf_floatwrite(&(array[0]),array.size(),file_);
    return *this;
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

// void 
// iRSF::get (const char* name, char* value, const char* defolt) const
// {
//     if (file_)  {
// 	if (0 == auxpar (name,"s", value, file_)) strcpy (value, defolt);
// 	if (0 != strcmp (file_, "in")) {
// 	    ostrstream buff (Buffer, MAXLEN);
// 	    buff << "From aux(" << file_  << "): " << file_ << "_" << name;
// 	    buff << '\0';
// 	    putch (Buffer, "s", value);
// 	}
//     } else {
// 	if (0 == getch (name, "s", value)) strcpy (value, defolt);
// 	ostrstream buff (Buffer, MAXLEN);
// 	buff << "From par: " << name<< '\0';
// 	putch (Buffer, "s", value);
//     }
// }

// void 
// iRSF::get (const char* name, char* value) const
// {
//     if (file_) {
// 	if (0 == auxpar(name,"s", value, file_)) {
// 	    cerr << "missing history value: " << name;
// 	    cerr << "from file: " << file_ << "\n";
// 	    exit (-1);
// 	}
// 	if (0 != strcmp (file_, "in")) {
// 	    ostrstream buff (Buffer, MAXLEN);
// 	    buff << "From aux(" << file_  << "): " << file_ << "_" << name;
// 	    buff << '\0';
// 	    putch (Buffer, "i", &value);
// 	}
//     } else {
// 	if (0 == getch(name,"s", value)) {
// 	    cerr << "missing parameter value: " << name << "\n";
// 	    exit (-1);
// 	}
// 	ostrstream buff (Buffer, MAXLEN);
// 	buff << "From par: " << name << '\0';
// 	putch (Buffer, "s", value);
//     }
// }

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

// 	$Id$	
