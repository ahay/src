// File: rsf.cc
// --------------
// RSF C++ interface
/////////////////////////////////////////////////

#include <vector>

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

// iRSF Constructor
// ----------------
// Standard input
iRSF::iRSF()
{
    file_ = sf_input("in");
}

// oRSF Constructor
// ----------------
// Standard output
oRSF::oRSF() 
{
    file_ = sf_output("out");
}

// oRSF Constructor
// ----------------
// Copies the file tag, writes out the grid
oRSF::oRSF(const char* file /* = "out" */)
{
    assert (file); // no NULLs (yet)
    file_ = sf_output(file);
}

// iRSF Destructor
// ----------------
iRSF::~iRSF()
{
    if (0 != file_) sf_fileclose(file_);
}

// oRSF Destructor
// ----------------
oRSF::~oRSF()
{
    sf_fileclose(file_);
}

// Reading data
// ------------
const iRSF&
iRSF::operator>> (std::vallarray<float> &array) const
{
    assert (0 != file_);
    sf_read(&(array[0]),sizeof(float),array.size(),file_);
    return *this;
}

// Writing data
// ------------
const oRSF&
oRSF::operator<< (std::valarray<float> &array) const
{
    sf_write(&(array[0]),sizeof(float),array.size(),file_);
    return *this;
}

// Reading parameters
/////////////////////
void 
iRSF::get (char* name, int& value, int default) const
{
    if (0 != file_)  {
	if (!sf_histint (file_,name,&value)) value = default;
    } else {
	if (!sf_getint (name, &value)) value = default;
    }
}

// void 
// iRSF::get (char* name, int& value) const
// {
//     if (0 != file_)  {
// 	if (0 == auxpar(name,"i", &value, file_)) {
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
// 	if (0 == getch(name,"i", &value)) { 
// 	    cerr << "missing parameter value: " << name << "\n";
// 	    exit (-1);
// 	}    
// 	ostrstream buff (Buffer, MAXLEN);
// 	buff << "From par: " << name << '\0';
// 	putch (Buffer, "i", &value);
//     }
// }
 
void 
iRSF::get (char* name, float& value, float default) const
{
    if (0 != file_)  {
	if (0 == auxpar (name,"r", &value, file_)) value = default;
	if (0 != strcmp (file_, "in")) {
	    ostrstream buff (Buffer, MAXLEN);
	    buff << "From aux(" << file_  << "): " << file_ << "_" << name;
	    buff << '\0';
	    putch (Buffer, "r", &value);
	}
    } else {
	if (0 == getch (name, "r", &value)) value = default;
	ostrstream buff (Buffer, MAXLEN);
	buff << "From par: " << name << '\0';
	putch (Buffer, "r", &value);
    }
}

void 
iRSF::get (char* name, float& value) const
{
    if (0 != file_)  {
	if (0 == auxpar(name,"r", &value, file_)) {
	    cerr << "missing history value: " << name;
	    cerr << "from file: " << file_ << "\n";
	    exit (-1);
	}
	if (0 != strcmp (file_, "in")) {
	    ostrstream buff (Buffer, MAXLEN);
	    buff << "From aux(" << file_  << "): " << file_ << "_" << name;
	    buff << '\0';
	    putch (Buffer, "r", &value);
	}
    } else {
	if (0 == getch(name,"r", &value)) {
	    cerr << "missing parameter value: " << name << "\n";
	    exit (-1);
	}
	ostrstream buff (Buffer, MAXLEN);
	buff << "From par: " << name << '\0';
	putch (Buffer, "r", &value);
    }
}

void 
iRSF::get (char* name, char* value, const char* default) const
{
    if (0 != file_)  {
	if (0 == auxpar (name,"s", value, file_)) strcpy (value, default);
	if (0 != strcmp (file_, "in")) {
	    ostrstream buff (Buffer, MAXLEN);
	    buff << "From aux(" << file_  << "): " << file_ << "_" << name;
	    buff << '\0';
	    putch (Buffer, "s", value);
	}
    } else {
	if (0 == getch (name, "s", value)) strcpy (value, default);
	ostrstream buff (Buffer, MAXLEN);
	buff << "From par: " << name<< '\0';
	putch (Buffer, "s", value);
    }
}

void 
iRSF::get (char* name, char* value) const
{
    if (0 != file_) {
	if (0 == auxpar(name,"s", value, file_)) {
	    cerr << "missing history value: " << name;
	    cerr << "from file: " << file_ << "\n";
	    exit (-1);
	}
	if (0 != strcmp (file_, "in")) {
	    ostrstream buff (Buffer, MAXLEN);
	    buff << "From aux(" << file_  << "): " << file_ << "_" << name;
	    buff << '\0';
	    putch (Buffer, "i", &value);
	}
    } else {
	if (0 == getch(name,"s", value)) {
	    cerr << "missing parameter value: " << name << "\n";
	    exit (-1);
	}
	ostrstream buff (Buffer, MAXLEN);
	buff << "From par: " << name << '\0';
	putch (Buffer, "s", value);
    }
}

// Reading parameter arrays
///////////////////////////
void 
iRSF::get (char* name, int size, int*   value, const int*   default) const
{
    if (0 != file_) {
	if (0 != auxpar(name, "i",value, file_)) return;
    } else if (0 != getch(name, "i",value)) return;
    int buflen = strlen(name) + size%10 + 2;
    char* buffer = new char[buflen];
    for (int i = 0; i < size; i++) {
	ostrstream buff (buffer, buflen);
	buff << name << i+1 << '\0';
	get (buffer, value[i], default[i]);
    }
    delete [] buffer;
}

void 
iRSF::get (char* name, int size, int*   value) const
{
    if (0 != file_) {
	if (0 != auxpar(name, "i",value, file_)) return;
    } else if (0 != getch(name, "i",value)) return;
    int buflen = strlen(name) + size%10 + 2;
    char* buffer = new char[buflen];
    for (int i = 0; i < size; i++) {
	ostrstream buff (buffer, buflen);
	buff << name << i+1 << '\0';
	get (buffer, value[i]);
    }
    delete [] buffer;
}

void 
iRSF::get (char* name, int size, float* value, const float* default) const
{
    if (0 != file_) {
	if (0 != auxpar(name, "r",value, file_)) return;
    } else if (0 != getch(name, "r",value)) return;
    int buflen = strlen(name) + size%10 + 2;
    char* buffer = new char[buflen];
    for (int i = 0; i < size; i++) {
	ostrstream buff (buffer, MAXLEN);
	buff << name << i+1 << '\0';
	get (buffer, value[i], default[i]);
    }
    delete [] buffer;
}

void 
iRSF::get (char* name, int size, float* value) const
{
    if (0 != file_) {
	if (0 != auxpar(name, "r",value, file_)) return;
    } else if (0 != getch(name, "r",value)) return;
    int buflen = strlen(name) + size%10 + 2;
    char* buffer = new char[buflen];
    for (int i = 0; i < size; i++) {
	ostrstream buff (buffer, MAXLEN);
	buff << name << i+1 << '\0';
	get (buffer, value[i]);
    }
    delete [] buffer;
}

// Writing parameters
/////////////////////

void 
oRSF::put (char* name, int value) const
{
    auxputch(name,"i",&value,file_);
}
 
void 
oRSF::put (char* name, float value) const
{
    auxputch(name,"r",&value,file_);
}

void 
oRSF::put (char* name, const char* value) const
{
    auxputch(name,"s",value,file_);
}

// Writing parameter arrays
///////////////////////////
void 
oRSF::put (char* name, int size, const int*   value) const
{ 
    int i;
    if (size <= MAXPAR && (0 == strcmp (file_, "out"))) { 
	// write "n=1,2,3,34,.."
	ostrstream buff (Buffer, MAXLEN);
	buff << "     " << name << '=';
	for (i = 0; i < size-1; i++) 
	    buff << value[i] << ',';
	buff << value[size-1] << '\0';
	putlin (Buffer);
    } else {             // write "n1=1 n2=2 ..."
	for (i = 0; i < size; i++) {
	    ostrstream buff (Buffer, MAXLEN);
	    buff << name << i+1 << '\0';
	    put (Buffer, value[i]);
	}
    }
}

void 
oRSF::put (char* name, int size, const float* value) const
{
    int i;
    if (size <= MAXPAR&& (0 == strcmp (file_, "out"))) { 
	// write "n=1,2,3,34,.."
	ostrstream buff (Buffer, MAXLEN);
	buff << "     " << name << '=';
	for (i = 0; i < size-1; i++) 
	    buff << value[i] << ',';
	buff << value[size-1] << '\0';
	putlin (Buffer);
    } else {             // write "n1=1 n2=2 ..."
	for (i = 0; i < size; i++) {
	    ostrstream buff (Buffer, MAXLEN);
	    buff << name << i+1 << '\0';
	    put (Buffer, value[i]);
	}
    }
}
