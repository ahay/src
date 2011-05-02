/* Command-line or parameter file parser. */
/*************************************************************************

Copyright Rice University, 2008.
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

/* 
parser.c
Igor Terentyev.
*/
/*============================================================================*/
/** @page parser IWAVE Parameter Parser 
    Author: Igor Terentyev
    <p>
    file: \ref parser.h
    <p>
    Decomposes file/string into pairs of strings: [key][value], and provides
    functions to convert value to prescribed type.
    Separator control symbol between key and value is '='
    (can be redefined in \ref utils.h file).
    Delimiters between key=value pairs are any non-graphical symbols.
    Keys CANNOT be empty: if key not found between separators,
    separator is ignored (treated as delimiter). Values can be empty.
    Quote control symbol is '"' (can be redefined in utils.h file).
    Everything between quotes is treated as graphical symbol.
    '"' symbol (non-control) is given by [""]. 
    Strings do not have to be null-terminated. String size is provided.
    For just in case strings produced by parser are by one symbol longer
    then their ireal sizes and have null-terminator in the extra space.
    <p>
    EXAMPLES: on left of "==>", input string (in square brackets). On right, [key] [value] produced by parser.
    <ul>
    <li>
    [a=b]
    ==> [a][b].
    </li>

    <li>
    [a = b]
    ==> [a][b].
    </li>

    <li>
    [a =b = c = d]
    ==> [a][b],[c][d]; second '=' ignored since no key to the left from it.
    </li>

    <li>
    ["  key w/ spaces" = "value w/ spaces and ""quotes"""]
    ==> [  key w/ spaces][value w/ spaces and "quotes"]
    </li>

    <li>
    [eq_sign="="]
    ==> [eq_sign][=].
    </li>

    <li>
    [a = 1 comment 1 b = 2 comment 2 c = 3]
    ==> [a][1],[b][2],[c][3].
    </li>

    <li>
    [ a

    <p>

    =

    <p>

    1]
    ==> [a][1].
    </li>

    <li>
    [a= = (empty parameter) b = 1 (non-empty parameter)]
    ==> [a][],[b][1].
    </li>
    </ul>
    <p>
    The package provides several types of functions which manipulate parameter
    arrays and access parameter values:
    <ol>
    <li>\ref create</li>
    <li>\ref print</li>
    <li>\ref ffaccess</li>
    <li>\ref flaccess</li>
    <li>\ref utility</li>
    </ol>
*/

#include "utils.h"
#include "parser.h"

#include "convert.h"
/*^*/

/* added 24.04.10 WWS */
#include "iwave_fopen.h"

#ifndef _sf_parser_h

/*
Parameter structure.
--------------------
int npar        :  number of parameters in the array.
SIZEDSTRING *ns :  keys array.
SIZEDSTRING *vs :  values array.
char *buffer    :  allocated buffer.
*/
typedef struct
{
  int npar;
  SIZEDSTRING *ns;
  SIZEDSTRING *vs;
  char *buffer;
	int buffersize;
} PARARRAY;
/*^*/

/*----------------------------------------------------------------------------*/
/*
Symbol type.
*/
#define GRAPHICAL 0
#define SEPARATOR 1
#define DELIMITER 2
/*----------------------------------------------------------------------------*/
/*^*/

#endif

SIZEDSTRING ps_name2z(const char *name)
/*< Convert zero terminated string to sized string (NO ALLOCATION). >*/
{
    SIZEDSTRING namez;
    namez.n = strlen(name);
    namez.s = (char*)name;
    return namez;
}
/*----------------------------------------------------------------------------*/

int ps_setnull(PARARRAY *parr)
/*< ** 
Set parameter array (all fields) to zeros - null initialization
@param[out] parr (PARARRAY *) - parameter array to initialize
@return (int) 0 if successful, else nonzero error code.
>*/
{
    memset((void*)parr, 0, sizeof(PARARRAY));    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ps_destroy(PARARRAY *parr)
/*< ** 
Destroy parameter array (STORAGE DEALLOCATION). 
@param[out] parr (PARARRAY *) - param array destroyed on successful return
@return (int) 0 if successful, else nonzero error code.
>*/
{
  if (parr->buffer) free(parr->buffer);
  if (parr->ns)  free(parr->ns);
  ps_setnull(parr);
  
  return 0;
}
/*----------------------------------------------------------------------------*/

int ps_concat(PARARRAY *parr, PARARRAY parr2)
/*< ** 
Concatenate parameter arrays (STORAGE REALLOCATION).
"parr = parr + parr2".
@param[out] parr (PARARRAY *) - first summand on call, sum on return
@param[in] parr2 (PARARRAY)  - second summand (const)
@return (int) 0 if successful, else nonzero error code.
>*/
{
    char *buffer, *ptr;
    SIZEDSTRING *ns, *vs;
    int npar, buffersize, i;
  
    if ( parr2.npar <= 0 ) return 0; /* nothing to attach */
    
    /* new sizes */
    buffersize = parr->buffersize + parr2.buffersize;
    npar = parr->npar + parr2.npar;
    
    /* new buffer */
    buffer = (char*)malloc(buffersize * sizeof(char));
    if ( buffer == NULL ) return E_ALLOC;
    
    /* new ns */
    ns = (SIZEDSTRING*)malloc(2 * npar * sizeof(SIZEDSTRING));
    if ( ns == NULL )
    {
        free(buffer);
        return E_ALLOC;
    }
    vs = ns + npar;
    
    /* copy data */
    ptr = buffer;
    for ( i = 0; i < parr->npar; ++i )
    {
        ns[i].s = ptr;
        ns[i].n = parr->ns[i].n;
        memcpy(ptr, parr->ns[i].s, (ns[i].n + 1) * sizeof(char));
        ptr += ns[i].n + 1;
    }
    for ( i = 0; i < parr2.npar; ++i )
    {
        ns[parr->npar + i].s = ptr;
        ns[parr->npar + i].n = parr2.ns[i].n;
        memcpy(ptr, parr2.ns[i].s, (ns[parr->npar + i].n + 1) * sizeof(char));
        ptr += ns[parr->npar + i].n + 1;
    }
    for ( i = 0; i < parr->npar; ++i )
    {
        vs[i].s = ptr;
        vs[i].n = parr->vs[i].n;
        memcpy(ptr, parr->vs[i].s, (vs[i].n + 1) * sizeof(char));
        ptr += vs[i].n + 1;
    }
    for ( i = 0; i < parr2.npar; ++i )
    {
        vs[parr->npar + i].s = ptr;
        vs[parr->npar + i].n = parr2.vs[i].n;
        memcpy(ptr, parr2.vs[i].s, (vs[parr->npar + i].n + 1) * sizeof(char));
        ptr += vs[parr->npar + i].n + 1;
    }

    /* delete old data */
    ps_destroy(parr);

    /* store new data */
    parr->npar = npar;
    parr->buffer = buffer;
    parr->ns = ns;
    parr->vs = vs;
    parr->buffersize = buffersize;
    
    return 0;
}

/*----------------------------------------------------------------------------*/
int ps_addpairstring(PARARRAY *parr, 
		     const char * key, int mk,
		     const char * value, int mv) 
/*< ** 
Add "key=value" line from string input (STORAGE REALLOCATION).
@param[out] parr (PARARRAY *) - contains additional key=val pair on return
@param[in] key (const char *)  - null-terminated char array = key
@param[in] mk (int) - allocated size of key, should be >= strlen(key)+1;
@param[in] value (const char *) - null-terminated char array = value
@param[in] mv (int) - allocated size of value, should be >= strlen(value)+1;
@return (int) 0 if successful, else nonzero error code.
>*/
{

  int nk;           /* key string len */
  int nv;           /* val string len */
  int err=0;        /* error flag */
  
  SIZEDSTRING sz;   /* workspace */
  PARARRAY addpar;  /* workspace */

  /* sanity-check inputs - this is probably sufficient */
  nk=strlen(key);
  nv=strlen(value);
  if ((nk<1)||(nv<1)||(mk<nk+1)||(mv<nv+1)) {
    fprintf(stderr,"ERROR: ps_addpairstring - bad input, mk=%d nk=%d mv=%d nv=%d\n",mk,nk,mv,nv);
    return E_BADINPUT;
  }

  /* build up sized string */
  sz.n=nk+nv+4;
  sz.s=(char*)malloc(sz.n*sizeof(char));
  strcpy(sz.s,key);
  strcat(sz.s," = ");
  strcat(sz.s,value);

  /* construct parray workspace */
  ps_setnull(&addpar);
  if ((err=ps_createstrz(&addpar,sz))) {
    fprintf(stderr,"ERROR: ps_addpairstring from ps_createstrz, err=%d\n",err);
    return err;
  }

  /* tack the PARARRAY data member of the FO onto the newly-created
     PARARRAY */
  if ((err=ps_concat(&addpar,*parr))) {
    fprintf(stderr,"ERROR: ps_addpairstring from ps_concat, err=%d\n",err);
    return err;
  }

  /* copy the newly created PARARRAY over the data member */
  if (ps_copy(addpar,parr)) return err;

  /* clean up */
  ps_destroy(&addpar);
  free(sz.s);

  return err;
}

/*----------------------------------------------------------------------------*/

int ps_createfile(PARARRAY *parr, const char *filename)
/*< ** 
Creates parameter array (STORAGE ALLOCATION) from file
@param[out] parr (PARARRAY *) - param array created on successful return
@param[in] filename (char *) - name of parfile containing key=value info
@return (int) 0 if successful, else nonzero error code.
>*/
{
    FILE *stream;             /* file stream */
    long size;                /* file size */
    SIZEDSTRING str;          /* string */
    int err;                  /* error code */
    
    /* clear parameter structure */
    ps_setnull(parr);

    /* open file */
    stream = iwave_const_fopen(filename, "r", NULL, stderr);
    if ( stream == NULL ) return E_FILEOPEN;

    /* get file size */
    if ( fseek(stream, 0L, SEEK_END) )
    {
        iwave_fclose(stream);
        return E_FILE;
    }
    size = ftell(stream);
    if ( size == -1L )
    {
        iwave_fclose(stream);
        return E_FILE;
    }
    else if ( size == 0L )
    {
        iwave_fclose(stream);
        return 0;
    }
    rewind(stream);

    /* allocate memory */
    str.s = (char*)malloc(size);
    if ( str.s == NULL ) 
    {
        iwave_fclose(stream);
        return E_ALLOC;
    }
    
    /* copy the file into the buffer */
    str.n = fread(str.s, 1L, size, stream);
    iwave_fclose(stream);

    if ( str.n != size )
    {
        free(str.s);
        return E_FILE;
    }
    
    err = ps_createstrz(parr, str);
    free(str.s);
    
    return err;
}
/*----------------------------------------------------------------------------*/

int ps_createargs(PARARRAY *parr, int argc, char **argv)
/*< ** 
Creates parameter array (STORAGE ALLOCATION) from command-line argument list
@param[out] parr (PARARRAY *) - param array created on successful return
@param[in] argc (int) - arg count
@param[in] argv (char**) - arg array
@return (int) 0 if successful, else nonzero error code.
>*/
{
    long *sizes;             /* string sizes including null-terminator */
    SIZEDSTRING str;         /* string */
    int i, err;              /* error code */
    char *buffer;
    
    /* clear parameter structure */
    ps_setnull(parr);

    /* calculate size */
    if ( argc <= 0 ) return 0;
        
    sizes = (long*)malloc(argc * sizeof(long));
    if ( sizes == NULL ) return E_ALLOC;
        
    str.n = 0L;
    for ( i = 0; i < argc; ++i )
    {
        sizes[i] = strlen(argv[i]) + 1L;
        str.n += sizes[i];
    }
    if ( str.n == 0L ) 
    {
        free(sizes);
        return 0;
    }
        
    /* allocate memory */
    str.s = (char*)malloc(str.n * sizeof(char));
    if ( str.s == NULL ) 
    {
        free(sizes);
        return E_ALLOC;
    }
    
    /* copy the strings into the buffer */
    buffer = str.s;
    for ( i = 0; i < argc; ++i ) 
    {
        strncpy(buffer, argv[i], sizes[i]);
        buffer += sizes[i];
        buffer[-1] = ' ';    /* substitute null with space symbol */
    }

    free(sizes);
    err = ps_createstrz(parr, str);
    free(str.s);
    
    return err;
}
/*----------------------------------------------------------------------------*/

int ps_createstrz(PARARRAY *parr, SIZEDSTRING str)
/*< ** 
Creates parameter array (STORAGE ALLOCATION) from string.
@param[out] parr (PARARRAY *) - param array created on successful return
@param[in] str (SIDEDSTRING) - sized string struct from which to draw 
key=value pairs
@return (int) 0 if successful, else nonzero error code.
>*/
{
    char *s = str.s;
    long n = str.n;
    long pr, ps, pe, pl;
    int npar, nsep, fquo, fnew, fside, ftype;
    SIZEDSTRING *ns, *vs;
    char *buffer, c;
    
    /*    fprintf(stderr,"ps_createstrz: str.n = %d str.s = %s\n",str.n,str.s);*/
    /* clear parameter structure, self-check */
    ps_setnull(parr);
    if ( PS_SEP == PS_QUO ) return E_INTERNAL; 
    
    /* count separators */
    nsep = fquo = 0;
    for ( pr = 0L; pr < n; ++pr )
    {
        c = s[pr];
        if ( c == PS_QUO )
        {
            if ( (pr + 1L < n) && (s[pr + 1L] == PS_QUO) ) ++pr; /* double */
            else fquo = 1 - fquo;                                /* single */
            continue;
        }
        if ( (c == PS_SEP) && (!fquo) ) ++nsep;
    }
    if ( fquo ) return E_PARSE; /* QUO not closed */
    if ( nsep == 0L ) return 0;

    /* allocate memory */
    buffer = (char*)malloc((n + 1L) * sizeof(char));
    ns = (SIZEDSTRING*)malloc(2L * nsep * sizeof(SIZEDSTRING));
    if ( (buffer == NULL) || (ns == NULL) )
    {
        free(buffer);
        free(ns);
        return E_ALLOC;
    }
    vs = ns + nsep;

    /* start in left mode */    
    pl = ps = pe = 0L;
    fnew = 1;
    npar = fquo = fside = 0;
    
    /* cycle through read buffer */
    for ( pr = 0L; pr <= n; ++pr )
    {
        c = (pr < n) ? s[pr] : PS_SEP;
        /* process QUO symbol */
        if ( c == PS_QUO )
        {
            if ( (pr + 1L < n) && (s[pr + 1L] == PS_QUO) ) ++pr; /* double */
            else { fquo = 1 - fquo; continue; }                  /* single */
        }
        /* determine symbol type */
        if      ( fquo       ) ftype = GRAPHICAL;
        else if ( c == PS_SEP) ftype = SEPARATOR;
        else if ( isgraph(c) ) ftype = GRAPHICAL;
        else                   ftype = DELIMITER;
        /* check type */
        if ( ftype == SEPARATOR ) --nsep;
        else if ( (ftype == DELIMITER) && (pe == ps) ) continue;
        /* process right mode */
        if ( fside )
        {
            if ( ftype == GRAPHICAL )
            {
                buffer[pe++] = c;
            }
            else /* DELIMITER/SEPARATOR after non-empty value*/
            {
                vs[npar - 1].s = buffer + ps;
                vs[npar - 1].n = pe - ps;
                buffer[pe] = 0;
                pl = ps = (++pe);
                fnew = 1;
                fside = 0;
                if ( nsep <= 0 ) break;
            }
        }
        /* process left mode */
        else
        {
            if ( ftype == SEPARATOR )
            {
                if ( pe > ps )
                {
                    ns[npar].s = buffer + ps;
                    ns[npar].n = pe - ps;
                    buffer[pe] = 0;
                    ps = (++pe);
                    fnew = 1;
                    fside = 1;
                    ++npar;
                }
                else if ( nsep <= 0 ) break;
            }
            else if ( ftype == GRAPHICAL )
            {
                if ( fnew )
                {
                    pe = ps;
                    fnew = 0;
                }
                buffer[pe++] = c;
            }
            else /* DELIMITER after non-empty name */
            {
                fnew = 1;
            }
        }
    }
    
    /* free and return if no parameters */
    if ( npar == 0 )
    {
        free(buffer);
        free(ns);
        return 0;
    }
    
    /* reallocate */
    memcpy(ns + npar, vs, npar * sizeof(SIZEDSTRING));
    parr->buffer = (char*)realloc(buffer, pl * sizeof(char));
    parr->ns = (SIZEDSTRING*)realloc(ns, 2L * npar * sizeof(SIZEDSTRING));
    if ( (parr->buffer == NULL) || (parr->ns == NULL) )
    {
        free(buffer);
        free(ns);
        return E_ALLOC;
    }
    parr->vs = parr->ns + npar;
    parr->npar = npar;
    parr->buffersize = pl;
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ps_countnamez(PARARRAY parr, SIZEDSTRING name)
/*< Return number of parameter key occurences. >*/
{
    long num = name.n * sizeof(char);
    int ip, count = 0;

    for ( ip = 0; ip < parr.npar; ++ip )
        if ( (parr.ns[ip].n == name.n) && (memcmp(name.s, parr.ns[ip].s, num) == 0) ) ++count;

    return count;
}
/*----------------------------------------------------------------------------*/

int ps_countname(PARARRAY parr, const char *name)
/*< Return number of parameter key occurences. >*/
{
    return ps_countnamez(parr, ps_name2z(name));
}
/*----------------------------------------------------------------------------*/

int ps_getvalz(PARARRAY parr, SIZEDSTRING name, int occ, SIZEDSTRING *val)
/*< 
Returns sized string value of the parameter from key.

int occ :  key occurence number.
>*/
{
    int ip = ps_nameindexz(parr, name, occ);
    if ( ip < 0 ) return E_PARSENONAME;
    *val = parr.vs[ip];
    return 0;
}
/*----------------------------------------------------------------------------*/

int ps_getval(PARARRAY parr, const char *name, int occ, SIZEDSTRING *val)
/*< 
Returns sized string value of the parameter from key.

int occ :  key occurence number.
>*/
{
    return ps_getvalz(parr, ps_name2z(name), occ, val);
}
/*----------------------------------------------------------------------------*/

int ps_nameindexz(PARARRAY parr, SIZEDSTRING name, int occ)
/*<
Return index of parameter occurence. 
Returns -1 if the occurence not found.

int occ :  key occurence number.
>*/
{
    long num = name.n * sizeof(char);
    int ip, count = -1;

    for ( ip = 0; ip < parr.npar; ++ip )
        if ( (parr.ns[ip].n == name.n) && (memcmp(name.s, parr.ns[ip].s, num) == 0) )
            if ( (++count) == occ ) return ip;
            
    return -1;
}
/*----------------------------------------------------------------------------*/

int ps_nameindex(PARARRAY parr, const char *name, int occ)
/*<
Return index of parameter occurence. 
Returns -1 if the occurence not found.

int occ :  key occurence number.
>*/
{
    return ps_nameindexz(parr, ps_name2z(name), occ);
}
/*----------------------------------------------------------------------------*/

int ps_val2type(SIZEDSTRING val, int type, void *p, SIZEDSTRING *valend)
/*<
Seeks next value and converts it.
MAY ALLOCATE (if type is cstring).
>*/
{
    char *str;
    long len;
    
    switch ( type )
    {
        case DT_CHAR:       return strz2char(val, (char*)p, valend);
        case DT_SHORT:      return strz2short(val, (short*)p, valend);
        case DT_INT:        return strz2int(val, (int*)p, valend);
        case DT_LONG:       return strz2long(val, (long*)p, valend);
        case DT_USHORT:     return strz2ushort(val, (unsigned short*)p, valend);
        case DT_UINT:       return strz2uint(val, (unsigned int*)p, valend);
        case DT_ULONG:      return strz2ulong(val, (unsigned long*)p, valend);
        case DT_FLOAT:      return strz2float(val, (float*)p, valend);
        case DT_DOUBLE:     return strz2double(val, (double*)p, valend);
        case DT_CSTRING:
            if ( valend != NULL ) *valend = val;
            if ( val.n == 0L ) return E_PARSENOVALUE; /* empty SIZEDSTRING */
            len = strlen(val.s) + 1L;  /* string lenght including terminator */
            if ( p != NULL )           /* allocate and copy */
            {
                str = (char*)calloc(len, sizeof(char));
                if ( str == NULL ) return E_ALLOC;
                strncpy(str, val.s, len);
                *((char**)p) = str;
            }
            if ( valend != NULL )      /* shift */
            {
                if ( len > val.n )
                {
                    val.s += val.n;
                    val.n = 0L;
                }
                else
                {
                    val.s += len;
                    val.n -= len;
                }
            }
            return 0;
    }
    return E_BADINPUT;
}
/*----------------------------------------------------------------------------*/

int ps_counttypez(PARARRAY parr, SIZEDSTRING name, int occ, int type)
/*<
Return number of parameter type value occurences in the parameter value. 

int occ :  key occurence number.
>*/
{
    SIZEDSTRING val;
    int err, count = 0;
    
    err = ps_getvalz(parr, name, occ, &val);
    if ( err ) return -1;
    
    while ( 1 )
    {
        if ( ps_val2type(val, type, NULL, &val) > 0 ) break;
        ++count;
    }
    return count;
}
/*----------------------------------------------------------------------------*/

int ps_counttype(PARARRAY parr, const char *name, int occ, int type)
/*<
Return number of parameter type value occurences in the parameter value. 

int occ :  key occurence number.
>*/
{
    return ps_counttypez(parr, ps_name2z(name), occ, type);
}
/*----------------------------------------------------------------------------*/

int ps_gettypez(PARARRAY parr, SIZEDSTRING name, int occ, int type, int n, void *p)
/*<
Extract type value.

int occ :  key occurence number.
int n   :  type occurence number in the parameter value.
>*/
{
    SIZEDSTRING val;
    int err, count;
    
    err = ps_getvalz(parr, name, occ, &val);
    if ( err ) return err;
    
    for ( count = 0; count < n; ++count) 
        if ( ps_val2type(val, type, NULL, &val) > 0 ) return E_PARSENOVALUE;
            
    return ps_val2type(val, type, p, NULL);
}
/*----------------------------------------------------------------------------*/

int ps_gettype(PARARRAY parr, const char *name, int occ, int type, int n, void *p)
/*<
Extract type value.

int occ :  key occurence number.
int n   :  type occurence number in the parameter value.
>*/
{
    return ps_gettypez(parr, ps_name2z(name), occ, type, n, p);
}
/*----------------------------------------------------------------------------*/
int ps_ffchar(PARARRAY parr, const char *name, char *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, 0, DT_CHAR, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffshort(PARARRAY parr, const char *name, short *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, 0, DT_SHORT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffint(PARARRAY parr, const char *name, int *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, 0, DT_INT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_fflong(PARARRAY parr, const char *name, long *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, 0, DT_LONG, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffushort(PARARRAY parr, const char *name, unsigned short *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, 0, DT_USHORT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffuint(PARARRAY parr, const char *name, unsigned int *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, 0, DT_UINT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffulong(PARARRAY parr, const char *name, unsigned long *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, 0, DT_ULONG, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/
int ps_fffloat(PARARRAY parr, const char *name, float *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, 0, DT_FLOAT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffdouble(PARARRAY parr, const char *name, double *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, 0, DT_DOUBLE, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffreal(PARARRAY parr, const char *name, ireal *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, 0, DT_REAL, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffcstring(PARARRAY parr, const char *name, char **p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, 0, DT_CSTRING, 0, (void*)p);
}

/*----------------------------------------------------------------------------*/
int ps_flchar(PARARRAY parr, const char *name, char *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_CHAR, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_flshort(PARARRAY parr, const char *name, short *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_SHORT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_flint(PARARRAY parr, const char *name, int *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_INT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_fllong(PARARRAY parr, const char *name, long *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_LONG, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_flushort(PARARRAY parr, const char *name, unsigned short *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_USHORT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_fluint(PARARRAY parr, const char *name, unsigned int *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_UINT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_flulong(PARARRAY parr, const char *name, unsigned long *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_ULONG, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_flfloat(PARARRAY parr, const char *name, float *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_FLOAT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_fldouble(PARARRAY parr, const char *name, double *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_DOUBLE, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_flreal(PARARRAY parr, const char *name, ireal *p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_REAL, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_flcstring(PARARRAY parr, const char *name, char **p)
/*< ** \defgroup ffaccess Parameter access - first occurence 

Get first occurence of data in the first occurence of the key. 

@param[in] parr (PARARRAY)  -  parameter array.
@param[in] key (char *)    -  key (null-terminated string).
@param[out] p (type *)       -  pointer to value of indicated type.
@return (int) 0 if successful, else nonzero error code.

Note: the string case allocates memory, which must be managed by the
calling unit; in other cases the last arg points to memory already 
allocated by the calling unit.
>*/
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_CSTRING, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_printall(PARARRAY parr, FILE *stream)
/*< **
Write contents of PARARRAY to stream
@param[in] parr (PARARRAY) - input param array
@param[out] stream (FILE *) - output stream
@return (int) 0 if successful, else nonzero error code.
>*/
{
    int i, n;

    n = parr.npar;
    /*    fprintf(stream, "Number of parameters: %d\n",n); */
    /*    for ( i = 0; i < n; ++i ) fprintf(stream, "%d [%s] [%s]\n", i, parr.ns[i].s, parr.vs[i].s); */
    for ( i = 0; i < n; ++i ) fprintf(stream, "%s %c %s\n", parr.ns[i].s, PS_SEP, parr.vs[i].s);
    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ps_copy(PARARRAY rhs, PARARRAY *lhs)
/*< **
Full copy function: lhs = rhs (STORAGE (RE)ALLOCATED)
@param lhs (PARARRAY *) - target parameter array
@param rhs (PARARRAY)   - source parameter array (const)
@return (int) 0 if successful, else nonzero error code.
>*/
{
  int i;
  
  ps_destroy(lhs);
  if ( rhs.buffersize == 0 ) return 0;
  
  /* Allocate and copy buffer. */
  lhs->buffer = (char*)malloc(rhs.buffersize * sizeof(char));
  if ( lhs->buffer == NULL ) return E_ALLOC;
  
  memcpy(lhs->buffer, rhs.buffer, rhs.buffersize * sizeof(char));
  
  /* Allocate and copy pointers. */
  lhs->ns = (SIZEDSTRING*)malloc(2L * rhs.npar * sizeof(SIZEDSTRING));
  if ( lhs->ns == NULL )
  {
    free(lhs->buffer);
    return E_ALLOC;
  }
  
  for ( i = 0; i < 2 * rhs.npar; ++i )
  {
    lhs->ns[i].n = rhs.ns[i].n;
    lhs->ns[i].s = lhs->buffer + (rhs.ns[i].s - rhs.buffer);
  }

  /* Copy other fields. */
  lhs->buffersize = rhs.buffersize;
  lhs->npar = rhs.npar;
  lhs->vs = lhs->ns + rhs.npar;

  return 0;
}
/*----------------------------------------------------------------------------*/

