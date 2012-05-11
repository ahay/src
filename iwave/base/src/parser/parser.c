/* 
parser.c
Igor Terentyev.
*/
/*============================================================================*/

#include "utils.h"
#include "parser.h"
#include "convert.h"
/* added 24.04.10 WWS */
#include "iwave_fopen.h"

/*----------------------------------------------------------------------------*/
/*
Symbol type.
*/
#define GRAPHICAL 0
#define SEPARATOR 1
#define DELIMITER 2
/*----------------------------------------------------------------------------*/

SIZEDSTRING ps_name2z(const char *name)
{
    SIZEDSTRING namez;
    namez.n = strlen(name);
    namez.s = (char*)name;
    return namez;
}
/*----------------------------------------------------------------------------*/

int ps_setnull(PARARRAY *parr)
{
    memset((void*)parr, 0, sizeof(PARARRAY));    
    return 0;
}
/*----------------------------------------------------------------------------*/

int ps_destroy(PARARRAY *parr)
{
  if (parr->buffer) free(parr->buffer);
  if (parr->ns)  free(parr->ns);
  ps_setnull(parr);
  
  return 0;
}
/*----------------------------------------------------------------------------*/

int ps_concat(PARARRAY *parr, PARARRAY parr2)
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
		     const char * value, int mv) {

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
// memory bust discovered by Muhong Zhou 03.12 - temporary
// fix creates memory leak
//    free(str.s);
    
    return err;
}
/*----------------------------------------------------------------------------*/

int ps_createargs(PARARRAY *parr, int argc, char **argv)
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
{
    long num = name.n * sizeof(char);
    int ip, count = 0;

    for ( ip = 0; ip < parr.npar; ++ip )
        if ( (parr.ns[ip].n == name.n) && (memcmp(name.s, parr.ns[ip].s, num) == 0) ) ++count;

    return count;
}
/*----------------------------------------------------------------------------*/

int ps_countname(PARARRAY parr, const char *name)
{
    return ps_countnamez(parr, ps_name2z(name));
}
/*----------------------------------------------------------------------------*/

int ps_getvalz(PARARRAY parr, SIZEDSTRING name, int occ, SIZEDSTRING *val)
{
    int ip = ps_nameindexz(parr, name, occ);
    if ( ip < 0 ) return E_PARSENONAME;
    *val = parr.vs[ip];
    return 0;
}
/*----------------------------------------------------------------------------*/

int ps_getval(PARARRAY parr, const char *name, int occ, SIZEDSTRING *val)
{
    return ps_getvalz(parr, ps_name2z(name), occ, val);
}
/*----------------------------------------------------------------------------*/

int ps_nameindexz(PARARRAY parr, SIZEDSTRING name, int occ)
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
{
    return ps_nameindexz(parr, ps_name2z(name), occ);
}
/*----------------------------------------------------------------------------*/

int ps_val2type(SIZEDSTRING val, int type, void *p, SIZEDSTRING *valend)
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
{
    return ps_counttypez(parr, ps_name2z(name), occ, type);
}
/*----------------------------------------------------------------------------*/

int ps_gettypez(PARARRAY parr, SIZEDSTRING name, int occ, int type, int n, void *p)
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
{
    return ps_gettypez(parr, ps_name2z(name), occ, type, n, p);
}
/*----------------------------------------------------------------------------*/
int ps_ffchar(PARARRAY parr, const char *name, char *p)
{
    return ps_gettype(parr, name, 0, DT_CHAR, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffshort(PARARRAY parr, const char *name, short *p)
{
    return ps_gettype(parr, name, 0, DT_SHORT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffint(PARARRAY parr, const char *name, int *p)
{
    return ps_gettype(parr, name, 0, DT_INT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_fflong(PARARRAY parr, const char *name, long *p)
{
    return ps_gettype(parr, name, 0, DT_LONG, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffushort(PARARRAY parr, const char *name, unsigned short *p)
{
    return ps_gettype(parr, name, 0, DT_USHORT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffuint(PARARRAY parr, const char *name, unsigned int *p)
{
    return ps_gettype(parr, name, 0, DT_UINT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffulong(PARARRAY parr, const char *name, unsigned long *p)
{
    return ps_gettype(parr, name, 0, DT_ULONG, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/
int ps_fffloat(PARARRAY parr, const char *name, float *p)
{
    return ps_gettype(parr, name, 0, DT_FLOAT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffdouble(PARARRAY parr, const char *name, double *p)
{
    return ps_gettype(parr, name, 0, DT_DOUBLE, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffreal(PARARRAY parr, const char *name, ireal *p)
{
    return ps_gettype(parr, name, 0, DT_REAL, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_ffcstring(PARARRAY parr, const char *name, char **p)
{
    return ps_gettype(parr, name, 0, DT_CSTRING, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/
int ps_flchar(PARARRAY parr, const char *name, char *p)
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_CHAR, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_flshort(PARARRAY parr, const char *name, short *p)
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_SHORT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_flint(PARARRAY parr, const char *name, int *p)
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_INT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_fllong(PARARRAY parr, const char *name, long *p)
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_LONG, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_flushort(PARARRAY parr, const char *name, unsigned short *p)
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_USHORT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_fluint(PARARRAY parr, const char *name, unsigned int *p)
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_UINT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_flulong(PARARRAY parr, const char *name, unsigned long *p)
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_ULONG, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/
int ps_flfloat(PARARRAY parr, const char *name, float *p)
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_FLOAT, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_fldouble(PARARRAY parr, const char *name, double *p)
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_DOUBLE, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_flreal(PARARRAY parr, const char *name, ireal *p)
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_REAL, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_flcstring(PARARRAY parr, const char *name, char **p)
{
    return ps_gettype(parr, name, ps_countname(parr, name) - 1, DT_CSTRING, 0, (void*)p);
}
/*----------------------------------------------------------------------------*/

int ps_printall(PARARRAY parr, FILE *stream)
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

