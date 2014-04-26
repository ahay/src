#include "parser.h"

/*#define VERBOSE */
#undef VERBOSE

/* max length of conversion string in set functions */
#define MAX_STR_LEN 1024

static char ps_sep_str[2] = { PS_SEP, '\0' };
/* static char ps_quo_str[2] = { PS_QUO, '\0' }; */

WORD * word_new() {
    WORD * w = (WORD *)usermalloc_(sizeof(WORD));
    if (!w) {
#ifdef VERBOSE
	fprintf(stderr,"Error: word_new - failed to allocate WORD\n");
#endif
	return w;
    }
    w->str=NULL;
    return w;
}

void word_delete(WORD ** w) {
    if ((*w)->str) userfree_((*w)->str);
    userfree_(*w);
    *w = NULL;
}

void word_reset(WORD * w) {
    if (w->str) userfree_(w->str);
    w->str = NULL;
}

/* post-construction initialization - use segment of a char array to
// initialize a word, which must have previously been constructed by
// word_new */
int word_assign(WORD * w, const char * str, int len) {
    int i;
    if (!w) {
#ifdef VERBOSE
	fprintf(stderr,"Error: word_assign - null word arg\n");
#endif
	return E_PARSE;
    }
    if (len<0) {
#ifdef VERBOSE
	fprintf(stderr,"Error: word_assign - len < 0\n");
#endif
	return E_PARSE;
    }
    if (w->str) userfree_(w->str);
    w->str = (char *)usermalloc_(sizeof(char)*(len+1));
    if (!(w->str)) {
#ifdef VERBOSE
	fprintf(stderr,"Error: word_assign - failed to allocate w->str\n");
#endif
	return E_ALLOC;
    }
    for (i=0;i<len;i++) (w->str)[i]=str[i];
    (w->str)[len]='\0';
    return 0;
}

int word_whitechar(char c) {
    if ((c == ' ')  ||
	(c == '\n') ||
	(c == '\t') ||
	(c == '\f')) return 1;
    return 0;
}

int word_copy(WORD * tgt, WORD src) {
    int err=0;
    if (!src.str) {
#ifdef VERBOSE 
	fprintf(stderr,"Error: word_copy - source WORD not defined, str=NULL\n");
#endif
	return E_PARSE;
    }
    if ((err=word_assign(tgt,src.str,strlen(src.str)+1))) {
#ifdef VERBOSE
	fprintf(stderr,"Error: word_copy from word_assign\n");
#endif
    }
    return err;
}

int word_read(WORD * w, char ** src) {
    /* finds a word in string str
    // allocates memory
    // on call, word should point to an allocated WORD 
    // on return, success indicated by w->str != NULL */

    int len=0;
    char * start;
    char * finger = *src;
    int qtoggle = 0;

#ifdef VERBOSE
    fprintf(stderr,"word_read: reset word\n");
#endif
    word_reset(w);

    /* if initial search turns up null char then no word can be found,
    // and word returns in its setnull initial state */
    if (*finger=='\0') return 0;

#ifdef VERBOSE
    fprintf(stderr,"word_read: on input src = \"%s\"\n",*src);
#endif

    /* find first non-whitespace char */
    while (word_whitechar(*finger)) finger++;

    /* have found a non-white, non-null char. Check to see if it is 
    //(a) the separator, or (b) the quote char */
    start=finger;
    /* quote charager toggles quoted string, which may have embedded 
    // whitespace */
    if (*finger == PS_QUO) qtoggle=!qtoggle;
#ifdef VERBOSE 
    fprintf(stderr,"word_read: quote toggled to %d\n",qtoggle);
#endif /* VERBOSE */
    if (*finger == PS_SEP) {
	word_assign(w,ps_sep_str,1);
	finger++;
#ifdef VERBOSE
	fprintf(stderr,"word_read, returning \"%s\"\n",w->str);
	fprintf(stderr,"word_read: on exit src = \"%s\"\n",*src);    
#endif
    }
    else {
	if (!qtoggle) {
	    /* keep scanning until you find either (a) whitespace or 
	    // (b) the separator char */ 
	    while (!(word_whitechar(*finger)   || 
		     (*finger == PS_SEP) ||
		     (*finger == PS_QUO) ||
		     (*finger == '\0'))) {
		finger++;
		len++;
	    }
	    word_assign(w,start,len);
#ifdef VERBOSE
	    fprintf(stderr,"word_read, returning \"%s\"\n",w->str);
	    fprintf(stderr,"word_read: on exit src = \"%s\"\n",*src);    
#endif
	}
	else {
#ifdef VERBOSE
	    fprintf(stderr,"word_read, quote branch, finger=%c\n",*finger);
#endif
	    finger++;
	    len=1;
	    while ((*finger != PS_QUO) &&
		   (*finger != '\0')) {
		finger++;
		len++;
	    }

	    /* check that we have found the end-of-quote, if so copy the 
	    // string */
	    if (*finger == PS_QUO) {
		/* then this char counts too */
		len++;
		word_assign(w,start,len);
		/* unset quote toggle */
		qtoggle=0;
		/* and move beyond final quote */
		finger++;
#ifdef VERBOSE
		fprintf(stderr,"word_read, returning %s\n",w->str);
		fprintf(stderr,"word_read: on exit src = \"%s\"\n",*src);    
#endif
	    }
	    else {
		/* otherwise, quote is not terminated */
#ifdef VERBOSE
		fprintf(stderr,"Error: word_read - quoted string not terminated with second quote\n");
#endif
		return E_PARSE;
	    }
	}
    }
    *src = finger;
    return 0;
}

KEYVAL * kv_new() {
    KEYVAL * kv = (KEYVAL *)usermalloc_(sizeof(KEYVAL));
    if (!kv) {
#ifdef VERBOSE
	fprintf(stderr,"Error: kv_new - failed to allocate new KEYVAL\n");
#endif
	return kv;
    }
    kv->key = word_new();
    if (!kv->key) {
#ifdef VERBOSE
	fprintf(stderr,"Error: kv_new - failed to allocate new KEYVAL key\n");
#endif
	userfree_(kv);
	return NULL;
    }
    kv->val = word_new();
    if (!kv->val) {
#ifdef VERBOSE
	fprintf(stderr,"Error: kv_new - failed to allocate new KEYVAL val\n");
#endif
	word_delete(&(kv->key));
	userfree_(kv);
	return NULL;
    }
    return kv;
}
  
void kv_delete(KEYVAL ** pair) {
    word_delete(&((*pair)->key));
    word_delete(&((*pair)->val));
    userfree_(*pair);
    *pair=NULL;
}

int kv_copy(KEYVAL * tgt, KEYVAL src) {
    int err=0;
    if ((err=word_copy(tgt->key, *(src.key)))) {
#ifdef VERBOSE
	fprintf(stderr,"Error: kv_copy from word_copy - key\n");
#endif
	return err;
    }
    if ((err=word_copy(tgt->val, *(src.val)))) {
#ifdef VERBOSE
	fprintf(stderr,"Error: kv_copy from word_copy - val\n");
#endif
	return err;
    }
    return err;
}

void kv_reset(KEYVAL * pair) {
    word_reset(pair->key);
    word_reset(pair->val);
}

int kv_check(KEYVAL src) {
    if ((src.key)->str && (src.val)->str) return 0;
    return 1;
}

int kv_read(KEYVAL * kv, char ** src) {
    int err=0;

#ifdef VERBOSE
    fprintf(stderr,"entering kv_read\n");
#endif
    WORD * lw = word_new();
    WORD * cw = word_new();
    WORD * nw = word_new();

#ifdef VERBOSE
    fprintf(stderr,"reset kv pair - erase key, value strings\n");
#endif
    kv_reset(kv);

#ifdef VERBOSE
    fprintf(stderr,"initialize - must find at least three word to play!\n");
#endif /* VERBOSE */
    if (word_read(lw,src) || !lw->str ||
	word_read(cw,src) || !cw->str ||
	word_read(nw,src) || !nw->str) {
	word_delete(&lw);
	word_delete(&cw);
	word_delete(&nw);
	return 0;
    }

    /* loop, look for cw = separator - if found, return implicit kv */
    while (kv_check(*kv)) {
#ifdef VERBOSE
	fprintf(stderr,"in kv_read: lw = %s\n",lw->str);
	fprintf(stderr,"in kv_read: cw = %s\n",cw->str);
	fprintf(stderr,"in kv_read: nw = %s\n",nw->str);
#endif
	/* have found a kv if (1) current word = separator, 
	// (2) neither previous nor next words = separator */
	if (!strcmp(cw->str,ps_sep_str) && 
	    strcmp(lw->str,ps_sep_str) &&
	    strcmp(nw->str,ps_sep_str)) {
	    if ((err=word_copy(kv->key,*lw))) {
#ifdef VERBOSE
		fprintf(stderr,"Error: kv_read from word_copy, key\n");
#endif
		return err;
	    }
	    if ((err=word_copy(kv->val,*nw))) {
#ifdef VERBOSE
		fprintf(stderr,"Error: kv_read from word_copy, val\n");
#endif
		return err;
	    }
	    word_delete(&lw);
	    word_delete(&cw);
	    word_delete(&nw);
	    return 0;
	}
	/* otherwise read another word */
	word_copy(lw,*cw);
	word_copy(cw,*nw);
	if (word_read(nw,src) || !nw->str) {
	    word_delete(&lw);
	    word_delete(&cw);
	    word_delete(&nw);
	    return 0;
	}
    }

#ifdef VERBOSE
    fprintf(stderr,"kv_read: deleting workspace, leave\n");
#endif
    word_delete(&lw);
    word_delete(&cw);
    word_delete(&nw);

    return 0;
}
      
void kv_fprint(KEYVAL kv, FILE * f) {
  if (!kv_check(kv)) fprintf(f,"%s%c%s\n",(kv.key)->str,PS_SEP,(kv.val)->str);
}

void kv_print(KEYVAL kv) { kv_fprint(kv,stdout); }

PSLINK * pslink_new() {
    PSLINK * par = (PSLINK *)usermalloc_(sizeof(PSLINK));
    if (!par) {
#ifdef VERBOSE
	fprintf(stderr,"Error: pslink_new - failed to allocate this\n");
#endif
	return NULL;
    }
    if (!(par->pair = kv_new())) {
#ifdef VERBOSE
	fprintf(stderr,"Error: pslink_new - failed to allocate KEYVAL data member\n");
#endif
	userfree_(par);
	return NULL;
    }
    par->prev = NULL;
    par->next = NULL;
    return par;
}

void pslink_delete(PSLINK ** p) {
    if (*p) {
	pslink_front(p);
	while ((*p)->next) {
	    *p=(*p)->next;
	    pslink_setnull(&((*p)->prev));
	}
	/* last one */
	pslink_setnull(p);
    }
}

void pslink_setnull(PSLINK ** p) {
    if (*p) {
	kv_delete(&((*p)->pair));
	userfree_(*p);
    }
    *p=NULL;
}

int pslink_front(PSLINK ** par) {
    if (!(*par)) {
#ifdef VERBOSE 
	fprintf(stderr,"Error: pslink_front - null pointer\n");
#endif
	return E_PARSE;
    }
    while ((*par)->prev) *par=(*par)->prev;
    return 0;
}

int pslink_back(PSLINK ** par) {
    if (!(*par)) {
#ifdef VERBOSE 
	fprintf(stderr,"Error: pslink_back - null pointer\n");
#endif
	return E_PARSE;
    }
    while ((*par)->next) *par=(*par)->next;
    return 0;
}

int pslink_read(PSLINK ** par, char ** str) {
    int stop = 0;
    if (!(*par)) {
#ifdef VERBOSE 
	fprintf(stderr,"Error: pslink_read - null pointer\n");
#endif
	return E_PARSE;
    }  
    /* reading always adds at the back, so make sure we are there 
    // first */
    if (pslink_back(par)) {
#ifdef VERBOSE
	fprintf(stderr,"Error: pslink_read from pslink_back\n");
#endif
	return E_PARSE;
    }     
    do {
	/* either read this link if it is blank, or advance to next link
	// and read there otherwise */
	if (kv_check(*((*par)->pair))) {
	    kv_read((*par)->pair,str);
	    /* check that read was successful */
	    stop = kv_check(*((*par)->pair));
	    /* if so move to next link */
	    if (!stop) {
		(*par)->next = pslink_new();
		(*par)->next->prev = *par;
		*par       = (*par)->next;
	    }
	    else {
		*par = (*par)->prev;
		pslink_setnull(&((*par)->next));  
	    }
	}
	else {
	    /* if current kv pair already initialized, do not overwrite, but 
	    // move to next link */
	    (*par)->next = pslink_new();
	    (*par)->next->prev = *par;
	    *par       = (*par)->next;
	    kv_read((*par)->pair,str);
	    stop = kv_check(*((*par)->pair));
	    /* if no read, we're at the end */
	    if (stop) {
		*par = (*par)->prev;
		pslink_setnull(&((*par)->next));  
	    }
	}	
    } while (!stop);
    return 0;
}

int pslink_findfirst(PSLINK * par, WORD skey, WORD * sval) {
    int err=0;
#ifdef VERBOSE
    fprintf(stderr,"pslink_findfirst\n");
#endif

    if ((err=pslink_front(&par))) {
	return err;
    }

    word_reset(sval);

#ifdef VERBOSE
    fprintf(stderr,"search key = %s\n",skey.str);
#endif
  
    do {
#ifdef VERBOSE
	fprintf(stderr,"compare key = %s\n",par->pair->key->str);
	fprintf(stderr,"search  key = %s\n",skey.str);
#endif
	if (par->pair->key->str && !strcmp(par->pair->key->str,skey.str)) {
#ifdef VERBOSE
	    fprintf(stderr,"found key match\n");
#endif
	    word_copy(sval,*(par->pair->val));
	    return 0;
	}
    } while ((par = par->next));

    return 0;
}

int pslink_findlast(PSLINK * par, WORD skey, WORD * sval) {
    int err=0;
    if ((err=pslink_back(&par))) {
	return err;
    }

    word_reset(sval);

#ifdef VERBOSE
    fprintf(stderr,"search key = %s\n",skey.str);
#endif
  
    do {
#ifdef VERBOSE
	fprintf(stderr,"key = %s\n",par->pair->key->str);
#endif
	if (par->pair->key->str && !strcmp(par->pair->key->str,skey.str)) {
	    word_copy(sval,*(par->pair->val));
	    return 0;
	}
    } while ((par = par->prev));

    return 0;
}

int pslink_setfirst(PSLINK ** par, WORD skey, WORD sval) {
    int err=0;
  
    if ((err=pslink_front(par))) {
#ifdef VERBOSE
	fprintf(stderr,"Error: pslink_setfirst - from pslink_front\n");
#endif
	return err;
    }

#ifdef VERBOSE
    fprintf(stderr,"pslink_setfirst: search key = %s\n",skey.str);
    fprintf(stderr,"pslink_setfirst: search val = %s\n",sval.str);
#endif
  
    while ((*par)->next) {
#ifdef VERBOSE
	fprintf(stderr,"key = %s ",(*par)->pair->key->str);
	if ((*par)->next) {
	    fprintf(stderr,"next != null\n");
	}
	else {
	    fprintf(stderr,"next == null\n");
	}
#endif
	if (!strcmp((*par)->pair->key->str,skey.str)) {
	    word_copy((*par)->pair->val,sval);
	    return 0;
	}
	*par = (*par)->next;
    } 
  
    /* back up one, since final "next" is null! */
    *par = (*par)->prev;

    /* if we get this far, no value for this key is present. 
    // add onto the list at the beginning */
    if ((err=pslink_front(par))) {
#ifdef VERBOSE
	fprintf(stderr,"Error: pslink_setfirst - null pointer at front\n");
#endif
	return err;
    }  
    if (!((*par)->prev = pslink_new())) {
#ifdef VERBOSE
	fprintf(stderr,"Error: pslink_setfirst - failed to allocate new pslink\n");
#endif
	return E_ALLOC;
    }
    (*par)->prev->next = (*par);
    if ((err=word_copy((*par)->prev->pair->key,skey))) {
#ifdef VERBOSE
	fprintf(stderr,"Error: pslink_setfirst - failed to copy key\n");
#endif
	return err;
    }
    if ((err=word_copy((*par)->prev->pair->val,sval))) {
#ifdef VERBOSE
	fprintf(stderr,"Error: pslink_setfirst - failed to copy val\n");
#endif
	return err;
    }
#ifdef VERBOSE
    fprintf(stderr,"pslink_setfirst: new pair = ");
    kv_fprint(*((*par)->prev->pair),stderr);
#endif
    return 0;
}

int pslink_setlast(PSLINK ** par, WORD skey, WORD sval) {
    int err=0;
  
    if ((err=pslink_back(par))) {
#ifdef VERBOSE
	fprintf(stderr,"Error: pslink_setlast - from pslink_back\n");
#endif
	return err;
    }

#ifdef VERBOSE
    fprintf(stderr,"pslink_setlast: search key = %s\n",skey.str);
    fprintf(stderr,"pslink_setlast: search val = %s\n",sval.str);
#endif
  
    while ((*par)->prev) {
	if (!strcmp((*par)->pair->key->str,skey.str)) {
	    word_copy((*par)->pair->val,sval);
	    return 0;
	}
	*par = (*par)->prev;
    } 
  
    /* back up one, since first "prev" is null! */
    if ((*par)->next) *par = (*par)->next;

    /* if we get this far, no value for this key is present. 
    // add onto the list at the end */
    if ((err=pslink_back(par))) {
#ifdef VERBOSE
	fprintf(stderr,"Error: pslink_setlast - null pointer at back\n");
#endif
	return err;
    }  
    if (!((*par)->next = pslink_new())) {
#ifdef VERBOSE
	fprintf(stderr,"Error: pslink_setlast - failed to allocate new pslink\n");
#endif
	return E_ALLOC;
    }
    (*par)->next->prev = (*par);
    if ((err=word_copy((*par)->next->pair->key,skey))) {
#ifdef VERBOSE
	fprintf(stderr,"Error: pslink_setlast - failed to copy key\n");
#endif
	return err;
    }
    if ((err=word_copy((*par)->next->pair->val,sval))) {
#ifdef VERBOSE
	fprintf(stderr,"Error: pslink_setlast - failed to copy val\n");
#endif
	return err;
    }
#ifdef VERBOSE
    fprintf(stderr,"pslink_setlast: new pair = ");
    kv_fprint(*((*par)->next->pair),stderr);
#endif
    return 0;
}

PARARRAY * ps_new() {
    PARARRAY * p = (PARARRAY *)usermalloc_(sizeof(PARARRAY));
    if (p) p->list = pslink_new();
    return p;
}

void ps_delete(PARARRAY ** p) {
    if (*p) {
	pslink_delete(&((*p)->list));
	userfree_(*p);
	*p=NULL;
    }
}

int ps_setnull(PARARRAY *par) {
    pslink_delete(&(par->list));
    par->list = pslink_new();
    return 0;
}

int ps_createfp(PARARRAY *parr, FILE * stream) {

    long size;                /* file size */
    char * str;               /* string */
    char * save;              /* start pointer */
    int n=0;                  /* counter */
    int err=0;                /* error code */
  
    /* get file size */
    if ( fseek(stream, 0L, SEEK_END) ) {
#ifdef VERBOSE
	fprintf(stderr,"Error: ps_createfp - seek failed\n");
#endif
	return E_FILE;
    }
    size = ftell(stream);
    if ( size == -1L ) {
#ifdef VERBOSE
	fprintf(stderr,"Error: ps_createfp - tell failed\n");
#endif
	return E_FILE;
    }
    else if ( size == 0L ) {
#ifdef VERBOSE
	fprintf(stderr,"NOTE: ps_createfp - zero length file\n");
#endif
	return 0;
    }

    rewind(stream);
  
    /* allocate memory, including space for trailing NUL */
    size++;
    str = (char*)usermalloc_(size);
    if ( !str ) { 
	return E_ALLOC;
    }
    /* reserve pointer to start of this segment */
    save=str;
  
    /* copy the file into the buffer */
    n = fread(str, 1L, size, stream);
  
    /* bail on read failure */
    if ( n != size-1 ) {
#ifdef VERBOSE
	fprintf(stderr,"Error: ps_createfp - read %d bytes expected %ld\n",n,size-1);
#endif

	userfree_(str);
	return E_FILE;
    }
  
    /* else tack on the trailing null */
    str[size-1] = '\0';
  
    /* read the string onto the end of the list */
    if ((err=pslink_read(&(parr->list),&str))) {
#ifdef VERBOSE
	fprintf(stderr,"Error: ps_createfp from pslink_read\n");
#endif
	userfree_(save);
	return err;
    }
  
    userfree_(save);
    return 0;
}

int ps_createfile_fproto(PARARRAY *parr, 
			 FILE ** fp, 
			 const char * proto,
			 const char * fname) {

  int err=0;                /* error code */
  
  /* sanity test - only on *fp = NULL! */
  if (*fp) {
#ifdef VERBOSE
    fprintf(stderr,"Error: ps_createfile_fproto - input stream non-null\n");
#endif
    return E_FILEOPEN;
  }
  /* open file */
  *fp = iwave_const_fopen(fname, "r", proto, stderr);
  if ( *fp == NULL ) {
#ifdef VERBOSE
    if (proto) 
      fprintf(stderr,"Error: ps_createfile_fproto - failed to open file=%s proto=%s\n",fname,proto);
    else
      fprintf(stderr,"Error: ps_createfile_fproto - failed to open file=%s proto=NULL\n",fname);
#endif
    return E_FILEOPEN;
  }
  
  err = ps_createfp(parr,*fp);
#ifdef VERBOSE
  if (err) 
    fprintf(stderr,"Error: ps_createfile_fproto from ps_createfp, filename=%s, err=%d\n",fname,err);
#endif
  
  return err;
}

int ps_createfile(PARARRAY *parr, const char *fname) {

  FILE * stream = NULL;
  char * proto = NULL;
  int err=0;                /* error code */

  err = ps_createfile_fproto(parr,&stream,proto,fname);
#ifdef VERBOSE
  if (err || !stream) 
    fprintf(stderr,"Error: ps_createfile from ps_createfile_fproto, filename=%s, err=%d\n",fname,err);
#endif

  iwave_fclose(stream);

  return err;
}

int ps_printall(PARARRAY parr, FILE *stream) {
    pslink_front(&(parr.list));
    while ((parr.list)->next) {
	kv_fprint(*((parr.list)->pair), stream);
	parr.list = parr.list->next;
    }
    kv_fprint(*((parr.list)->pair), stream);
    return 0;
}
/*@}*/

/* this construction is dangerous if let out, as the behaviour is 
// unpredictable if the void* pointer does not actually point to the
// indicated type - so should be used only locally, whence no mention
// in header file. */
int ps_get(PARARRAY * par, int last, const char * type, const char * key, void *p) {
    int err  = 0;   /* flag */
    /* workspace for stripping quotes */
    char * snq;     /* finger */
    int lnq;        /* word length */
  
    /* set up key, value buffers */
    WORD * skey = word_new();
    WORD * sval = word_new();
    if (!skey || !sval) err = E_ALLOC;
    if (!err) err=word_assign(skey,key,strlen(key)); 
  
    /* choose first or last */
    if (!err) {
	if (last) err=pslink_findlast(par->list,*skey,sval);
	else err=pslink_findfirst(par->list,*skey,sval);
    }
    if (err || !(sval->str)) err = iwave_max(err,1);
  
    if (!err && !strcmp(type,"cstring")) {
	char **s = (char **)p;
	/* modification 01.07.12: strip quotes, if any */
	if (((sval->str)[0]==PS_QUO) && ((sval->str)[strlen(sval->str)-1]==PS_QUO)) {
	    snq=&((sval->str)[1]);
	    lnq=strlen(sval->str)-2;
	}
	else {
	    snq=sval->str;
	    lnq=strlen(sval->str);
	}
	if (!(*s = (char *)usermalloc_(sizeof(char)*(strlen(sval->str)+1)))) err=E_ALLOC;    
	if (!err) strcpy(*s,snq);
	/* kill off last quote */
	(*s)[lnq]='\0';
    }
    else if (!err && !strcmp(type,"char")) {
	char * l = (char *)p;
	if (1!=sscanf(sval->str,"%c",l)) err=E_PARSECONVERT;
    }
    else if (!err && !strcmp(type,"long")) {
	long * l = (long *)p;
	if (1!=sscanf(sval->str,"%ld",l)) err=E_PARSECONVERT;
    }
    else if (!err && !strcmp(type,"ulong")) {
	unsigned long * l = (unsigned long *)p;
	if (1!=sscanf(sval->str,"%lu",l)) err=E_PARSECONVERT;
    }
    else if (!err && !strcmp(type,"int")) {
	int * l = (int *)p;
	if (1!=sscanf(sval->str,"%d",l)) err=E_PARSECONVERT;
    }
    else if (!err && !strcmp(type,"uint")) {
	unsigned int * l = (unsigned int *)p;
	if (1!=sscanf(sval->str,"%u",l)) err=E_PARSECONVERT;
    }
    else if (!err && !strcmp(type,"short")) {
	short * l = (short *)p;
	if (1!=sscanf(sval->str,"%hd",l)) err=E_PARSECONVERT;
    }
    else if (!err && !strcmp(type,"ushort")) {
	unsigned short * l = (unsigned short *)p;
	if (1!=sscanf(sval->str,"%hu",l)) err=E_PARSECONVERT;
    }    
    else if (!err && !strcmp(type,"float")) {
	float * l = (float *)p;
	if (1!=sscanf(sval->str,"%g",l)) err=E_PARSECONVERT;
    }
    else if (!err && !strcmp(type,"double")) {
	double * l = (double *)p;
	if (1!=sscanf(sval->str,"%lg",l)) err=E_PARSECONVERT;
    }
    else if (!err && !strcmp(type,"ireal")) {
#if DT_REAL == DT_DOUBLE
	double * l = (double *)p;
	if (1!=sscanf(sval->str,"%lg",l)) err=E_PARSECONVERT;
#else
	float * l = (float *)p;
	if (1!=sscanf(sval->str,"%g",l)) err=E_PARSECONVERT;
#endif
    }    
    else {
	err=E_PARSECONVERT;
    }
    word_delete(&skey);
    word_delete(&sval);  
    return err;
}

int ps_set(PARARRAY * par, int last, const char * type, const char * key, const void * val) {
    int err  = 0;     /* flag */
    char * sq;        /* finger, quoted string */
    char * l = (char *)usermalloc_(MAX_STR_LEN*sizeof(char));

    /* set up key, value buffers */
    WORD * skey = word_new();
    WORD * sval = word_new();
    if (!skey || !sval) err = E_ALLOC;
    if (!err) err=word_assign(skey,key,strlen(key)); 

    /* mod of 01.07.12: add quotes to any string, unless it alread begins and ends with
    // quotes */
    if (!err && !strcmp(type,"cstring")) {
	memset(l,'\0',MAX_STR_LEN);
	sq = (char *)val;
	if (strlen(sq) > MAX_STR_LEN) return E_PARSECONVERT;
	if ((sq[0]==PS_QUO) && (sq[strlen(sq)-1]==PS_QUO)) {
	    strcpy(l,sq);
	}
	else {
	    strcpy(&(l[1]),sq);
	    l[0]=PS_QUO;
	    l[strlen(sq)+1]=PS_QUO;
	}
	if (!err) err=word_assign(sval,l,strlen(l));
    }
    else {
	if (!err && !strcmp(type,"char")) {
	    memset(l,'\0',MAX_STR_LEN);
	    err=sprintf(l,"%c",*((char *)val));
	    if (0==err || MAX_STR_LEN==err) {
		err=E_PARSECONVERT;
	    }
	    else err=0;
	}
	else if (!err && !strcmp(type,"long")) {
	    memset(l,'\0',MAX_STR_LEN);
	    err=sprintf(l,"%ld",*((long *)val));
	    if (0==err || MAX_STR_LEN==err) {
		err=E_PARSECONVERT;
	    }
	    else err=0;
	}
	else if (!err && !strcmp(type,"ulong")) {
	    memset(l,'\0',MAX_STR_LEN);
	    err=sprintf(l,"%lu",*((unsigned long *)val));
	    if (0==err || MAX_STR_LEN==err) {
		err=E_PARSECONVERT;
	    }
	    else err=0;
	}
	else if (!err && !strcmp(type,"int")) {
	    /*      fprintf(stderr,"in ps_set - int\n"); */
	    memset(l,'\0',MAX_STR_LEN);
	    err=sprintf(l,"%d",*((int *)val));
	    /*      fprintf(stderr,"output = %s err=%d\n",l,err); */
	    if (0==err || MAX_STR_LEN==err) {
		err=E_PARSECONVERT;
	    }
	    else err=0;
	}
	else if (!err && !strcmp(type,"uint")) {
	    memset(l,'\0',MAX_STR_LEN);
	    err=sprintf(l,"%u",*((unsigned int *)val));
	    if (0==err || MAX_STR_LEN==err) {
		err=E_PARSECONVERT;
	    }
	    else err=0;
	}
	else if (!err && !strcmp(type,"short")) {
	    memset(l,'\0',MAX_STR_LEN);
	    err=sprintf(l,"%hd",*((short *)val));
	    if (0==err || MAX_STR_LEN==err) {
		err=E_PARSECONVERT;
	    }
	    else err=0;
	}
	else if (!err && !strcmp(type,"ushort")) {
	    memset(l,'\0',MAX_STR_LEN);
	    err=sprintf(l,"%hu",*((unsigned short *)val));
	    if (0==err || MAX_STR_LEN==err) {
		err=E_PARSECONVERT;
	    }
	    else err=0;
	}    
	else if (!err && !strcmp(type,"float")) {
	    memset(l,'\0',MAX_STR_LEN);
	    err=sprintf(l,"%g",*((float *)val));
	    if (0==err || MAX_STR_LEN==err) {
		err=E_PARSECONVERT;
	    }
	    else err=0;
	}
	else if (!err && !strcmp(type,"double")) {
	    memset(l,'\0',MAX_STR_LEN);
	    err=sprintf(l,"%g",*((double *)val));
	    if (0==err || MAX_STR_LEN==err) {
		err=E_PARSECONVERT;
	    }
	    else err=0;
	}
	else if (!err && !strcmp(type,"ireal")) {
	    memset(l,'\0',MAX_STR_LEN);
#if DT_REAL == DT_DOUBLE
	    err=sprintf(l,"%g",*((double *)val));
#else
	    err=sprintf(l,"%g",*((float *)val));
#endif
	    if (0==err || MAX_STR_LEN==err) {
		err=E_PARSECONVERT;
	    }
	    else err=0;
	}    
	else {
	    err=E_PARSECONVERT;
	}
	/*    fprintf(stderr,"assigning value = %s err=%d\n",l,err); */

	if (!err) err=word_assign(sval,l,strlen(l));

    }
    /* choose first or last
    //  fprintf(stderr,"ps_set: %s =  %s\n",skey->str,sval->str); */
    if (!err) {
	if (last) err=pslink_setlast(&(par->list),*skey,*sval);
	else err=pslink_setfirst(&(par->list),*skey,*sval);
    }

    userfree_(l);
    word_delete(&skey);
    word_delete(&sval);

    return err;
}

/* this is the main public post-construction initialization */
int ps_createargs(PARARRAY *par, int argc, char **argv) {

    int n;                   /* length of arg string */
    int i;                   /* counter, error code */
    int err = 0;
    char *buffer;            /* arg buffer */
    char *save;              /* reserve start of buffer */
    char *tmp;               /* another reserve pointer */
    KEYVAL * kv;             /* workspace for decoding kv pairs */
  
    /* bail if no args */
    if ( argc <= 0 ) return 0;
  
    /* build string out of args */
    n = 0;
    for ( i = 0; i < argc; ++i ) n += strlen(argv[i]) + 2;
    /* trailing null */
    n++;

    /* allocate arg buffer */
    buffer = (char*)usermalloc_(n * sizeof(char));
    if ( buffer == NULL ) return E_ALLOC;
    memset(buffer,'\0',n);
    save=buffer;
   
    if (!(kv=kv_new())) {
	userfree_(buffer);
	return E_ALLOC;
    }
    /* copy the args into the buffer
       in course of copy, look for "par=" and 
       if found initialize PARARRAY 
    */
    for ( i = 0; i < argc; ++i ) {
	tmp=argv[i];
	kv_reset(kv);
	if (!kv_read(kv, &tmp) &&
	    !kv_check(*kv)           &&
	    !strcmp(kv->key->str,"par")) {
	    /* have found par file */
	    err=ps_createfile(par,kv->val->str);
	    if (err) {
#ifdef VERBOSE
		fprintf(stderr,"Error: ps_createargs from ps_createfile\n");
#endif
		return err;
	    }
	}
	else {
#ifdef VERBOSE
	    fprintf(stderr,"ps_createargs: adding %s\n",argv[i]);
#endif
	    strcat(buffer,argv[i]);
	    strcat(buffer," ");
	}
    }
#ifdef VERBOSE 
    fprintf(stderr,"ps_createargs: arg buffer =\n");
    fprintf(stderr,"%s\n",buffer);
#endif

    /* now add on the arg pairs - these are appended, and do 
    // not replace already initialized keys. So to give primacy to 
    // the arg list over a par file, use find-last - other way around,
    // use find-first. */
    if (!err) err=pslink_read(&(par->list),&buffer);
    if (err) {
#ifdef VERBOSE
	fprintf(stderr,"ERROR: ps_createargs from pslink_read\n");
#endif
	return err;
    }

#ifdef VERBOSE 
    ps_printall(*par,stderr);
#endif

    kv_delete(&kv);
 
    userfree_(save);
    return err;
}

/* public functions */

/* find first */

int ps_ffcstring(PARARRAY par, const char *key, char **p) {
    return ps_get(&par,0,"cstring",key,p);
}
int ps_ffchar(PARARRAY par, const char *key, char *p) {
    return ps_get(&par,0,"char",key,p);
}
int ps_ffshort(PARARRAY par, const char *key, short *p) {
    return ps_get(&par,0,"short",key,p);
}
int ps_ffint(PARARRAY par, const char *key, int *p) {
    return ps_get(&par,0,"int",key,p);
}

int ps_fflong(PARARRAY par, const char *key, long *p) {
    return ps_get(&par,0,"long",key,p);
}
int ps_ffushort(PARARRAY par, const char *key, unsigned short *p) {
    return ps_get(&par,0,"ushort",key,p);
}
int ps_ffuint(PARARRAY par, const char *key, unsigned int *p) {
    return ps_get(&par,0,"uint",key,p);
}
int ps_ffulong(PARARRAY par, const char *key, unsigned long *p) {
    return ps_get(&par,0,"ulong",key,p);
}
int ps_fffloat(PARARRAY par, const char *key, float *p) {
    return ps_get(&par,0,"float",key,p);
}
int ps_ffdouble(PARARRAY par, const char *key, double *p) {
    return ps_get(&par,0,"double",key,p);
}
int ps_ffreal(PARARRAY par, const char *key, ireal *p) {
    return ps_get(&par,0,"ireal",key,p);
}

/* find-last */

int ps_flcstring(PARARRAY par, const char *key, char **p) {
  int err = ps_get(&par,1,"cstring",key,p);
  return err;
}
int ps_flchar(PARARRAY par, const char *key, char *p) {
    return ps_get(&par,1,"char",key,p);
}
int ps_flshort(PARARRAY par, const char *key, short *p) {
    return ps_get(&par,1,"short",key,p);
}
int ps_flint(PARARRAY par, const char *key, int *p) {
    return ps_get(&par,1,"int",key,p);
}
int ps_fllong(PARARRAY par, const char *key, long *p) {
    return ps_get(&par,1,"long",key,p);
}
int ps_flushort(PARARRAY par, const char *key, unsigned short *p) {
    return ps_get(&par,1,"ushort",key,p);
}
int ps_fluint(PARARRAY par, const char *key, unsigned int *p) {
    return ps_get(&par,1,"uint",key,p);
}
int ps_flulong(PARARRAY par, const char *key, unsigned long *p) {
    return ps_get(&par,1,"ulong",key,p);
}
int ps_flfloat(PARARRAY par, const char *key, float *p) {
    return ps_get(&par,1,"float",key,p);
}
int ps_fldouble(PARARRAY par, const char *key, double *p) {
    return ps_get(&par,1,"double",key,p);
}
int ps_flreal(PARARRAY par, const char *key, ireal *p) {
    return ps_get(&par,1,"ireal",key,p);
}

/* assign first */

int ps_sfcstring(PARARRAY par, const char *key, const char *p) {
    return ps_set(&par,0,"cstring",key,p);
}
int ps_sfchar(PARARRAY par, const char *key, char p) {
    return ps_set(&par,0,"char",key,&p);
}
int ps_sfshort(PARARRAY par, const char *key, short p) {
    return ps_set(&par,0,"short",key,&p);
}
int ps_sfint(PARARRAY par, const char *key, int p) {
    return ps_set(&par,0,"int",key,&p);
}
int ps_sflong(PARARRAY par, const char *key, long p) {
    return ps_set(&par,0,"long",key,&p);
}
int ps_sfushort(PARARRAY par, const char *key, unsigned short p) {
    return ps_set(&par,0,"ushort",key,&p);
}
int ps_sfuint(PARARRAY par, const char *key, unsigned int p) {
    return ps_set(&par,0,"uint",key,&p);
}
int ps_sfulong(PARARRAY par, const char *key, unsigned long p) {
    return ps_set(&par,0,"ulong",key,&p);
}
int ps_sffloat(PARARRAY par, const char *key, float p) {
    return ps_set(&par,0,"float",key,&p);
}
int ps_sfdouble(PARARRAY par, const char *key, double p) {
    return ps_set(&par,0,"double",key,&p);
}
int ps_sfreal(PARARRAY par, const char *key, ireal p) {
    return ps_set(&par,0,"ireal",key,&p);
}

/* assign last */

int ps_slcstring(PARARRAY par, const char *key, const char *p) {
    return ps_set(&par,1,"cstring",key,p);
}
int ps_slchar(PARARRAY par, const char *key, char p) {
    return ps_set(&par,1,"char",key,&p);
}
int ps_slshort(PARARRAY par, const char *key, short p) {
    return ps_set(&par,1,"short",key,&p);
}
int ps_slint(PARARRAY par, const char *key, int p) {
    return ps_set(&par,1,"int",key,&p);
}
int ps_sllong(PARARRAY par, const char *key, long p) {
    return ps_set(&par,1,"long",key,&p);
}
int ps_slushort(PARARRAY par, const char *key, unsigned short p) {
    return ps_set(&par,1,"ushort",key,&p);
}
int ps_sluint(PARARRAY par, const char *key, unsigned int p) {
    return ps_set(&par,1,"uint",key,&p);
}
int ps_slulong(PARARRAY par, const char *key, unsigned long p) {
    return ps_set(&par,1,"ulong",key,&p);
}
int ps_slfloat(PARARRAY par, const char *key, float p) {
    return ps_set(&par,1,"float",key,&p);
}
int ps_sldouble(PARARRAY par, const char *key, double p) {
    return ps_set(&par,1,"double",key,&p);
}
int ps_slreal(PARARRAY par, const char *key, ireal p) {
    return ps_set(&par,1,"ireal",key,&p);
}

int ps_copy(PARARRAY ** tgt, PARARRAY src) {

    PSLINK * slst;      /* source list object */
    PSLINK * tlst;      /* target list object */

    if (*tgt) ps_delete(tgt);
  
    slst = src.list;
    tlst = NULL;

    /* start at the beginning... */
    if (pslink_front(&slst)) {
	return E_OTHER;
    }
  
    while (slst) {
	/* first pair: tgt has been deleted, assign new PARARRAY */
	if (!(*tgt)) *tgt=ps_new();
	else {
	    /* target PARARRAY has been initialized, create next link
	    // move up */
	    (*tgt)->list->next = pslink_new();
	    (*tgt)->list->next->prev = (*tgt)->list;
	    (*tgt)->list = (*tgt)->list->next;
	}
	tlst = (*tgt)->list;
	word_copy(tlst->pair->key,*(slst->pair->key));
	word_copy(tlst->pair->val,*(slst->pair->val));
	slst = slst->next;
    }

    return 0;

}

namespace RVL {

  bool parse(PARARRAY const & par, std::string name, string & val) {
    char * p;
    int res=ps_flcstring(par,name.c_str(),&p);
    if (!res) {
      val=p;
      userfree_(p);
      return true;
    }
    return false;
  }

  bool parse(PARARRAY const & par, std::string name, bool & val) {
    int ival; 
    if (parse(par, name, ival)) {
      if (ival==0) val=false;
      else val=true;
      return true;
    }
    return false;
  }

  bool parse(PARARRAY const & par, std::string name, char & val) {
    return !ps_flchar(par,name.c_str(),&val);
  }

  bool parse(PARARRAY const & par, std::string name, short & val) {
    return !ps_flshort(par,name.c_str(),&val);
  }

  bool parse(PARARRAY const & par, std::string name, int & val) {
    return !ps_flint(par,name.c_str(),&val);
  }

  bool parse(PARARRAY const & par, std::string name, long & val) {
    return !ps_fllong(par,name.c_str(),&val);
  }

  bool parse(PARARRAY const & par, std::string name, unsigned short & val) {
    return !ps_flushort(par,name.c_str(),&val);
  }

  bool parse(PARARRAY const & par, std::string name, unsigned int & val) {
    return !ps_fluint(par,name.c_str(),&val);
  }

  bool parse(PARARRAY const & par, std::string name, unsigned long & val) {
    return !ps_flulong(par,name.c_str(),&val);
  }

  bool parse(PARARRAY const & par, std::string name, float & val) {
    return !ps_flfloat(par,name.c_str(),&val);
  }

  bool parse(PARARRAY const & par, std::string name, double & val) {
    return !ps_fldouble(par,name.c_str(),&val);
  }

}


