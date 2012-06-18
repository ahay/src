#include "keyval.h"

//#define VERBOSE
#undef VERBOSE

static char ps_sep_str[2] = { PS_SEP, '\0' };
//static char ps_quo_str[2] = { PS_QUO, '\0' };

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

// post-construction initialization - use segment of a char array to
// initialize a word, which must have previously been constructed by
// word_new
int word_assign(WORD * w, char * str, int len) {
    int i;
    if (!w) {
#ifdef VERBOSE
	fprintf(stderr,"Error: word_assign - null word arg\n");
#endif
	return E_PARSE;
    }
    if (len<1) {
#ifdef VERBOSE
	fprintf(stderr,"Error: word_assign - len < 1\n");
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
    if ((c == ' ') ||
	(c == '\n')) return 1;
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
    // finds a word in string str
    // allocates memory
    // on call, word should point to an allocated WORD 
    // on return, success indicated by w->str != NULL

    int len=0;
    char * start;
    char * finger = *src;
    int qtoggle = 0;

#ifdef VERBOSE
    fprintf(stderr,"word_read: reset word\n");
#endif
    word_reset(w);

    // if initial search turns up null char then no word can be found,
    // and word returns in its setnull initial state
    if (*finger=='\0') return 0;

#ifdef VERBOSE
    fprintf(stderr,"word_read: on input src = \"%s\"\n",*src);
#endif

    // find first non-whitespace char 
    while (word_whitechar(*finger)) finger++;

    // have found a non-white, non-null char. Check to see if it is 
    //(a) the separator, or (b) the quote char
    start=finger;
    // quote charager toggles quoted string, which may have embedded 
    // whitespace
    if (*finger == PS_QUO) qtoggle=!qtoggle;
#ifdef VERBOSE 
    fprintf(stderr,"word_read: quote toggled to %d\n",qtoggle);
#endif /* VERBOSE */
    if (*finger == PS_SEP) {
	word_assign(w,ps_sep_str,1);
	finger++;
#ifdef VERBOSE
	fprintf(stderr,"word_read, returning %s\n",w->str);
	fprintf(stderr,"word_read: on exit src = \"%s\"\n",*src);    
#endif
    }
    else {
	if (!qtoggle) {
	    // keep scanning until you find either (a) whitespace or 
	    // (b) the separator char 
	    while (!(word_whitechar(*finger)   || 
		     (*finger == PS_SEP) ||
		     (*finger == PS_QUO) ||
		     (*finger == '\0'))) {
		finger++;
		len++;
	    }
	    word_assign(w,start,len);
#ifdef VERBOSE
	    fprintf(stderr,"word_read, returning %s\n",w->str);
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

	    // check that we have found the end-of-quote, if so copy the 
	    // string
	    if (*finger == PS_QUO) {
		// then this char counts too
		len++;
		word_assign(w,start,len);
		// unset quote toggle
		qtoggle=0;
		// and move beyond final quote
		finger++;
#ifdef VERBOSE
		fprintf(stderr,"word_read, returning %s\n",w->str);
		fprintf(stderr,"word_read: on exit src = \"%s\"\n",*src);    
#endif
	    }
	    else {
		// otherwise, quote is not terminated
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
	word_read(nw,src) || !nw->str) return 0;

    // loop, look for cw = separator - if found, return implicit kv
    while (kv_check(*kv)) {
#ifdef VERBOSE
	fprintf(stderr,"in kv_read: lw = %s\n",lw->str);
	fprintf(stderr,"in kv_read: cw = %s\n",cw->str);
	fprintf(stderr,"in kv_read: nw = %s\n",nw->str);
#endif
	// have found a kv if (1) current word = separator, 
	// (2) neither previous nor next words = separator
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
	    return 0;
	}
	// otherwise read another word
	word_copy(lw,*cw);
	word_copy(cw,*nw);
	if (word_read(nw,src) || !nw->str) return 0;
    }
    return 0;
}
      
void kv_print(KEYVAL kv) {
    if (!kv_check(kv)) printf("%s %c %s\n",(kv.key)->str,PS_SEP,(kv.val)->str);
}

void kv_fprint(KEYVAL kv, FILE * f) {
    if (!kv_check(kv)) fprintf(f,"%s %c %s\n",(kv.key)->str,PS_SEP,(kv.val)->str);
}

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
	kv_delete(&((*p)->pair));
	userfree_(*p);
	*p = NULL;
    }
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
    do {
	kv_read((*par)->pair,str);
	stop = kv_check(*((*par)->pair));
	if (!stop) {
	    (*par)->next = pslink_new();
	    (*par)->next->prev = *par;
	    *par       = (*par)->next;
	}
	else {
	    *par = (*par)->prev;
	    pslink_delete(&((*par)->next));  
	}
    } while (!stop);
    return 0;
}

int pslink_findfirst(PSLINK * par, WORD skey, WORD * sval) {
    int err=0;
    if ((err=pslink_front(&par))) {
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
	if (!strcmp(par->pair->key->str,skey.str)) {
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
	if (!strcmp(par->pair->key->str,skey.str)) {
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
  
    // back up one, since final "next" is null!
    *par = (*par)->prev;

    // if we get this far, no value for this key is present. 
    // add onto the list at the beginning
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

int pslink_setlast(PSLINK * par, WORD skey, WORD sval) {
    int err=0;
    if ((err=pslink_back(&par))) {
	return err;
    }

#ifdef VERBOSE
    fprintf(stderr,"search key = %s\n",skey.str);
    fprintf(stderr,"search val = %s\n",sval.str);
#endif
  
    do {
#ifdef VERBOSE
	fprintf(stderr,"key = %s\n",par->pair->key->str);
#endif
	if (!strcmp(par->pair->key->str,skey.str)) {
	    word_copy(par->pair->val,sval);
	    return 0;
	}
    } while ((par = par->prev));

    return 0;
}



