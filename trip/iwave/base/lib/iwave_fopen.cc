#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif

#include "iwave_fopen.h"

#define NLEN 16
#define BLEN 1000
#define CONST_FP

static char buf[BLEN];

static struct filestat {
  struct filestat *nextfpr;         /* linked list pointer */
  FILE * fp;                        /* FILE */
  char * nm;                        /* NAME */
  char * pr;                        /* PROTO */
  char * md;                        /* MODE */
  int istmp;                        /* temp flag */
  int inuse;                        /* availability flag */
} * filestatlist = (struct filestat *)NULL;

static struct filestat *fpr, **oldfpr;

FILE * iwave_fopen(char ** name, 
		   const char * mode, 
		   const char * proto,
		   FILE * stream) {
  
  FILE * retfp = (FILE *)NULL; /* return value */
  FILE * fph   = (FILE *)NULL; /* proto file pointer */
  int fd;                      /* file descr for return from mkstemp */
  int nr;                      /* #words read/written in copy of proto */
  int chunk;                   /* size of buffer read in copy of proto */
  off_t plen;                  /* size of prototype file */   
#ifdef CONST_FP
  off_t sav;                   /* saved offset prior to reopen */
#endif

  /* sanity check: if proto file given, then previously opened for read */
  if (proto) {
    /*    fprintf(stderr,"begin proto - search for %s\n",proto);*/
    oldfpr=&filestatlist;
    for (fpr=filestatlist; fpr != ((struct filestat *)NULL); 
	 fpr = fpr->nextfpr) {
      if (!(strcmp(proto,fpr->nm)) &&
	  (fpr->istmp==0) && 
	  /*
	  (fpr->md[0]=='r')*/
	  (strcmp("w",fpr->md))) break;
      oldfpr=&fpr->nextfpr;
    }

    if (fpr==((struct filestat *)NULL)) {
      /* didn't find it */
      fprintf(stream,"Error: iwave_fopen\n");
      fprintf(stream,"requested proto file %s not in database\n",proto);
      iwave_fprintall(stream);
      fflush(stream);
      return retfp;
    }      

    /* successful search for proto file; assign pointer */
    fph=fpr->fp;
    /*    fprintf(stderr,"found proto file %s\n",fpr->nm);*/
  }

  /* CASE I - request for temp file

  if name==NULL, then a temp file is being requested. Co-conditions:
  - a prototype file exists - this must be a previously-opened archival 
  file, and proto!=NULL with file permission permitting read
  - w(+) file permissiion being requested
  if these conditions hold, then check for existing filename with same
  proto and temp status - if you find one not in use, return its FILE*,
  else open a new file, copy the prototype file to it, and return the 
  pointer.
  */

  if (*name==NULL) {

    /*    fprintf(stderr,"begin temp\n");*/
    
    /* sanity check: access mode = w+ */
    if (strcmp(mode,"w+")) {
      fprintf(stream,"Error: iwave_fopen\n");
      fprintf(stream,"the only legal permission for a temp file is w+\n");
      return retfp;
    }

    /* sanity check: prototype file legal */
    if (proto==NULL || fph==NULL) {
      fprintf(stream,"Error: iwave_fopen\n");
      fprintf(stream,"cannot request temp file without file prototype\n");
      return retfp;
    }
      
    /*    fprintf(stream,"  proto file %s\n",proto);*/
    /*    iwave_fprintall(stream);*/
    /* having found legal prototype file, search database for unused
       temp file modeled on same proto */

    oldfpr=&filestatlist;
    for (fpr=filestatlist; fpr != ((struct filestat *)NULL);fpr = fpr->nextfpr) {

      /* break if (1) prototypes match, (2) filestat has tmp status, and
	 (3) associated file not in use */
      /*      fprintf(stream,"-- check filename=%s\n",fpr->nm);*/
      if ((fpr->pr) && 
	  !(strcmp(proto,fpr->pr)) &&
	  (fpr->istmp==1) &&
	  (fpr->inuse==0)) break;
      /*      fprintf(stream,"-- next\n");*/
      oldfpr=&(fpr->nextfpr);
    }

    /*    fprintf(stream,"  search finished\n");*/
    /* CASE IA: no existing temp file available, open new one */

    if (fpr==((struct filestat *)NULL)) {

      /* first identify DATAPATH, if it exists*/
      char * dpath = getenv("DATAPATH");
      int namelen = NLEN;
      if (dpath) namelen += strlen(dpath)+1;

      /*      fprintf(stream,"new temp file construction\n");*/
      /* first, generate name - NOTE: THIS IS MEMORY THAT
	 MUST BE MANAGED BY THE CALLING UNIT 
	 Note that DATAPATH is used as a path, if it is defined, else the
	 current working directory
      */
      *name = (char *)usermalloc_(namelen*sizeof(char));
      memset(*name,'\0',namelen);
      if (dpath) {
	strcat(*name,dpath);
	if (strlen(dpath) && (dpath[strlen(dpath)-1] != '/')) strcat(*name,"/");
      }
      else strcat(*name,"./");
      strcat(*name,"tmp.XXXXXX");

      fd=mkstemp(*name);
      if (fd<0) {
	fprintf(stream,"Error: iwave_fopen\n");
	fprintf(stream,"failed to open temp file - error from mkstemp\n");
	return retfp;
      }

      /*      fprintf(stream,"iwave_fopen: temp file name = %s\n",*name); */

      /* open stream - always in w+ mode */
      if (!(retfp=fdopen(fd,"w+"))) {
	fprintf(stream,"Error: iwave_fopen\n");
	fprintf(stream,"failed to open stream - error from fdopen\n");
	return retfp;
      }

      /* set filestat struct params */
      *oldfpr = (struct filestat *)usermalloc_(sizeof(struct filestat));
      fpr=*oldfpr;
      fpr->fp = retfp;
      fpr->nm = (char *)usermalloc_((strlen(*name)+1)*sizeof(char));
      strcpy(fpr->nm,*name);
      fpr->pr = (char *)usermalloc_((strlen(proto)+1)*sizeof(char));
      strcpy(fpr->pr,proto);
      fpr->md = (char *)usermalloc_(3*sizeof(char));
      strcpy(fpr->md,"w+");
      fpr->istmp=1;
      fpr->inuse=1;
      fpr->nextfpr = (struct filestat *)NULL;

      /* determine length of proto file */
      fseeko(fph,0L,SEEK_END);
      plen=ftello(fph);

      /* now copy prototype file onto temp file */
      fseeko(fpr->fp,0L,SEEK_SET);
      fseeko(fph,0L,SEEK_SET);

      chunk=iwave_min(BLEN,plen);
      while (chunk>0) {

	/*	fprintf(stderr,"t loop chunk=%d\n",chunk);*/

	nr=fread(&(buf[0]),sizeof(char),chunk,fph);
	if (nr!=chunk) {
	  fprintf(stream,"Error: iwave_fopen\n");
	  fprintf(stream,"in initialization of temp file %s\n",fpr->nm);
	  fprintf(stream,"failed to read %d chars from proto file %s\n",chunk,fpr->pr);
	  return NULL;
	}
	nr=fwrite(&(buf[0]),sizeof(char),chunk,fpr->fp);
	if (nr!=chunk) {
	  fprintf(stream,"Error: iwave_fopen\n");
	  fprintf(stream,"in initialization of temp file %s\n",fpr->nm);
	  fprintf(stream,"failed to write %d chars to temp file %s\n",chunk,fpr->nm);
	  return NULL;
	}
	/* having read chunk words, there are plen-chunk left to read - update
	   plen, recompute chunk */
	plen=plen-chunk;
	chunk=iwave_min(BLEN,plen);
      }
    
      /* reset to begin-of-file to emulate w+ behaviour */
      fseeko(fpr->fp,0L,SEEK_SET);
      fseeko(fph,0L,SEEK_SET);
      
    }

    /* CASE IB: found unused temp file with matching prototype, 
       return its pointer, update name */
    
    else {
      /* found suitable unused temp file modeled on same proto. set 
	 inuse flag, copy pointer to return value, copy name to arg. */
      fpr->inuse=1;
      /* modification 23.08.11:

	 use freopen file to ensure that buffer is synched - otherwise
	 updated file may not be reflected in buffered data!

	 According to ISO standard, freopen closes the stream, then reopens it on 
	 the file in the first arg. The pointer in the third arg continues to point 
	 to the stream, and is the return value on success - on failure, returns a NULL.

	 ALWAYS OPEN WITH "r+" - this works because file has been written if we get 
	 to this point, in its complete length, from the prototype.
      */
#ifdef CONST_FP
      if ((sav = ftello(fpr->fp) < 0)) {
	fprintf(stream,"ERROR: iwave_fopen\n");
	fprintf(stream,"-- failed to record current offset before reopen\n");
	fflush(stream);
	return NULL;;
      }
#endif
      retfp = freopen(fpr->nm,"r+",fpr->fp);
      if (!retfp) {
	fprintf(stream,"ERROR: iwave_fopen\n");
	fprintf(stream,"-- failed to reopen stream on file %s mode r+ with same file pointer\n",fpr->nm);
	return retfp;
      }
#ifdef CONST_FP
      if (fseeko(retfp,sav,SEEK_SET))  {
	fprintf(stream,"ERROR: iwave_fopen\n");
	fprintf(stream,"-- failed to seek to saved offset\n");
	fflush(stream);
	return NULL;;
      }
#endif
      /* MEMORY WHICH MUST BE MANAGED BY CALLING UNIT */
      *name=(char *)usermalloc_((strlen(fpr->nm)+1)*sizeof(char));
      strcpy(*name,fpr->nm);

      /*      fprintf(stream,"iwave_fopen - re-use temp file %s\n",*name); */
    }

    /* in either case, return file pointer */
    return retfp;

    /*    fprintf(stderr,"end temp case\n");*/
  }

  /* CASE II: archival file requested */

  else {

    /*    fprintf(stderr,"begin archival case, name = %s\n",*name);*/

    /* match if 
       - filenames match, and
       - modes match, and
       - either no prototype, or prototypes match. and
       - inuse flag unset
       MOD OF 15.01.13: disregard inuse - return copy of pointer regardless
       MOD OF 15.01.13: r+ and w+ are mode wildcards for already opened file
       MOD OF 26.01.13: NULL proto is a prototype wildcard for already opened file 
    */
    oldfpr=&filestatlist;
    for (fpr=filestatlist; fpr != (struct filestat *)NULL; 
	 fpr = fpr->nextfpr) {
      /* break if filenames match */
      if (!(strcmp(*name,fpr->nm)) &&
	  (
	  !(strcmp(mode,fpr->md)) ||
	  !(strcmp(fpr->md,"r+")) ||
	  !(strcmp(fpr->md,"w+")) 
	   ) &&
	  (
	   ((proto) && (fpr->pr) &&
	   !(strcmp(proto,fpr->pr))) ||
	   /*	   ((!proto) && (!(fpr->pr)))*/
	   !proto
	   )
	  /* &&
	  (!fpr->inuse)
	  */
	  )
	break;
      oldfpr=&fpr->nextfpr;
    }
    
    /*    fprintf(stderr,"begin new archival file\n");*/
    
    /* CASE IIA didn't find one -have to build a new one */
    
    if (fpr==((struct filestat *)NULL)) {
      
      /* open stream */
      if (!(retfp=fopen(*name,mode))) {
	/* NOT NECESSARILY AN ERROR - so don't print! 
	fprintf(stream,"NOTE: iwave_fopen\n"); 
	fprintf(stream,"-- failed to open stream on file %s mode %s\n",*name,mode); 
	*/
	return retfp;
      }
      
      /* CASE IIA-1: if r mode, check against length of existing file */
      
      /* if r mode, check that file has same length as prototype - no
	 other tests make sense at this level. */
      /* 31.03.10: no error checking at this level for r access - defer
	 to DC classes */
      /*
	if (mode[0]=='r' && proto && fph) {
	
	fseeko(fph,0L,SEEK_END);
	plen=ftello(fph);
	fseeko(retfp,0L,SEEK_END);
	if (plen!=ftello(retfp)) {
	fprintf(stream,"Error: file %s does not have same length \n",*name);
	fprintf(stream,"as proto file %s\n",proto);
	fclose(retfp);
	  retfp=NULL;
	  return retfp;
	}
	fpr->inuse=1;
	

	}
      */

      /* CASE IIA-2: (write mode) copy prototype file to claim disk space,
         if a prototype is provided */
      /* note that this is done ONLY to claim disk space. The contents
	 of this file may be modified by users eg. DC classes */

      if (mode[0]=='w' && proto && fph) {

	/*fprintf(stderr,"old file case with prototype\n");*/

	/* determine length of proto file */
	fseeko(fph,0L,SEEK_END);
	plen=ftello(fph);

	fseeko(retfp,0L,SEEK_SET);
	fseeko(fph,0L,SEEK_SET);

	chunk=iwave_min(BLEN,plen);
	while (chunk>0) {

	  /*	  fprintf(stderr,"chunk=%d\n",chunk);*/

	  nr=fread(&(buf[0]),sizeof(char),chunk,fph);

	  if (nr!=chunk) {
	    fprintf(stream,"Error: iwave_fopen\n");
	    fprintf(stream,"in initialization of new archival file %s\n",fpr->nm);
	    fprintf(stream,"failed to read %d chars from proto file %s\n",chunk,fpr->pr);
	    fclose(retfp);
	    retfp=NULL;
	    return retfp;
	  }
	  nr=fwrite(&(buf[0]),sizeof(char),chunk,retfp);
	  if (nr!=chunk) {
	    fprintf(stream,"Error: iwave_fopen\n");
	    fprintf(stream,"in initialization of new archival file %s\n",fpr->nm);
	    fprintf(stream,"failed to write %d chars to temp file %s\n",chunk,fpr->nm);
	    fclose(retfp);
	    retfp=NULL;
	    return retfp;
	  }

	  /* having read chunk words, there are plen-chunk left to read - update
	     plen, recompute chunk */
	  plen=plen-chunk;
	  chunk=iwave_min(BLEN,plen);
	}

      }
	
      /* reset to begin-of-file to emulate fopen behaviour */
      fseeko(retfp,0L,SEEK_SET);
      if (fph) fseeko(fph,0L,SEEK_SET);

      /* set filestat struct params */
      *oldfpr = (struct filestat *)usermalloc_(sizeof(struct filestat));
      fpr=*oldfpr;
      fpr->fp = retfp;
      fpr->nm = (char *)usermalloc_((strlen(*name)+1)*sizeof(char));
      strcpy(fpr->nm,*name);
      if (proto) {
	fpr->pr = (char *)usermalloc_((strlen(proto)+1)*sizeof(char));
	strcpy(fpr->pr,proto);
      }
      else {
	fpr->pr = NULL;
      }
      fpr->md = (char *)usermalloc_((strlen(mode)+1)*sizeof(char));
      strcpy(fpr->md,mode);
      fpr->istmp=0;
      fpr->inuse=1;
      fpr->nextfpr = (struct filestat *)NULL;
    }

    else {

      /* CASE IIB: found already-opened archival *** or temp *** file */
      /* return to start-of-file, to imitate setup with freshly opened
	 file. */
      /* no, don't - very useful for the behaviour of iwave_fopen to be
	 DIFFERENT from that of fopen! if the file is already opened, 
	 simply return the pointer, at whatever position it was last
	 left, and leave any further positioning decisions to the calling
	 unit. For a sequential read-by-chunks through a file, that's 
	 exactly what is needed. If there is any risk that the pointer 
	 has moved since the last access by the unit in question, then 
	 of course that unit must possess overall positioning information.
	 Otherwise, no need! This way, the file pointers act like a 
	 static array.
      if ((fpr->md[0]=='w')||(fpr->md[0]=='r')) fseeko(fpr->fp,0L,SEEK_SET);
      */
      /* furthermore, if open for write, reopen/truncate */
      /* not clear that this is a good idea - hold off
	 if (fpr->md[0]=='w') fpr->fp=freopen(fpr->nm,fpr->md,fpr->fp);
      */

      /* UPSHOT: set in-use flag, simply pass the pointer back!!! */
      /* MODIFICATION 23.08.11: use freopen to force synch of
	 file data with buffered data 

	 if another fp is used to write to the file which is entirely
	 buffered, then buffer associated with this fp is not
	 updated. Must force update - only way to do this within stdio
	 is to close the stream, thus flushing and destroying all
	 buffers. Apparently, input buffers are not flushed to the
	 disk, but are discarded, as they should be.

	 According to ISO standard, freopen closes the stream, then
	 reopens it on the file in the first arg. The pointer in the
	 third arg continues to point to the stream, and is the return
	 value on success - on failure, returns a NULL.

	 ALWAYS REOPEN WITH r+ mode - this works because if file was
	 previously opened with w or w+ mode, then data was copied in
	 to fully populate it.
      */
#ifdef CONST_FP
      if ((sav = ftello(fpr->fp) < 0)) {
	fprintf(stream,"ERROR: iwave_fopen\n");
	fprintf(stream,"-- failed to record current offset before reopen\n");
	fflush(stream);
	return NULL;;
      }
#endif
      fpr->inuse++;      
      retfp = freopen(fpr->nm,"r+",fpr->fp);
      if (!retfp) {
	fprintf(stream,"ERROR: iwave_fopen\n");
	fprintf(stream,"-- failed to reopen stream on file %s mode r+ with same file pointer\n",fpr->nm);
	return retfp;
      }     
#ifdef CONST_FP
      if (fseeko(retfp,sav,SEEK_SET))  {
	fprintf(stream,"ERROR: iwave_fopen\n");
	fprintf(stream,"-- failed to seek to saved offset\n");
	fflush(stream);
	return NULL;;
      }
#endif
    }
  }

  /*
  if (retfp)
    fprintf(stderr,"opening fp=%x name=%s mode=%s istmp=%d inuse=%d\n",
      retfp,fpr->nm,fpr->md,fpr->istmp,fpr->inuse);
  */

  return retfp;

}

FILE * iwave_const_fopen(const char * name, 
			 const char * mode, 
			 const char * proto,
			 FILE * stream) {
  char * tmpname;
  FILE * fp = NULL;
  if (!name) {
    fprintf(stream,"Error: iwave_const_fopen\n");
    fprintf(stream,"called with null filename\n");
  }
  else {
    tmpname=(char *)usermalloc_((1+strlen(name))*sizeof(char));
    strcpy(tmpname,name);
    fp=iwave_fopen(&tmpname,mode,proto,stream);
    userfree_(tmpname);
  }

  /*
  if (!fp) {
    fprintf(stream,"NOTE: returning from iwave_const_fopen with null ptr\n");
    fprintf(stream,"NOTE: input params name=%s mode=%s proto=%s\n",name,mode,proto);
    fprintf(stream,"NOTE: state of file system:\n");
    iwave_fprintall(stream);
  }
  */

  return fp;
}

/* close merely resets the in-use flag */
void iwave_fclose(FILE * fp) {
  if (!fp) return;

  for (fpr=filestatlist; fpr != ((struct filestat *)NULL); 
       fpr = fpr->nextfpr) {
    if (fp==fpr->fp) break;
  }
  if (fpr!=((struct filestat *)NULL)) {
      (fpr->inuse)--;
      fpr->inuse=iwave_max(0,fpr->inuse);
    /*
    fprintf(stderr,"closing fp=%x name=%s mode=%s istmp=%d inuse=%d\n",
      fp,fpr->nm,fpr->md,fpr->istmp,fpr->inuse);
    */
  }
}

void iwave_fprintall(FILE * stream) {

  int iptr=0; int istmp=0; int inuse=0;

  fprintf(stream,"IWAVE file manager database\n");
  
  for (fpr=filestatlist; fpr != ((struct filestat *)NULL); fpr = fpr->nextfpr) {
    fprintf(stream,"FILE*=%p ",(void *)fpr->fp);
    if (fpr->nm) fprintf(stream,"name=%s ",fpr->nm);
    else fprintf(stream,"name=NULL ");
    if (fpr->md) fprintf(stream,"mode=%s ",fpr->md);
    else fprintf(stream,"mode=(none) ");
    if (fpr->pr) fprintf(stream,"proto=%s ",fpr->pr);
    else fprintf(stream,"proto=NULL ");
    fprintf(stream,"istmp = %d inuse = %d\n",fpr->istmp,fpr->inuse);
    if (fpr->istmp) istmp++;
    if (fpr->inuse) inuse++;
    iptr++;
  }
  fprintf(stream,"-- total FILE * = %d, open as tmp = %d, in use = %d\n",
	  iptr,istmp,inuse);
}

/* destructor */
void iwave_fdestroy() {
  /* workspace to hold temp copy of next pointer */
  struct filestat * tmpnext = filestatlist->nextfpr;

  for (fpr=filestatlist; fpr != ((struct filestat *)NULL); 
       fpr = tmpnext) {
    if (fpr->fp)  { fclose(fpr->fp); fpr->fp=NULL; }
#ifdef UNLINK_TMPS
    if (fpr->istmp && fpr->nm) unlink(fpr->nm);
#endif
    if (fpr->nm) { userfree_(fpr->nm); fpr->nm=NULL; }
    if (fpr->pr) { userfree_(fpr->pr); fpr->pr=NULL; }
    userfree_(fpr->md); fpr->md=NULL;
    fpr->istmp=0;
    fpr->inuse=0;
    tmpnext=fpr->nextfpr;
    userfree_(fpr);
    fpr=(struct filestat *)NULL;
  }
  filestatlist=(struct filestat *)NULL;
}

const char * iwave_getproto(const char * name) { 
  for (fpr=filestatlist; fpr != ((struct filestat *)NULL); fpr = fpr->nextfpr)
    if (!(strcmp(fpr->nm,name))) return fpr->pr;
  return NULL;
}
