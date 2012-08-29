#include "traceio.h"

/* axis indices throughout: 0=z, 1=x, 2=y */

/** helper function to determine whether an index tuple is within rarray */
int ingrid(int ndim, IPNT n, IPNT gs, IPNT itr) {
  int ret=0;
  if (((ndim==1) &&
       (itr[0]>gs[0]-1)     &&
       (itr[0]<gs[0]+n[0])) ||
      ((ndim==2) &&
       (itr[0]>gs[0]-1)     &&
       (itr[0]<gs[0]+n[0])  &&
       (itr[1]>gs[1]-1)     &&
       (itr[1]<gs[1]+n[1])) ||
      ((ndim==3) &&
       (itr[0]>gs[0]-1)     &&
       (itr[0]<gs[0]+n[0])  &&
       (itr[1]>gs[1]-1)     &&
       (itr[1]<gs[1]+n[1])  &&
       (itr[2]>gs[2]-1)     &&
       (itr[2]<gs[2]+n[2]))) ret=1;
  return ret;
}

/** helper function to determine whether a coord tuple is within grid */
/*int ringrid(int ndim, IPNT n, RPNT d, RPNT o, RPNT x,FILE * stream) {*/
int ringrid(int ndim, IPNT n, RPNT d, RPNT o, RPNT x) {
  int ret=0;
  RPNT e;
  int i;

  for (i=0;i<RARR_MAX_NDIM;i++) e[i]=o[i]+(n[i]-1)*d[i];

  if (((ndim==1) &&
       (x[0]>=o[0])     &&
       (x[0]<=e[0])) ||
      ((ndim==2) &&
       (x[0]>=o[0])     &&
       (x[0]<=e[0])  &&
       (x[1]>=o[1])     &&
       (x[1]<=e[1])) ||
      ((ndim==3) &&
       (x[0]>=o[0])     &&
       (x[0]<=e[0])  &&
       (x[1]>=o[1])     &&
       (x[1]<=e[1])  &&
       (x[2]>=o[2])     &&
       (x[2]<=e[2]))) ret=1;
  
  /*
    fprintf(stream,"x=[%e %e] ",x[0],x[1]);
    fprintf(stream,"o=[%e %e] ",o[0],o[1]);
    fprintf(stream,"e=[%e %e] ",e[0],e[1]);
    fprintf(stream,"ringrid=%d\n",ret);
  */      

  return ret;
}

#ifdef IWAVE_USE_MPI 

/** Courtesy of Bee Bednar: light modification of his Build_MPI_SegY_tr.c 
    15.01.09: include file offset as first word, stored as off_t. */

void buildMPI_OFFSEGY(MPI_Datatype * ptr, int ns) {

  int i;
  const int     nblks = 11;
  int           lengs[nblks];
  MPI_Datatype  types[nblks];
  MPI_Aint      disps[nblks];
  size_t        sizes[nblks];
  
  lengs[0]=1;  types[0]=MPI_LONG_LONG;   sizes[0]=sizeof(off_t);
  lengs[1]=7;  types[1]=MPI_INT;    sizes[1]=sizeof(int);
  lengs[2]=4;  types[2]=MPI_SHORT;  sizes[2]=sizeof(short);
  lengs[3]=8;  types[3]=MPI_INT;    sizes[3]=sizeof(int);
  lengs[4]=2;  types[4]=MPI_SHORT;  sizes[4]=sizeof(short);
  lengs[5]=4;  types[5]=MPI_INT;    sizes[5]=sizeof(int);
  lengs[6]=46; types[6]=MPI_SHORT;  sizes[6]=sizeof(short);
  lengs[7]=6;  types[7]=MPI_FLOAT;  sizes[7]=sizeof(float);
  lengs[8]=1;  types[8]=MPI_INT;    sizes[8]=sizeof(int);
  lengs[9]=16; types[9]=MPI_SHORT;  sizes[9]=sizeof(short);
  lengs[10]=ns; types[10]=MPI_FLOAT;  sizes[10]=sizeof(float);

  disps[0]=0;
  for (i=1;i<nblks;i++) 
    disps[i]=disps[i-1]+lengs[i-1]*sizes[i-1];
  /* deprecated */
  MPI_Type_struct(nblks,lengs,disps,types,ptr); 
  /* supposedly more up-to-date */
  /*  MPI_Type_create_struct(nblks,lengs,disps,types,ptr); */
  MPI_Type_commit(ptr);

}

#endif

/* compute group index, set return value if irec >= last index
   for this group, else increment irec. Note - only for group 
   root 

   20.11.10: group id computation now in iwave, only responsibility
   here is to compute first and last record indices for current group
*/

void calc_group(int * first, int * last, 
		int nrec) {

  div_t res;
  int sg;

  /* retrieve MPI info */
  int ng=retrieveNumGroups();
  int g =retrieveGroupID();

  res=div(nrec,ng);
  sg=res.quot;

  /* 13.03.11: comment out 
  if (res.rem) sg++;

  *first=g*sg;
  *last =iwave_min((g+1)*sg-1,nrec-1);
  */
  
  /* D.S.: 13.03.11*/
  if(g >= res.rem){
    *first = res.rem*(sg+1) + (g - res.rem)*sg;
    *last = *first + sg - 1;
  }
  else {
    *first = g*(sg+1);
    *last = *first + sg;
  }
  
}

int traceserver_init(FILE ** fpin, char * fin,
		     FILE ** fpout, char * fout,
		     int * irec, int * xrec, 
		     int * first, int * last,
		     int * nrec, 
		     int * ntr, off_t * recoff, RPNT * src,
#ifdef IWAVE_USE_MPI
		     MPI_Datatype * p,
#endif
		     float tol,
		     FILE * stream) {
  int ir;
  Value val;
  ireal scalco;
  ireal scalel;

  float sz,sx,sy;
  float ds;
  segy tr;

  short ns;
  int ins;
  int err=0;

  off_t offtmp;     /* buffer for offset of next trace */

  int rkw;          /* global rank */
  int rkl;          /* local rank */
  int szw;          /* global size */
  int szl;          /* local size */

#ifdef IWAVE_USE_MPI
  MPI_Comm cmw;     /* global communicator */
  MPI_Comm cml;     /* local communicator */
#endif

  rkw=retrieveGlobalRank();
  rkl=retrieveRank();
  szw=retrieveGlobalSize();
  szl=retrieveSize();

#ifdef IWAVE_USE_MPI
  cml=retrieveComm();
  cmw=retrieveGlobalComm();
#endif

  if (rkl==0) {

    /* fprintf(stderr,"TRACEIO::TRACESERVER_INIT: start with fin = %s, fout= %s \n",fin, fout);*/
 
    *fpin=NULL;
    *fpout=NULL;
    
    /* open files on i/o process of each group
     NOTE: 010510 use of file mgr requres seeks, since already-opened
     files may have pointers positioned anywhere in file.
    */
    /* input file must exist - acts as model for output */
#ifdef IWAVE_USE_FMGR
    /* iwave_fopen version */
    if (!(*fpin=iwave_const_fopen(fin,"r",NULL,stream))) {
      fprintf(stream,"Error: traceserver_init\n");
      fprintf(stream,"failed to open input file\n");
      err=E_FILE;
      fflush(stream);
#ifdef IWAVE_USE_MPI
      MPI_Abort(cmw,err);
#else
      return err;
#endif
    }
    
#else
    /* stdio version */
    if (!(*fpin=fopen(fin,"r"))) {
      fprintf(stream,"Error: traceserver_init\n");
      fprintf(stream,"failed to open input file\n");
      err=E_FILE;
      fflush(stream);
#ifdef IWAVE_USE_MPI
      MPI_Abort(cmw,err);
#else
      return err;
#endif
    }
#endif /* end file mgr branch */
    
    if (fseeko(*fpin,0L,SEEK_SET)) {
      fprintf(stream,"Error: traceserver_init\n");
      fprintf(stream,"failed to reset input file\n");
      err=E_FILE;
      fflush(stream);
#ifdef IWAVE_USE_MPI
      MPI_Abort(cmw,err);
#else
      return err;
#endif
    }
    
    /* output file may or may not exist. If it does, then either it
       has the same length as the input file, or it is closed and
       reopened with the "w" flag, in particular truncated. In the
       latter case the input is copied onto the output, to create a
       file of the same length. We also enable reading of the output
       file, as this may be necessary to accumulate output traces over
       subdomain boundaries.
    */
#ifdef IWAVE_USE_FMGR
    /* iwave_fopen version - fin is the proto file. In w+ case,
       iwave_fopen copies fin to fout, so copy stuff later in this
       function can be skipped. In r+ case, fin, fout must have same
       file length - nothing else is checked here.
    */
    if (!(*fpout=iwave_fopen(&fout,"r+",fin,stream))) { 
      if (!(*fpout=iwave_fopen(&fout,"w+",fin,stream))) {
	fprintf(stream,"Error: traceserver_init\n");
	fprintf(stream,"failed to open output file (w+)\n");
	err=E_FILE;
	fflush(stream);
#ifdef IWAVE_USE_MPI
	MPI_Abort(cmw,err);
#else
	return err;
#endif
      }
    }

#else
    /* stdio version */
    if (!(*fpout=fopen(fout,"r+"))) { 
      if (!(*fpout=fopen(fout,"w+"))) {
	fprintf(stream,"Error: traceserver_init\n");
	fprintf(stream,"failed to open output file (w+)\n");
	err=E_FILE;
	fflush(stream);
#ifdef IWAVE_USE_MPI
	MPI_Abort(cmw,err);
#else
	return err;
#endif
      }
    }
#endif // end file mgr branch

    if (fseeko(*fpout,0L,SEEK_SET)) {
      fprintf(stream,"Error: traceserver_init\n");
      fprintf(stream,"failed to reset output file\n");
      err=E_FILE;
      fflush(stream);
#ifdef IWAVE_USE_MPI
      MPI_Abort(cmw,err);
#else
      return err;
#endif
    }
    
  }

  /* initialize src array, number of traces on all processes */
  for (ir=0;ir<MAX_RECS;ir++) {
    RASN(src[ir],RPNT_0);
    ntr[ir]=0;
  }

  /* on root process of global comm, read input file to determine record
     info. If necessary also copy traces to output file.*/

  if (rkw==0) {

    /* in this case, local rank must also be zero (check!), and 
       thus files open in this process */
    if (rkl!=0) {
      fprintf(stream,"Error: traceserver_init\n");
      fprintf(stream,"world root must be root in its group\n");
      fprintf(stream,"world rank = %d, group rank = %d\n",rkw,rkl);
      err=E_OTHER;
      fflush(stream);
#ifdef IWAVE_USE_MPI
      MPI_Abort(cmw,err);
#else
      return err;
#endif
    }

    /* read first trace */
    ir=0;
    *nrec=0;
    offtmp=0;
    recoff[ir]=0;
    ntr[ir]=0;
    
    if (fgettr(*fpin,&tr)) {
      /* tentative begin of next rec */
      offtmp=ftello(*fpin);
      
      /* at least one record, at least one trace in it - presumed
	 not to exceed bounds for either! */
      (*nrec)++;
      ntr[ir]++;
      
      /* record ns for first trace, both as short and as int.*/
      ns=tr.ns;
      ins=tr.ns;
      
      /* read scale information */
      gethdval(&tr,"scalco",&val);
      scalco = vtof(hdtype("scalco"),val);
      gethdval(&tr,"scalel",&val);
      scalel = vtof(hdtype("scalel"),val);
      
      /* read source position */
      gethdval(&tr,"sx",&val);
      src[ir][1] = vtof(hdtype("sx"),val);
      if (scalco > 0) { src[ir][1] *=  scalco; }
      if (scalco < 0) { src[ir][1] /= -scalco; }
      
      gethdval(&tr,"sy",&val);
      src[ir][2] = vtof(hdtype("sy"),val);
      if (scalco > 0) { src[ir][2] *=  scalco; }
      if (scalco < 0) { src[ir][2] /= -scalco; }
      
      gethdval(&tr,"selev",&val);
      src[ir][0] = vtof(hdtype("selev"),val);
      if (scalel > 0) { src[ir][0] *=  scalel; }
      if (scalel < 0) { src[ir][0] /= -scalel; }
      src[ir][0]=-src[ir][0];
      
      /* set tracr, write trace to output file */
      tr.tracr=ntr[ir];
      fputtr(*fpout,&tr);
      
      /* trace reading loop */
      while (fgettr(*fpin,&tr)) {
	
	if (ntr[ir]>MAX_TRACES-1) {
	  fprintf(stream,"Error: traceserver_init\n");
	  fprintf(stream,"read too many traces in record %d\n",ir);
	  fprintf(stream,"exceeded MAX_TRACES set in trace/include/traceio.h\n");
	  fprintf(stream,"current value of MAX_TRACES = %d\n",MAX_TRACES);
	  err=E_FILE;
	  fflush(stream);
#ifdef IWAVE_USE_MPI
	  /* MPI_Abort(cml,err);*/
	  MPI_Abort(cmw,err);
#else
	  return err;
#endif
	} 
        
	/* compare ns with first trace - MUST BE SAME */
	if (tr.ns != ns) {
	  fprintf(stream,"Error: traceserver_init\n");
	  fprintf(stream,"ns on trace %d rec %d = %d\n",ntr[ir],ir,tr.ns);
	  fprintf(stream,"does not match ns on first trace = %d\n",ns);
	  err=E_FILE;
	  fflush(stream);
#ifdef IWAVE_USE_MPI
	  /* MPI_Abort(cml,err);*/
	  MPI_Abort(cmw,err);
#else
	  return err;
#endif	
	}

	/* extract source position */
	 
	/* read scale information */
	gethdval(&tr,"scalco",&val);
	scalco = vtof(hdtype("scalco"),val);
	gethdval(&tr,"scalel",&val);
	scalel = vtof(hdtype("scalel"),val);
	
	/* read source position */
	gethdval(&tr,"sx",&val);
	sx = vtof(hdtype("sx"),val);
	if (scalco > 0) { sx *=  scalco; }
	if (scalco < 0) { sx /= -scalco; }
	
	gethdval(&tr,"sy",&val);
	sy = vtof(hdtype("sy"),val);
	if (scalco > 0) { sy *=  scalco; }
	if (scalco < 0) { sy /= -scalco; }
	
	gethdval(&tr,"selev",&val);
	sz = vtof(hdtype("selev"),val);
	if (scalel > 0) { sz *=  scalel; }
	if (scalel < 0) { sz /= -scalel; }
	sz=-sz;
	
	/* compare to previous */
	ds = fmaxf(fabsf(sz-src[ir][0]),
		   fmaxf(fabsf(sx-src[ir][1]),
			 fabsf(sy-src[ir][2])));
	
	if (ds<tol) {
	  /* found another trace in current record */
	  ntr[ir]++;
	  /* the following line added by D.S. 10.10.10 */
	  offtmp=ftello(*fpin);
	}
	else {
	  /* found first trace of next record */
	  ir++;
	  if (ir>MAX_RECS-1) {
	    fprintf(stream,"Error: traceserver_init\n");
	    fprintf(stream,"read too many records\n");
	    fprintf(stream,"exceeded MAX_RECS set in trace/include/traceio.h\n");
	    fprintf(stream,"current value of MAX_RECS = %d\n",MAX_RECS);
	    err = E_FILE;
	    fflush(stream);
#ifdef IWAVE_USE_MPI
	  /* MPI_Abort(cml,err);*/
	  MPI_Abort(cmw,err);
#else
	    return err;
#endif
	  }

	  //	  recoff[ir]=offtmp;
	  recoff[ir]=recoff[ir-1]+ntr[ir-1]*(HDRBYTES+ns*sizeof(float));
	  ntr[ir]=1;
	  (*nrec)++;
	  src[ir][0]=sz;
	  src[ir][1]=sx;
	  src[ir][2]=sy;
	  /* tentative begin of next rec */
	  /* comment out by D.S. 10.10.10 (REF L436) */
	  /* offtmp=ftello(*fpin); */
	}
	/* set tracr, write trace to output */
	tr.tracr=ntr[ir];
	fputtr(*fpout,&tr);
      }
      /* set tracr, write trace to output */
      
      /* (D.S. 06.10.10) move the following two lines into the while-loop above */
      /*
	tr.tracr=ntr[ir];
	fputtr(*fpout,&tr);
      */
    }
    else {
      fprintf(stream,"Error: traceserver_init\n");
      fprintf(stream,"failed to read first trace on input file %s\n",fin);
      err=E_FILE;
      fflush(stream);
#ifdef IWAVE_USE_MPI
      /* MPI_Abort(cml,err);*/
      MPI_Abort(cmw,err);      
#else
      return err;
#endif
    } 
  }
  
#ifdef IWAVE_USE_MPI

  if(szw>1) {
    /* broadcast record info */
    MPI_Bcast(nrec,1,MPI_INT,0,cmw);
    MPI_Bcast(ntr,*nrec,MPI_INT,0,cmw);
    /* CAUTION: the following line assumes that the PEs in the comm all
       agree in word structure, i.e. that an off_t is the same sequence
     of bytes in all machines */
    MPI_Bcast((char *)(&(recoff[0])),(*nrec)*sizeof(off_t),MPI_CHAR,0,cmw);
    for (ir=0;ir<*nrec;ir++) 
      MPI_Bcast(src[ir],RARR_MAX_NDIM,MPI_FLOAT,0,cmw);
    
    /* broadcast number of samples as int, so that MPI_Datatype can be initialized */
    MPI_Bcast(&ins,1,MPI_INT,0,cmw);
  }
  
  /* construct MPI datatype in every process */
  buildMPI_OFFSEGY(p,ins);

#endif

  /* at this point we have initialized nrec, ntr, and recoff, and
     broadcast all of these. so every process can determine its own
     first and last record numbers
  */
  calc_group(first,last,*nrec);

  /* seek to first record (local to each group) */

  if (rkl==0) {
    fflush(*fpout);
    err=err||fseeko(*fpin,recoff[*first],SEEK_SET);
    err=err||fseeko(*fpout,recoff[*first],SEEK_SET);
    if (err) {
      fprintf(stream,"ERROR: traceserver_init\n");
      fprintf(stream,"seek failed\n");
      fflush(stream);
#ifdef IWAVE_USE_MPI
      MPI_Abort(cmw,err);
#else
      return err;
#endif
    }
  }
  
  /* initial set for irec, xrec */
  *irec=*first;
  *xrec=*first;
  
  /* which will surely be 0 */
  return err;
}
  
int traceserver_rec(int * irec,
		    int * xrec,
		    int last,
		    FILE * stream) {
  int err=0;

  int lsize;         /* size of the local communicator */
  int lrank;         /* rank in the local communicator  */

#ifdef IWAVE_USE_MPI
  MPI_Comm wcomm;    /* global communicator */
  MPI_Comm lcomm;    /* local communicator */

  /* local variables for mpi info */
  wcomm=retrieveGlobalComm();
  lcomm=retrieveComm();
#endif

  lrank=retrieveRank();
  lsize=retrieveSize();


  /* only on rank 0  */

  if (lrank==0) {
    *irec=*xrec;
    (*xrec)++;    
    if (*irec > last) {
      err=W_EOF;
      fprintf(stream,"ERROR: traceserver_rec\n");
      fprintf(stream,"record index out of bounds\n");
      fflush(stream);
#ifdef IWAVE_USE_MPI
      MPI_Abort(wcomm,err);
#else
      return err;
#endif
    }
  }    
    
#ifdef IWAVE_USE_MPI

  /* broadcast new record index, error */
  if (lsize>1) {
    if (MPI_SUCCESS != MPI_Bcast(irec,1,MPI_INT,0,lcomm)) {
      fprintf(stream,"ERROR: traceserver_rec\n");
      fprintf(stream,"failed to broadcast record index\n");
      fflush(stream);
      MPI_Abort(wcomm,E_OTHER);
    }
  
  /* 30.08.09: must keep xrec coherent across comm as well, because 
     it is used in test BEFORE this function is called in traceterm_init
     see notes in asg driver
  */
    if (MPI_SUCCESS != MPI_Bcast(xrec,1,MPI_INT,0,lcomm)) {
      fprintf(stream,"ERROR: traceserver_rec\n");
      fprintf(stream,"failed to broadcast next record index\n");
      fflush(stream);
      MPI_Abort(wcomm,E_OTHER);
    }
  }
  /*  fprintf(stream,"tr_rec: rkw=%d rk=%d ig=%d\n",rkw,rk,ig);*/
#endif
  /*  fprintf(stream,"tr_rec: last=%d irec=%d\n",last,*irec);*/
  return err;
}

int traceserver_seek(FILE * fp,
		     off_t * n) {
  int res=0;

  int lrank;         /* rank in the local communicator  */

  /* local variables for mpi info */
  lrank=retrieveRank();

#ifdef IWAVE_USE_MPI 
  if (lrank==0)  {
#endif
    res=fseeko(fp,*n,SEEK_SET);
#ifdef IWAVE_USE_MPI 
  }
#endif
  return res;
}

int traceserver_get(FILE * fp,
#ifdef IWAVE_USE_MPI
		    MPI_Datatype *p,
#endif
		    offsegy * otr) {

  int err=0;

  int lsize;         /* size of the local communicator */
  int lrank;         /* rank in the local communicator  */

#ifdef IWAVE_USE_MPI
  MPI_Comm wcomm;    /* global communicator */
  MPI_Comm lcomm;    /* local communicator */

  /* local variables for mpi info */
  wcomm=retrieveGlobalComm();
  lcomm=retrieveComm();
#endif
  lrank=retrieveRank();
  lsize=retrieveSize();

  if (lrank==0) {

    otr->m=ftello(fp);
    if (!fgettr(fp,&(otr->tr))) {
      fprintf(stderr,"PANIC: traceserver_get\n");
      fprintf(stderr,"failed to read trace on input unit\n");
#ifdef IWAVE_USE_MPI      
      MPI_Abort(wcomm,W_EOF);
#else
      return W_EOF;
#endif
    }
  }

  /* broadcast trace */
#ifdef IWAVE_USE_MPI
  if(lsize>1)
    MPI_Bcast(otr,1,*p,0,lcomm);
#endif

  return err;
}

int traceserver_put(FILE * fp,
#ifdef IWAVE_USE_MPI
		    MPI_Datatype *p,
#endif
		    offsegy * otr) {

  int err=0;

  int lsize;         /* size of the local communicator */
  int lrank;         /* rank in the local communicator  */
#ifdef IWAVE_USE_MPI
  MPI_Comm wcomm;    /* global communicator */
  MPI_Comm lcomm;    /* local communicator */

  /* local variables for mpi info */
  wcomm=retrieveGlobalComm();
  lcomm=retrieveComm();
#endif
  lrank=retrieveRank();
  lsize=retrieveSize();

#ifdef IWAVE_USE_MPI

  MPI_Status mpistat;
  //  int i;
  if (lsize>1) {
    if (MPI_SUCCESS != MPI_Send(otr,1,*p,0,lrank,lcomm)) {
      fprintf(stderr,"PANIC: traceserver_put\n");
      fprintf(stderr,"failed to send trace from non-root\n");
      MPI_Abort(wcomm,E_OTHER);
    }
    if (lrank==0) {
      if (MPI_SUCCESS != MPI_Recv(otr,1,*p,MPI_ANY_SOURCE,MPI_ANY_TAG,lcomm,&mpistat)) {
	fprintf(stderr,"PANIC: traceserver_put\n");
	fprintf(stderr,"failed to recv trace at root\n");
	MPI_Abort(wcomm,E_OTHER);
      }
    }
  }

  /*
  for (i=0;i<retrieveSize();i++) {

    if (retrieveRank()==i) {
      if (MPI_SUCCESS != MPI_Send(otr,1,*p,0,i,retrieveComm())) {
	fprintf(stderr,"PANIC: traceserver_put\n");
	fprintf(stderr,"failed to send trace from non-root\n");
	MPI_Abort(retrieveGlobalComm(),E_OTHER);
      }
    }

    if (retrieveRank()==0) {
      if (MPI_SUCCESS != MPI_Recv(otr,1,*p,i,i,retrieveComm(),&mpistat)) {
	fprintf(stderr,"PANIC: traceserver_put\n");
	fprintf(stderr,"failed to recv trace at root\n");
	MPI_Abort(retrieveGlobalComm(),E_OTHER);
      }
    }
  }
  */
#endif

  if (lrank==0) {
    
      if ((err=fseeko(fp,otr->m,SEEK_SET))) {
	  fprintf(stderr,"Error: traceserver_put from fseeko\n");
#ifdef IWAVE_USE_MPI
	  MPI_Abort(wcomm,err);
#else
	  return err;
#endif
      }
    fputtr(fp,&(otr->tr));
    
  }
  
  return err;
}

int construct_tracegeom(tracegeom * tg,
			char * hdr,
			char * data,
			float tol,
			FILE * stream) {

  int err=0;

#ifdef IWAVE_USE_MPI
  MPI_Comm wcomm=retrieveGlobalComm();  /*global communicator*/
#endif

  /*Igor - destroy previous one*/
  /*destroy_tracegeom(tg); */
  /*WWS 14.12.08 - clear out buffers, do not destroy tg - this object depends on it!*/
  /*  fprintf(stderr,"TRACEIO::CONSTRUCT_TRACEGEOM hdr=%s data=%s\n",hdr,data);*/
  /*  fprintf(stderr,"TRACEIO::CONSTRUCT_TRACEGEOM->setnull\n");*/
  setnull_tracegeom(tg);

  /*  fprintf(stderr,"TRACEIO::CONSTRUCT_TRACEGEOM->allocate source arrays\n");*/

  if (!(tg->ntr=(int*)usermalloc_(MAX_RECS*sizeof(int)))) {
    fprintf(stream,"Error: contruct_tracegeom\n");
    fprintf(stream,"failed to allocate ntr array\n");
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_ALLOC);
#else
    return E_ALLOC;
#endif
  }

  if (!(tg->recoff=(off_t*)usermalloc_(MAX_RECS*sizeof(off_t)))) {
    fprintf(stream,"Error: contruct_tracegeom\n");
    fprintf(stream,"failed to allocate recoff array\n");
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_ALLOC);
#else
    return E_ALLOC;
#endif
  }

  if (!(tg->src=(RPNT*)usermalloc_(MAX_RECS*sizeof(RPNT)))) {
    fprintf(stream,"Error: contruct_tracegeom\n");
    fprintf(stream,"failed to allocate ntr array\n");
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_ALLOC);
#else
    return E_ALLOC;
#endif
  }

  /* D.S. 06.10.10 */
  /*
  fprintf(stderr,"TRACEIO::CONSTRUCT_TRACEGEOM->traceserver_init (before call), err=%d \n",err);
  fprint_tracegeom(tg, stderr);
  scanf("%d",&err);
  */

  err=traceserver_init(&(tg->fpin),hdr,
		       &(tg->fpout),data,
		       &(tg->irec),&(tg->xrec),
		       &(tg->first),&(tg->last),
		       &(tg->nrec),
		       tg->ntr,tg->recoff,tg->src,
#ifdef IWAVE_USE_MPI
		       &(tg->p),
#endif
		       tol,
		       stream);

  /* D.S. 06.10.10 */
  /*
    fprintf(stderr,"TRACEIO::CONSTRUCT_TRACEGEOM->traceserver_init (after call), err=%d \n",err);
    fprint_tracegeom(tg, stderr);
    scanf("%d",&err);
  */

  if (err) {
    fprintf(stream,"Error: contruct_tracegeom\n");
    fprintf(stream,"from traceserver_init\n");
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,err);
#else
    return err;
#endif
  }

  /*  fprintf(stderr,"TRACEIO::CONSTRUCT_TRACEGEOM->allocate receiver arrays\n");*/

  if (!(tg->ig=(IPNT*)usermalloc_(MAX_TRACES*sizeof(IPNT)))) {
    fprintf(stream,"Error: contruct_tracegeom\n");
    fprintf(stream,"failed to allocate ig array\n");
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_ALLOC);
#else
    return E_ALLOC;
#endif
  }

  if (!(tg->rg=(RPNT*)usermalloc_(MAX_TRACES*sizeof(RPNT)))) {
    fprintf(stream,"Error: contruct_tracegeom\n");
    fprintf(stream,"failed to allocate rg array\n");
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_ALLOC);
#else
    return E_ALLOC;
#endif
  }

  if (!(tg->tracl=(int*)usermalloc_(MAX_TRACES*sizeof(int)))) {
    fprintf(stream,"Error: contruct_tracegeom\n");
    fprintf(stream,"failed to allocate tracl array\n");
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_ALLOC);
#else
    return E_ALLOC;
#endif
  }

  if (!(tg->tracr=(int*)usermalloc_(MAX_TRACES*sizeof(int)))) {
    fprintf(stream,"Error: contruct_tracegeom\n");
    fprintf(stream,"failed to allocate tracr array\n");
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_ALLOC);
#else
    return E_ALLOC;
#endif
  }

  if (!(tg->tracf=(int*)usermalloc_(MAX_TRACES*sizeof(int)))) {
    fprintf(stream,"Error: contruct_tracegeom\n");
    fprintf(stream,"failed to allocate tracf array\n");
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_ALLOC);
#else
    return E_ALLOC;
#endif
  }

  if (!(tg->fldr=(int*)usermalloc_(MAX_TRACES*sizeof(int)))) {
    fprintf(stream,"Error: contruct_tracegeom\n");
    fprintf(stream,"failed to allocate fldr array\n");
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_ALLOC);
#else
    return E_ALLOC;
#endif
  }

  if (!(tg->troff=(off_t*)usermalloc_(MAX_TRACES*sizeof(off_t)))) {
    fprintf(stream,"Error: contruct_tracegeom\n");
    fprintf(stream,"failed to allocate troff array\n");
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_ALLOC);
#else
    return E_ALLOC;
#endif
  }

  /*  fprintf(stderr,"TRACEIO::CONSTRUCT_TRACEGEOM->allocate arrays err=%d\n",err);*/
  return err;

}



int init_tracegeom(tracegeom * tg,
		   RPNT og,
		   IPNT n, RPNT d, RPNT o,
		   IPNT axord,
		   int order,
		   float dt, int ndim,
		   int usernt, float usert0, 
		   int initbuf,
		   FILE * stream) {

  /******************* BEGIN LOCAL DECLARATIONS **********************/

  int i;                        /* counter */
  offsegy otr;                  /* segy workspace aug. with off_t */
  Value val;                    /* header val workspace */

  float tmpx;                   /* temp storage for recvr coords */
  float tmpy;
  float tmpz;
  float tmpsx;                  /* temp storage for source coords */
  float tmpsy;
  float tmpsz;
  float tmpdt;                  /* temp storage for dt */
  int tmpnt;                    /* temp storage for number of samples */
  float tmpt0;                  /* temp storage for first sample time */

  float ozg;                    /* GLOBAL grid origin coordinates */ 
  float oxg; 
  float oyg;
  int nz; int nx; int ny;       /* LOCAL grid parameters */
  float oz; float ox; float oy;
  float dz; float dx; float dy;
  RPNT ol;                      /* parameters for extended grid */
  IPNT nl; 
  RPNT x;                       /* workspace to pass recv coords to ringrid */
  int marg;                     /* margin implied by sample order - added to grid
				   axis on each end */
  int * init;                   /* flag array for convenience on re-reading */
  int iinit;                    /* initialized trace counter */
  int itr;                      /* trace counter */
  float * work;                 /* workspace for spline adjoint interpolation */
  int wlen;                     /* length of allocated spline workspace */
  int iend=1;                   /* endpoint code for splines (linear) */
  int err=0;                    /* error flag */

  int rk=0;                     /* rank in local MPI comm */

  MPI_Comm wcomm;               /* global communicator */
  /******************* END LOCAL DECLARATIONS **********************/

  rk=retrieveRank();
  wcomm=retrieveGlobalComm();

  /* sanity check */
  if (tg->nrec==0) {
    fprintf(stream,"Error: init_tracegeom\n");
    fprintf(stream,"number of records = 0 on call\n");
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_BADINPUT);
#else
    return E_BADINPUT;
#endif
  }
 
  /* allocate init array */
  if (!(init=(int *)usermalloc_(MAX_TRACES*sizeof(int)))) {
    fprintf(stream,"Error: init_tracegeom\n");
    fprintf(stream,"failed to allocate %d ints in init array\n",MAX_TRACES);
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_ALLOC);
#else
    return E_ALLOC;
#endif
  }
    
  /* obtain correct record number - if no more records, set return flag */
  err=traceserver_rec(&(tg->irec),&(tg->xrec),tg->last,stream);
  /*
  fprintf(stream,"tg_init: irec=%d err=%d\n",tg->irec,err);
  fflush(stream);
  */

  /* note that this is NOT an error! */
  if (err) {
    fprintf(stream,"NOTE: init_tracegeom from traceserver_rec\n");
    fprintf(stream,"no more records!\n");
    fflush(stream);
    return err;
  }

  /* read trace info from fpout rather than fpin
     ONLY because tracr is GUARANTEED TO BE SET CORRECTLY
     in fpout - done in traceserver_init. Otherwise all 
     header and offset info is the same between the two.
  */

  /*  if (rk==0) fprintf(stderr,"in init: call tr_seek\n");*/
  //  fprintf(stream,"rk=%d call traceserver_seek\n",retrieveRank()); 
  /* seek to correct location */
  err=traceserver_seek(tg->fpin,&(tg->recoff[tg->irec]));
  // err=traceserver_seek(tg->fpout,&(tg->recoff[tg->irec]));
  //  fprintf(stream,"rk=%d return traceserver_seek\n",retrieveRank());
  fflush(stream);
  if (err) {
    fprintf(stream,"Error: init_tracegeom from traceserver_seek\n");
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,err);
#else
    return err;
#endif
  }
  /* common initializations for record */

  /*  if (rk==0) fprintf(stderr,"in init: initializations\n");*/
  /* initialize flag array FOR RECORD */
  for (itr=0;itr<MAX_TRACES;itr++) init[itr]=0;
  /* initialize trace counter FOR RECORD */
  itr=0;
   
  /* initialize grid parameters - default outside dimn */
  nz=1; dz=1.0; oz=0.0; ozg = 0.0;
  nx=1; dx=1.0; ox=0.0; oxg = 0.0;
  ny=1; dy=1.0; oy=0.0; oyg = 0.0;
  
  /*  fprintf(stream,"axord: 0=%d 1=%d 2=%d\n",axord[0],axord[1],axord[2]);*/

  /* sanity check axis ordering array */
  if (ndim==1) {
    if (axord[0]!=0) { 
      fprintf(stream,"Error: [init_tracegeom] axis order munged\n");
      fflush(stream);
#ifdef IWAVE_USE_MPI
      MPI_Abort(wcomm,E_BADARRINDEX);
#else
      return E_BADARRINDEX; 
#endif
    }
  }
  else if (ndim==2) {
    if (axord[0] < 0 || 
	axord[1] < 0 || 
	axord[0]+axord[1] != 1 ||
	axord[0]*axord[1] != 0) { 
      fprintf(stream,"Error: axis order munged\n");
      fflush(stream);
#ifdef IWAVE_USE_MPI
      MPI_Abort(wcomm,E_BADARRINDEX);
#else
      return E_BADARRINDEX; 
#endif
    }
  }
  else if (ndim==3) {
    if (axord[0] < 0 || 
	axord[1] < 0 || 
	axord[2] < 0 ||
	axord[0]+axord[1]+axord[2] != 3 ||
	axord[0]*axord[1]+axord[0]*axord[2]+axord[1]*axord[2] != 2 ||
	axord[0]*axord[1]*axord[2] != 0) { 
      fprintf(stream,"Error: axis order munged\n");
      fflush(stream);
#ifdef IWAVE_USE_MPI
      MPI_Abort(wcomm,E_BADARRINDEX);
#else
      return E_BADARRINDEX; 
#endif
    }
  }
  else {
    fprintf(stream,"Error: infeasible ndim=%d\n",ndim); 
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_OTHER);
#else
    return E_OTHER; 
#endif
  }

  /* in computing local grid params for testing, include points
     at which interpolation will add contributions from the grid
     passed as argument - this means adding marg points on each
     end of each axis, shifting the origin left by marg*d[i] on 
     each axis. For symmetric interpolation (assumed), marg = 
     (order+1)/2.
     
     Added 18.04.09: even for order 0 sampling, permit search in 
     immediate nbhd of gridpoint so marg should always be >= 1. 
     Otherwise 
  */
  /*  marg=iwave_max(1,(order+1)/2);*/
  marg=(order+1)/2;
  if (marg<0 || marg>2) {
    fprintf(stream,"Error: infeasible order=%d\n",order);
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_OTHER);
#else
    return E_OTHER; 
#endif
  }

  RASN(ol,RPNT_0);
  IASN(nl,IPNT_1);

  if (ndim > 0) { 
    nz=n[axord[0]]; 
    dz=d[axord[0]]; 
    oz=o[axord[0]]; 
    ozg=og[axord[0]];
    nl[axord[0]]=nz+2*marg;
    ol[axord[0]]=oz-marg*dz;
  }
  if (ndim > 1) { 
    nx=n[axord[1]]; 
    dx=d[axord[1]]; 
    ox=o[axord[1]]; 
    oxg=og[axord[1]];
    nl[axord[1]]=nx+2*marg;
    ol[axord[1]]=ox-marg*dx;
  }
  if (ndim > 2) { 
    ny=n[axord[2]]; 
    dy=d[axord[2]]; 
    oy=o[axord[2]]; 
    oyg=og[axord[2]];
    nl[axord[2]]=ny+2*marg;
    ol[axord[2]]=oy-marg*dy;
  }

  //  fprintf(stderr,"nz=%d nx=%d dz=%g dx=%g oz=%g ox=%g dim=%d\n",nz,nx,dz,dx,oz,ox,ndim);
  /*
  fprintf(stream,"nl = %d %d\n",nl[0],nl[1]);
  fprintf(stream,"ol = %e %e\n",ol[0],ol[1]);
  fprintf(stream,"d  = %e %e\n",d[0],d[1]);
  */

  /* record grid dimension in tg struct */
  tg->ndim = ndim;
  
  /* record axis order for use on output */
  for (i=0;i<RARR_MAX_NDIM;i++) tg->axord[i]=axord[i];

  /* accumulate cell volume */
  tg->dvol=REAL_ONE;
  for (i=0;i<ndim;i++) tg->dvol*=d[axord[i]];

  /* loop over traces in record, reading until either no more traces are
     read or until the source location is no longer the same */

  /*  fprintf(stream,"rk=%d start loop over %d traces\n",retrieveRank(),(tg->ntr)[tg->irec]);*/

  itr=0;          /* number of traces processed from this record */

  tg->ntraces=0;  /* number of traces located in grid */

  /*  if (rk==0) fprintf(stderr,"in init: read loop\n");*/

  while ( (itr < (tg->ntr)[tg->irec]) && (!err) ) {

    /*    fprintf(stream,"rk=%d itr=%d\n",retrieveRank(),itr);*/

    /*    if (rk==0) fprintf(stderr,"in init: call tr_get\n");*/
    //    err=traceserver_get(tg->fpout,
    err=traceserver_get(tg->fpin,
#ifdef IWAVE_USE_MPI
			&(tg->p),
#endif
			&otr); 
    /* any error returned from traceserver functions in MPI case
       is impossible - has resulted in abort */
    if (err) {
      fprintf(stream,"Error: init_tracegeom from traceserver_get\n");
      fflush(stream);
      return err;
    }

    /*    if (rk==0) fprintf(stderr,"in init: process trace\n");*/
  
    /* extract info common TO RECORD from first trace IN RECORD */
    /* bug fix - 28.05.10 WWS - if no traces on processor, ntraces
       never updated and every trace is checked! */
    /*      if (tg->ntraces==0) */
    if (itr==0) {
      /* time step is always that passed as arg */
      tg->dt   = dt;
      
      /* take trace time data from first trace read */
      gethdval(&(otr.tr),"ns",&val);
      tmpnt = vtoi(hdtype("ns"),val);
      
      /* TIME UNITS HARDWIRED TO MS */
      gethdval(&(otr.tr),"dt",&val);
      tmpdt = (1.e-3)*vtof(hdtype("dt"),val);      
      
      /* save int value of header word = dt in musec, for use on output */
      tg->dtmus = vtoi(hdtype("dt"),val);

      gethdval(&(otr.tr),"delrt",&val);
      tmpt0   = vtof(hdtype("delrt"),val);
      
      /* time step and output time data.
	 case 1: user specified output time grid */
      if ( usernt > 0 ) {
	tg->interp=0;
	tg->nt = usernt;
	tg->ntout = tg->nt;
	tg->dtout = tg->dt;
	/* adjust t0 so that 0 is a sample point */
	tg->t0 = dt*((int)floor(usert0/dt));
	tg->t0out = tg->t0;
      }
      /* case 2: output time grid read from file */ 
      else {
	tg->interp=1;
	tg->ntout = tmpnt;
	tg->dtout = tmpdt;
	tg->t0out = tmpt0;
	tg->nt = (int)((tmpnt-1)*tmpdt/dt)+1;
	tg->t0 = dt*((int)floor(tmpt0/dt));
      }
      
      /* in any case, tmax (simulator stopping time) computed
	 from simulation time grid: */
      tg->tmax = tg->t0 + (tg->nt-1)*(tg->dt);	  
      
      /* read scale information */
      gethdval(&(otr.tr),"scalco",&val);
      tg->scalco = vtof(hdtype("scalco"),val);
      gethdval(&(otr.tr),"scalel",&val);
      tg->scalel = vtof(hdtype("scalel"),val);
      
      /* read source position */
      gethdval(&(otr.tr),"sx",&val);
      tmpsx = vtof(hdtype("sx"),val);
      if (tg->scalco > 0) { tmpsx *=  tg->scalco; }
      if (tg->scalco < 0) { tmpsx /= -tg->scalco; }
      
      gethdval(&(otr.tr),"sy",&val);
      tmpsy = vtof(hdtype("sy"),val);
      if (tg->scalco > 0) { tmpsy *=  tg->scalco; }
      if (tg->scalco < 0) { tmpsy /= -tg->scalco; }
      
      gethdval(&(otr.tr),"selev",&val);
      tmpsz = vtof(hdtype("selev"),val);
      if (tg->scalel > 0) { tmpsz *=  tg->scalel; }
      if (tg->scalel < 0) { tmpsz /= -tg->scalel; }
      tmpsz=-tmpsz;
      
      /* store source info, assumed same for all traces
	 25.02.08: change is* to absolute rather than relative int
	 coord 
	 10.07.08: change it back.
	 22.10.08: introduce axis order permutation - only two allowed
	 23.10.08: sensible dimension options
      */
      
      IASN(tg->is,IPNT_0);
      RASN(tg->rs,RPNT_0);
      if ( ndim == 3 ) {
	tg->is[axord[1]]=(int)((tmpsx-oxg)/dx);
	tg->rs[axord[1]]=(tmpsx-oxg-dx*tg->is[axord[1]])/dx;
	tg->is[axord[2]]=(int)((tmpsy-oyg)/dy);
	tg->rs[axord[2]]=(tmpsy-oyg-dy*tg->is[axord[2]])/dy;
	tg->is[axord[0]]=(int)((tmpsz-ozg)/dz);
	tg->rs[axord[0]]=(tmpsz-ozg-dz*tg->is[axord[0]])/dz;
      }
      else if ( ndim == 2 ) {
	tg->is[axord[1]]=(int)((tmpsx-oxg)/dx);
	tg->rs[axord[1]]=(tmpsx-oxg-dx*tg->is[axord[1]])/dx;
	tg->is[axord[2]]=0;
	tg->rs[axord[2]]=0.0;
	tg->is[axord[0]]=(int)((tmpsz-ozg)/dz);
	tg->rs[axord[0]]=(tmpsz-ozg-dz*tg->is[axord[0]])/dz;
      }
      else if ( ndim == 1 ) {
	tg->is[axord[1]]=0;
	tg->rs[axord[1]]=0.0;
	tg->is[axord[2]]=0;
	tg->rs[axord[2]]=0.0;
	tg->is[axord[0]]=(int)((tmpsz-ozg)/dz);
	tg->rs[axord[0]]=(tmpsz-ozg-dz*tg->is[axord[0]])/dz;
      }
      else { 
	fprintf(stream,"Error: tracegeom_init - infeasible ndim=%d\n",ndim);
	fflush(stream);
#ifdef IWAVE_USE_MPI
	MPI_Abort(wcomm,E_BADINPUT);
#else
	return E_BADINPUT; 
#endif 
      }
    }

    tmpz=0.0;
    tmpx=0.0;
    tmpy=0.0;
    
    gethdval(&(otr.tr),"gx",&val);
    tmpx = vtof(hdtype("gx"),val);
    if (tg->scalco > 0) { tmpx *=  tg->scalco;}
    if (tg->scalco < 0) { tmpx /= -tg->scalco;}
    
    gethdval(&(otr.tr),"gy",&val);
    tmpy = vtof(hdtype("gy"),val);
    if (tg->scalco > 0) { tmpy *=  tg->scalco;}
    if (tg->scalco < 0) { tmpy /= -tg->scalco;}
    
    gethdval(&(otr.tr),"gelev",&val);
    tmpz = vtof(hdtype("gelev"),val);
    if (tg->scalel > 0) { tmpz *=  tg->scalel; }
    if (tg->scalel < 0) { tmpz /= -tg->scalel; }
    tmpz =- tmpz; 
    
    /* For all subsequent traces read, check that traces 
       are all same length */
    gethdval(&(otr.tr),"ns",&val);
    tmpnt=vtoi(hdtype("ns"),val);
    /* if usernt not overriding trace nt, make sure all traces agree
       otherwise who cares */
    if (!usernt && tmpnt != tg->ntout) {
      fprintf(stream,"Error: tracegeom_init - change in nt from %d to %d\n",
	      tg->ntout,tmpnt);
      fflush(stream);
#ifdef IWAVE_USE_MPI
      MPI_Abort(wcomm,E_OTHER);
#else
      return E_OTHER;
#endif
    }
    
    /* store trace-dep info for each trace properly inside the grid */
    
    IASN(tg->ig[tg->ntraces],IPNT_0);
    RASN(tg->rg[tg->ntraces],RPNT_0);
    x[axord[0]]=tmpz;
    x[axord[1]]=tmpx;
    x[axord[2]]=tmpy;
    
    init[itr]=0;
    
    if ( (ndim ==3) &&	
	 ringrid(ndim,nl,d,ol,x) ) {
      tg->ig[tg->ntraces][axord[1]]=(int)((tmpx-oxg)/dx);
      tg->rg[tg->ntraces][axord[1]]=(tmpx-oxg-dx*tg->ig[tg->ntraces][axord[1]])/dx;
      tg->ig[tg->ntraces][axord[2]]=(int)((tmpy-oyg)/dy);
      tg->rg[tg->ntraces][axord[2]]=(tmpy-oyg-dy*tg->ig[tg->ntraces][axord[2]])/dy;
      tg->ig[tg->ntraces][axord[0]]=(int)((tmpz-ozg)/dz);
      tg->rg[tg->ntraces][axord[0]]=(tmpz-ozg-dz*tg->ig[tg->ntraces][axord[0]])/dz;
      init[itr] = 1;
    }
    if ( (ndim == 2) &&  
	 ringrid(ndim,nl,d,ol,x) ) {
      (tg->ig)[tg->ntraces][axord[1]]=(int)((tmpx-oxg)/dx);
      (tg->rg)[tg->ntraces][axord[1]]=(tmpx-oxg-dx*(tg->ig)[tg->ntraces][axord[1]])/dx;
      (tg->ig)[tg->ntraces][axord[2]]=0;
      (tg->rg)[tg->ntraces][axord[2]]=0.0;
      (tg->ig)[tg->ntraces][axord[0]]=(int)((tmpz-ozg)/dz);
      (tg->rg)[tg->ntraces][axord[0]]=(tmpz-ozg-dz*(tg->ig)[tg->ntraces][axord[0]])/dz;
      init[itr] = 1;
    }
    if ( (ndim == 1) &&
	 ringrid(ndim,nl,d,ol,x) ) {
      (tg->ig)[tg->ntraces][axord[1]]=0;
      (tg->rg)[tg->ntraces][axord[1]]=0.0;
      (tg->ig)[tg->ntraces][axord[2]]=0;
      (tg->rg)[tg->ntraces][axord[2]]=0.0;
      (tg->ig)[tg->ntraces][axord[0]]=(int)((tmpz-ozg)/dz);
      (tg->rg)[tg->ntraces][axord[0]]=(tmpz-ozg-dz*tg->ig[tg->ntraces][axord[0]])/dz;
      init[itr] = 1;
    }
    
    if (init[itr]) {
      gethdval(&(otr.tr),"tracl",&val);
      tg->tracl[tg->ntraces]=vtoi(hdtype("tracl"),val);
      gethdval(&(otr.tr),"tracr",&val);
      tg->tracr[tg->ntraces]=vtoi(hdtype("tracr"),val);
      gethdval(&(otr.tr),"fldr",&val);
      tg->fldr[tg->ntraces]=vtoi(hdtype("fldr"),val);
      gethdval(&(otr.tr),"tracf",&val);
      tg->tracf[tg->ntraces]=vtoi(hdtype("tracf"),val);
      // WWS 30.11.10
      tg->troff[tg->ntraces]=otr.m;
      //      tg->troff[tg->ntraces]=tg->recoff[tg->irec]+itr*(HDRBYTES + tg->ntout * sizeof(float));
      
      /*	fprintf(stream,"rk=%d ntraces=%d offset=%ld\n",retrieveRank(),tg->ntraces,tg->troff[tg->ntraces]);*/
      tg->ntraces++;
      if (tg->ntraces > MAX_TRACES) {
	fprintf(stream,"Error: tracegeom_init - ntraces exceeds MAX_TRACES = %d\n",MAX_TRACES);
	fflush(stream);
#ifdef IWAVE_USE_MPI
	MPI_Abort(wcomm,E_ALLOC);
#else
	return E_ALLOC;
#endif	
      }
    }      
    /* one more trace read... */
    itr++; /* in record */
    
  }
   
  /* allocate data buffer, optionally store data */
#ifdef IWAVE_VERBOSE
  fprintf(stream,"NOTE: allocating %d traces of %d samples in tracegeom buffer\n",tg->ntraces,tg->nt);
  fflush(stream);
#endif

  /*  if (tg->ntraces) { */
  /* fix of 26.07.10: since this is SPMD code, have to call
     traceserver_get (which contains Bcasts) whether there are any
     traces or not */
  /* fix of 20.10.10: HUGE memory leak! Resolution: first time, buf=NULL
     since follows construct call. Subsequent shots - buf non-null and
     should be freed.
  */
  if (tg->buf) {
    userfree_(tg->buf);
    tg->buf = NULL;
  }
  if (tg->ntraces) tg->buf =(float *)usermalloc_(tg->ntraces * tg->nt * sizeof(float));
  if (tg->ntraces && !tg->buf) {
    fprintf(stream,"Error: tracegeom_init - failed to allocate trace\n");
    fprintf(stream,"buffer of length %ld\n", 
	    (long) (tg->ntraces*tg->nt*sizeof(float)));
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_ALLOC);
#else
    return E_ALLOC;
#endif
  }
  if (tg->ntraces) {
    for (i=0; i < tg->ntraces * tg->nt; i++) tg->buf[i] = 0.0;
  }
    
  /* if data is to be initialized, re-read file, adjoint-interpolate
     active traces onto simulation time grid.
  */
  if (initbuf) {
    if (traceserver_seek(tg->fpin,&(tg->recoff[tg->irec]))) return E_FILE;
    iinit=0;
    wlen=cubicadj_getworksize(tg->nt,tmpnt);
    work=(float *)usermalloc_(wlen*sizeof(float));
    if (!work) {
      fprintf(stream,"Error: tracegeom_init - failed to allocate work\n");
      fprintf(stream,"buffer of length %ld\n",(long) (wlen*sizeof(float)));
      fflush(stream);
#ifdef IWAVE_USE_MPI
      MPI_Abort(wcomm,E_ALLOC);
#else
      return E_ALLOC;
#endif 
    }

    /*  fprintf(stream,"loop over %d traces\n",tg->ntraces); fflush(stream);*/
    for (i=0;i<(tg->ntr)[tg->irec];i++) {
      /*      fprintf(stream,"get tr=%d\n"); fflush(stream);*/
      err=traceserver_get(tg->fpin,
#ifdef IWAVE_USE_MPI
			  &(tg->p),
#endif
			  &otr); 
      if (err) {
	fprintf(stream,"init_tracegeom - error from traceserver_get\n");
	fflush(stream);
#ifdef IWAVE_USE_MPI
	MPI_Abort(wcomm,err);
#else
	return err;
#endif
      }
      /*      fprintf(stream,"interp tr=%d\n"); fflush(stream);*/
      /* init[i] is only set if trace i recvr location lies in the domain,
	 in particular if ntraces=0 for this domain then cubicadj is 
	 never called */
      if (init[i]) {
	err=cubicadj_(&tmpt0,   &tmpdt,   (otr.tr).data,             &tmpnt,
		      &(tg->t0),&(tg->dt),&((tg->buf)[iinit*tg->nt]),&(tg->nt),
		      &iend,    work,     &wlen);
	if (err) {
	  fprintf(stream,"Error: tracegeom_init from cubicadj, err=%d\n",
		  err);
	  fflush(stream);
#ifdef IWAVE_USE_MPI
	  MPI_Abort(wcomm,err);
#else
	  return err;
#endif
	}
	iinit++;
      }
      /*      fprintf(stream,"iinit=%d\n",iinit); fflush(stream);*/
    }
    userfree_(work);
    /*    fprintf(stream,"init_tracegeom: exit load loop\n"); fflush(stream);*/
  }	       

  /*  fprintf(stream,"init_tracegeom: return with err=%d\n",err); fflush(stream);*/

  userfree_(init);

  return err;
}

void destroy_tracegeom(tracegeom * tg) {
  tg->nrec=0;
  tg->irec=0;
  tg->xrec=0;
  tg->first=0;
  tg->last=0;
  fflush(tg->fpout);
#ifdef IWAVE_USE_FMGR
  if (tg->fpin) iwave_fclose(tg->fpin);
  if (tg->fpout) iwave_fclose(tg->fpout);
#else 
  if (tg->fpin) fclose(tg->fpin);
  if (tg->fpout) fclose(tg->fpout);
#endif
  tg->fpin=NULL;
  tg->fpout=NULL;
  /* comment out setnull_tracegeom() (D.S. 06.10.10) */
  /* setnull_tracegeom(tg); */

  /* moved from setnull, which is a default constructor */
  if (tg->buf) userfree_(tg->buf); tg->buf=NULL;
  if (tg->ntr) userfree_(tg->ntr); tg->ntr=NULL;
  if (tg->recoff) userfree_(tg->recoff); tg->recoff=NULL;
  if (tg->src) userfree_(tg->src); tg->src=NULL;
  if (tg->ig) userfree_(tg->ig); tg->ig=NULL;
  if (tg->rg) userfree_(tg->rg); tg->rg=NULL;
  if (tg->tracl) userfree_(tg->tracl); tg->tracl=NULL;
  if (tg->tracr) userfree_(tg->tracr); tg->tracr=NULL;
  if (tg->tracf) userfree_(tg->tracf); tg->tracf=NULL;
  if (tg->fldr) userfree_(tg->fldr); tg->fldr=NULL;
  if (tg->troff) userfree_(tg->troff); tg->troff=NULL;

#ifdef IWAVE_USE_MPI
  MPI_Type_free(&(tg->p));
#endif
}

void setnull_tracegeom(tracegeom * tg) {
  /* nullify file pointers - permits all procs other than
     rk=0 to avoid file ops by testing
  */
  tg->fpin=NULL;
  tg->fpout=NULL;

  /* Since struct members are not initialized when the struct
     is allocated, it's not possible to assume that buf pointing
     to anything other than NULL is a sign of previous allocation.
     Rewrite: make sure that setnull is called ONLY ONCE as a 
     de facto default constructor. Never after construction.
     If this pattern is followed, then the destructor (above) can
     safely free buf.
  fprintf(stderr,"TRACEIO::SETNULL_TRACEGEOM -> free\n");
  if (tg->buf) free(tg->buf); 
  fprintf(stderr,"TRACEIO::SETNULL_TRACEGEOM -> freed\n");
  */
  tg->buf=NULL;
  /* added in move do dynamic mem 15.09.10 */
  tg->ig=NULL;
  tg->rg=NULL;
  tg->tracl=NULL;
  tg->tracr=NULL;
  tg->tracf=NULL;
  tg->fldr=NULL;
  tg->troff=NULL;
  /* end addition */
  tg->ntraces=0;
  tg->nt=0;
  tg->ntout=0;
  tg->dt=0.0;
  tg->dtout=0.0;
  tg->t0=0.0;
  tg->t0out=0.0;
  tg->scalel=0;
  tg->scalco=0;
  tg->ndim=0;
  /* added by D.S. 06.10.10 */
  tg->nrec=0;
  tg->xrec=0;
  tg->irec=0;
  tg->first=0;
  tg->last=0;
  tg->dtmus=0.0;
  tg->interp=0;
  tg->dvol=0.0;
  tg->ntr=NULL;
  tg->recoff=NULL;
  tg->src=NULL;
  /* end addition */
  /*  fprintf(stderr,"TRACEIO::SETNULL_TRACEGEOM -> exit\n");*/

}

void fprint_tracegeom(tracegeom const * tg,FILE *fp) {
  int i;
  fprintf(fp,"TRACE GEOMETRY STRUCT\n");
  fprintf(fp,"  -- global quantities -- \n");
  fprintf(fp,"nrec      = %12.1d\n",tg->nrec);
  fprintf(fp,"irec      = %12.1d\n",tg->irec);
  fprintf(fp,"xrec      = %12.1d\n",tg->xrec);
  fprintf(fp,"first     = %12.1d\n",tg->first);
  fprintf(fp,"last      = %12.1d\n",tg->last);
  fprintf(fp,"recoff    = %12.1ld\n",tg->recoff[tg->irec]);
  fprintf(fp,"ndim      = %12.1d\n",tg->ndim);
  fprintf(fp,"dvol      = %12.4e\n",tg->dvol);
  if (tg->ndim > 0) fprintf(fp,"sz        = %12.4e\n",tg->src[tg->irec][0]);
  if (tg->ndim > 1) fprintf(fp,"sx        = %12.4e\n",tg->src[tg->irec][1]);
  if (tg->ndim > 2) fprintf(fp,"sy        = %12.4e\n",tg->src[tg->irec][2]);
  fprintf(fp,"  -- record quantities -- \n");
  if (tg->nrec>0) {
    /*fprintf(fp,"ntr[irec] = %12.1d\n",(tg->ntr)[tg->irec]);*/
    for(i=0; i < tg->nrec; i++)
      fprintf(fp,"ntr[%d] = %12.1d\n",i,(tg->ntr)[i]);
  }
  fprintf(fp,"ntraces   = %12.1d\n",tg->ntraces);
  fprintf(fp,"scalel    = %12.1d\n",tg->scalel);
  fprintf(fp,"scalco    = %12.1d\n",tg->scalco);
  if (tg->ndim > 0) fprintf(fp,"isz       = %12.1d\n",tg->is[tg->axord[0]]);
  if (tg->ndim > 1) fprintf(fp,"isx       = %12.1d\n",tg->is[tg->axord[1]]);
  if (tg->ndim > 2) fprintf(fp,"isy       = %12.1d\n",tg->is[tg->axord[2]]);
  if (tg->ndim > 0) fprintf(fp,"rsz       = %12.4e\n",tg->rs[tg->axord[0]]);
  if (tg->ndim > 1) fprintf(fp,"rsx       = %12.4e\n",tg->rs[tg->axord[1]]);
  if (tg->ndim > 2) fprintf(fp,"rsy       = %12.4e\n",tg->rs[tg->axord[2]]);
  fprintf(fp,"nt        = %12.1d\n",tg->nt);
  fprintf(fp,"ntout     = %12.1d\n",tg->ntout);
  fprintf(fp,"dt        = %12.4e\n",tg->dt);
  fprintf(fp,"dtout     = %12.4e\n",tg->dtout);
  fprintf(fp,"dt (musec)= %12.1d\n",tg->dtmus);
  fprintf(fp,"t0        = %12.4e\n",tg->t0);
  fprintf(fp,"t0out     = %12.4e\n",tg->t0out);
  fprintf(fp,"interp    = %12.1d\n",tg->interp);
  fprintf(fp,"    itr");
  if (tg->ndim > 0) fprintf(fp,"    igz");
  if (tg->ndim > 1) fprintf(fp,"    igx");
  if (tg->ndim > 2) fprintf(fp,"    igy");
  if (tg->ndim > 0) fprintf(fp,"      rgz    ");
  if (tg->ndim > 1) fprintf(fp,"      rgx    ");
  if (tg->ndim > 2) fprintf(fp,"      rgy    ");
  if (tg->ndim > 0) fprintf(fp,"  tracr");
  if (tg->ndim > 0) fprintf(fp,"   offset");
  fprintf(fp,"\n");
  for (i=0;i<tg->ntraces;i++) {
    fprintf(fp,"%6.d",i);
    if (tg->ndim > 0) fprintf(fp," %6.1d",tg->ig[i][tg->axord[0]]);
    if (tg->ndim > 1) fprintf(fp," %6.1d",tg->ig[i][tg->axord[1]]);
    if (tg->ndim > 2) fprintf(fp," %6.1d",tg->ig[i][tg->axord[2]]);
    if (tg->ndim > 0) fprintf(fp,"   %10.4e",tg->rg[i][tg->axord[0]]);
    if (tg->ndim > 1) fprintf(fp,"   %10.4e",tg->rg[i][tg->axord[1]]);
    if (tg->ndim > 2) fprintf(fp,"   %10.4e",tg->rg[i][tg->axord[2]]);
    if (tg->ndim > 0) fprintf(fp," %6.1d",tg->tracr[i]);
    if (tg->ndim > 0) fprintf(fp," %10.1ld",tg->troff[i]);
    fprintf(fp,"\n");
  }
  fprintf(fp, "\n");
  fflush(fp);
}

void print_tracegeom(tracegeom const * tg) {
  fprint_tracegeom(tg,stdout);
}

void sampletraces(tracegeom * tg,
		  int order,
		  int load,
		  int it, 
		  IPNT n0,
		  IPNT gs0,
		  IPNT n,
		  IPNT gs,
		  ireal * d,
		  //		  IPNT n0m,
		  //		  IPNT gs0m,
		  //		  IPNT nm,
		  //		  IPNT gsm,		  
		  //		  ireal * m) {
		  ireal mult) {
  /*****************************************
   * gs input added 25.02.08: must treat 
   * ix,iy,iz as GLOBAL grid indices!!! 
   *****************************************
   * only orders 0,1 implemented: 06.21.08 
   *****************************************
   * input pointer, avoid copy - 15.01.09 
   *****************************************
   * test for being in grid - let accum. 
   * in writetraces handle
   * sampling by multiple domains. use 
   * ingrid function to determine in/out.
   * 11.03.09
   *****************************************
   * observe that test for grid should refer 
   * to COMPUTATIONAL grid, whereas offset
   * has to be computed for ALLOCATED grid
   * so BOTH have to be passed. This would
   * be simpler if it were an rarray app, as
   * then could use rget.
   * 14.3.09
   *****************************************
   * add load option 08.12.09
   *****************************************
   * add multiplier field option 25.01.10
   *****************************************
   * add global grid info for mult field -
   * cannot assume that size, shape are same
   * so must do independent offset comp.
   * 26.01.10
   *****************************************
   * remove multiplier field option - corr. 
   * to change in adjoint state approach - 
   * go back to "pure" load and save, which
   * are adjoints of each other, but with
   * on-the-fly scale factor
   *****************************************/

  int itr;    /* trace counter */
  int ioff=0;   /* offset into sampled array */
  //  int moff;   /* offset into multiplier array */
  int ndim;   /* problem dimension */
  IPNT ind;   /* integer part of sample coords */
  
  /* convenience storage for multiplier */
  //  ireal mult;
  ireal fac;

  if (it<0 || it>tg->nt-1) return;
  
  ndim=tg->ndim;
  IASN(ind,IPNT_0);
  //  /* default value for multiplier, just to make code simple */
  // this is now an input - 03.12
  //  mult=REAL_ONE;
  /* generic scale factor - accounts for weights on L2 norms for
     viewing load option as adjoint of save
  */
  // big change 03.12: scale by mult; also by ratio of cell vols for adjoint
  if (load) fac=mult*tg->dt/tg->dvol;
  else fac=mult;

  for (itr=0;itr<tg->ntraces;itr++) {
    
    /* set up index tuple, grid offset */
    if (ndim > 0) {
      ind[0]=(tg->ig)[itr][0];
      ioff = ind[0]-gs0[0];
    }
    if (ndim > 1) {
      ind[1]=(tg->ig)[itr][1];
      ioff+= (ind[1]-gs0[1] )*n0[0];
    }
    if (ndim > 2) {
      ind[2]=(tg->ig)[itr][2];
      ioff+= (ind[2]-gs0[2] )*n0[0] *n0[1];      
    }

    if (load) {

      // printf("sampletraces: trace=%d ioff=%d moff=%d ind[0]=%d ind[1]=%d\n gsm[0]=%d nm[0]=%d gsm[1]=%d nm[1]=%d gs0m[0]=%d n0m[0]=%d gs0m[1]=%d n0m[1]=%d\n",itr,ioff,moff,ind[0],ind[1],gsm[0],nm[0],gsm[1],nm[1],gs0m[0],n0m[0],gs0m[1],n0m[1]);
       
      if (order==0) {
	if (ingrid(ndim,n,gs,ind)) {
	  //	  mult=fac;
	  //	  if (m && ingrid(ndim,nm,gsm,ind)) mult*=m[moff];
	  d[ioff]+=fac*(tg->buf)[it+itr*tg->nt];
	}
      }

      else if (order==1) {
	if (ndim==1) {
	  /* 0,0,0 */
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    moff=ind[0]-gs0m[0];
	    //	    mult=fac;
	    //	    if (m && ingrid(ndim,nm,gsm,ind)) mult*=m[moff];
	    d[ioff]+=fac*
	      (1.0-(tg->rg)[itr][0])*
	      (tg->buf)[it+itr*tg->nt];
	  }
	  /* 1,0,0 */
	  ind[0]++;
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    moff=ind[0]-gs0m[0];
	    //	    mult=fac;
	    //	    if (m && ingrid(ndim,nm,gsm,ind)) mult*=m[moff+1];
	    d[ioff+1]+=fac*
	      (tg->rg)[itr][0]*
	      (tg->buf)[it+itr*tg->nt];
	  }
	  ind[0]--;
	}
	else if (ndim == 2) {
	  /* 0,0,0 */
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    moff=ind[0]-gs0m[0] + (ind[1]-gs0m[1])*n0m[0];
	    //	    mult=fac;
	    //	    if (m && ingrid(ndim,nm,gsm,ind)) mult*=m[moff];	    
	    d[ioff]+=fac*
	      (1.0-(tg->rg)[itr][0])*
	      (1.0-(tg->rg)[itr][1])*
	      (tg->buf)[it+itr*tg->nt];
	    //	    printf("sampletraces: itr=%d ig=[0,0] ind=[%d,%d] moff=%d m=%20.14e\n",itr,ind[0],ind[1],moff,m[moff]);
	  }
	  /* 1,0,0 */
	  ind[0]++;
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    moff=ind[0]-gs0m[0] + (ind[1]-gs0m[1])*n0m[0];
	    //	    mult=fac;
	    //	    if (m && ingrid(ndim,nm,gsm,ind)) mult*=m[moff];
	    d[ioff+1]+=fac*
	      (tg->rg)[itr][0]*
	      (1.0-(tg->rg)[itr][1])*
	      (tg->buf)[it+itr*tg->nt];
	    //	    printf("sampletraces: itr=%d ig=[1,0] ind=[%d,%d] moff=%d m=%20.14e\n",itr,ind[0],ind[1],moff,m[moff]);
	  }
	  ind[0]--;
	  /* 0,1,0 */
	  ind[1]++;
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    mult=fac;
	    // 	    moff=ind[0]-gs0m[0] + (ind[1]-gs0m[1])*n0m[0];
	    //	    if (m && ingrid(ndim,nm,gsm,ind)) mult*=m[moff]; 
	    d[ioff+n0[0]]+=fac*
	      (1.0-(tg->rg)[itr][0])*
	      (tg->rg)[itr][1]*
	      (tg->buf)[it+itr*tg->nt];
	    //	    printf("sampletraces: itr=%d ig=[0,1] ind=[%d,%d] moff=%d m=%20.14e\n",itr,ind[0],ind[1],moff,m[moff]);
	  }
	  ind[1]--;
	  /* 1,1,0 */
	  ind[0]++;
	  ind[1]++;
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    mult=fac;
	    // 	    moff=ind[0]-gs0m[0] + (ind[1]-gs0m[1])*n0m[0];
	    //	    if (m && ingrid(ndim,nm,gsm,ind)) mult*=m[moff+1+n0[0]]; 
	    d[ioff+1+n0[0]]+=fac*
	      (tg->rg)[itr][0]*
	      (tg->rg)[itr][1]*
	      (tg->buf)[it+itr*tg->nt];
	    //	    printf("sampletraces: itr=%d ig=[1,1] ind=[%d,%d] moff=%d m=%20.14e\n",itr,ind[0],ind[1],moff,m[moff]);
	  }
	  ind[0]--;
	  ind[1]--;
	}
	else if (ndim == 3) {
	  /* 0,0,0 */
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    mult=fac;
	    //	    moff=ind[0]-gs0m[0] + (ind[1]-gs0m[1])*n0m[0] + (ind[2]-gs0m[2])*n0m[0]*n0m[1];
	    //	    if (m && ingrid(ndim,nm,gsm,ind)) 
	    //	      mult*=m[moff]; 
	    d[ioff]+=fac*
	      (1.0-(tg->rg)[itr][0])*
	      (1.0-(tg->rg)[itr][1])*
	      (1.0-(tg->rg)[itr][2])*
	      (tg->buf)[it+itr*tg->nt];
	  }
	  /* 1,0,0 */
	  ind[0]++;
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    mult=fac;
	    //	    moff=ind[0]-gs0m[0] + (ind[1]-gs0m[1])*n0m[0] + (ind[2]-gs0m[2])*n0m[0]*n0m[1];
	    //	    if (m && ingrid(ndim,nm,gsm,ind)) 
	    //	      mult*=m[moff+1]; 
	    d[ioff+1]+=fac*
	      (tg->rg)[itr][0]*
	      (1.0-(tg->rg)[itr][1])*
	      (1.0-(tg->rg)[itr][2])*
	      (tg->buf)[it+itr*tg->nt];
	  }
	  ind[0]--;
	  /* 0,1,0 */
	  ind[1]++;
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    mult=fac;
	    //	    moff=ind[0]-gs0m[0] + (ind[1]-gs0m[1])*n0m[0] + (ind[2]-gs0m[2])*n0m[0]*n0m[1];
	    //	    if (m && ingrid(ndim,nm,gsm,ind)) 
	    //	      mult*=m[moff+n0[0]];
	    d[ioff+n0[0]]+=fac*
	      (1.0-(tg->rg)[itr][0])*
	      (tg->rg)[itr][1]*
	      (1.0-(tg->rg)[itr][2])*
	      (tg->buf)[it+itr*tg->nt];
	  }
	  ind[1]--;
	  /* 0,0,1 */
	  ind[2]++;
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    mult=fac;
	    //	    moff=ind[0]-gs0m[0] + (ind[1]-gs0m[1])*n0m[0] + (ind[2]-gs0m[2])*n0m[0]*n0m[1];
	    //	    if (m && ingrid(ndim,nm,gsm,ind)) 
	    //	      mult*=m[moff+n0[0]*n0[1]];
	    d[ioff+n0[0]*n0[1]]+=fac*
	      (1.0-(tg->rg)[itr][0])*
	      (1.0-(tg->rg)[itr][1])*
	      (tg->rg)[itr][2]*
	      (tg->buf)[it+itr*tg->nt];
	  }
	  ind[2]--;
	  /* 1,1,0 */
	  ind[0]++;
	  ind[1]++;
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    mult=fac;
	    //	    moff=ind[0]-gs0m[0] + (ind[1]-gs0m[1])*n0m[0] + (ind[2]-gs0m[2])*n0m[0]*n0m[1];
	    //	    if (m && ingrid(ndim,nm,gsm,ind)) 
	    //	      mult*=m[moff+1+n0[0]];
	    d[ioff+1+n0[0]]+=fac*
	      (tg->rg)[itr][0]*
	      (tg->rg)[itr][1]*
	      (1.0-(tg->rg)[itr][2])*
	      (tg->buf)[it+itr*tg->nt];
	  }
	  ind[0]--;
	  ind[1]--;
	  /* 1,0,1 */
	  ind[0]++;
	  ind[2]++;
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    mult=fac;
	    //	    moff=ind[0]-gs0m[0] + (ind[1]-gs0m[1])*n0m[0] + (ind[2]-gs0m[2])*n0m[0]*n0m[1];
	    //	    if (m && ingrid(ndim,nm,gsm,ind)) 
	    //	      mult*=m[moff+1+n0[0]*n0[1]]; 
	    d[ioff+1+n0[0]*n0[1]]+=fac*
	      (tg->rg)[itr][0]*
	      (1.0-(tg->rg)[itr][1])*
	      (tg->rg)[itr][2]*
	      (tg->buf)[it+itr*tg->nt];
	  }
	  ind[0]--;
	  ind[2]--;
	  /* 0,1,1 */
	  ind[1]++;
	  ind[2]++;
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    mult=fac;
	    //	    moff=ind[0]-gs0m[0] + (ind[1]-gs0m[1])*n0m[0] + (ind[2]-gs0m[2])*n0m[0]*n0m[1];
	    //	    if (m && ingrid(ndim,nm,gsm,ind)) 
	    //	      mult*=m[moff+n[0]+n0[0]*n0[1]];
	    d[ioff+n0[0]+n0[0]*n0[1]]+=fac*
	      (1.0-(tg->rg)[itr][0])*
	      (tg->rg)[itr][1]*
	      (tg->rg)[itr][2]*
	      (tg->buf)[it+itr*tg->nt];
	  }
	  ind[1]--;
	  ind[2]--;
	  /* 1,1,1 */
	  ind[0]++;
	  ind[1]++;
	  ind[2]++;
	  if (ingrid(ndim,n,gs,ind)) {
	    //	    mult=fac;
	    //	    moff=ind[0]-gs0m[0] + (ind[1]-gs0m[1])*n0m[0] + (ind[2]-gs0m[2])*n0m[0]*n0m[1];
	    //	    if (m && ingrid(ndim,nm,gsm,ind))
	    //	      mult*=m[moff+1+n0[0]+n0[0]*n0[1]];	    
	    d[ioff+1+n0[0]+n0[0]*n0[1]]+=fac*
	      (tg->rg)[itr][0]*
	      (tg->rg)[itr][1]*
	      (tg->rg)[itr][2]*
	      (tg->buf)[it+itr*tg->nt];
	  }
	  ind[0]--;
	  ind[1]--;
	  ind[2]--;
	}
	else {
	  return;
	}
      }
      else {
	return;
      }
    }
    else {
      
      /* this should be done ahead of time!!! */
      /* (tg->buf)[it+itr*tg->nt]=REAL_ZERO;*/
      
      if (order==0) {
	if (ingrid(ndim,n,gs,ind)) 
	  (tg->buf)[it+itr*tg->nt]+=fac*d[ioff];
/*
	fprintf(stderr,"sample=%e buffer index=%d\n",(tg->buf)[it+itr*tg->nt],it+itr*tg->nt);
*/
      }

      else if (order==1) {
	if (ndim==1) {
	  /* 0,0,0 */
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+=
	      (1.0-(tg->rg)[itr][0])*fac*d[ioff];
	  /* 1,0,0 */
	  ind[0]++;
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+= (tg->rg)[itr][0]*fac*d[ioff+1];
	  ind[0]--;
	}
	else if (ndim == 2) {
	  /* 0,0,0 */
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+=
	      (1.0-(tg->rg)[itr][0])*(1.0-(tg->rg)[itr][1])*fac*d[ioff];
	  /* 1,0,0 */
	  ind[0]++;
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+=
	      (tg->rg)[itr][0]*(1.0-(tg->rg)[itr][1])*fac*d[ioff+1];
	  ind[0]--;
	  /* 0,1,0 */
	  ind[1]++;
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+=
	      (1.0-(tg->rg)[itr][0])*(tg->rg)[itr][1]*fac*d[ioff+n0[0]];
	  ind[1]--;
	  /* 1,1,0 */
	  ind[0]++;
	  ind[1]++;
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+=
	      (tg->rg)[itr][0]*(tg->rg)[itr][1]*fac*d[ioff+1+n0[0]];
	  ind[0]--;
	  ind[1]--;
	}
	else if (ndim == 3) {
	  /* 0,0,0 */
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+=
	      (1.0-(tg->rg)[itr][0])*(1.0-(tg->rg)[itr][1])*(1.0-(tg->rg)[itr][2])*fac*d[ioff];
	  /* 1,0,0 */
	  ind[0]++;
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+=
	      (tg->rg)[itr][0]*(1.0-(tg->rg)[itr][1])*(1.0-(tg->rg)[itr][2])*fac*d[ioff+1];
	  ind[0]--;
	  /* 0,1,0 */
	  ind[1]++;
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+=
	      (1.0-(tg->rg)[itr][0])*(tg->rg)[itr][1]*(1.0-(tg->rg)[itr][2])*fac*d[ioff+n0[0]];
	  ind[1]--;
	  /* 0,0,1 */
	  ind[2]++;
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+=
	      (1.0-(tg->rg)[itr][0])*(1.0-(tg->rg)[itr][1])*(tg->rg)[itr][2]*fac*d[ioff+n0[0]*n0[1]];
	  ind[2]--;
	  /* 1,1,0 */
	  ind[0]++;
	  ind[1]++;
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+=
	      (tg->rg)[itr][0]*(tg->rg)[itr][1]*(1.0-(tg->rg)[itr][2])*fac*d[ioff+1+n0[0]];
	  ind[0]--;
	  ind[1]--;
	  /* 1,0,1 */
	  ind[0]++;
	  ind[2]++;
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+=
	      (tg->rg)[itr][0]*(1.0-(tg->rg)[itr][1])*(tg->rg)[itr][2]*fac*d[ioff+1+n0[0]*n0[1]];
	  ind[0]--;
	  ind[2]--;
	  /* 0,1,1 */
	  ind[1]++;
	  ind[2]++;
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+=
	      (1.0-(tg->rg)[itr][0])*(tg->rg)[itr][1]*(tg->rg)[itr][2]*fac*d[ioff+n0[0]+n0[0]*n0[1]];
	  ind[1]--;
	  ind[2]--;
	  /* 1,1,1 */
	  ind[0]++;
	  ind[1]++;
	  ind[2]++;
	  if (ingrid(ndim,n,gs,ind)) 
	    (tg->buf)[it+itr*tg->nt]+=
	      (tg->rg)[itr][0]*(tg->rg)[itr][1]*(tg->rg)[itr][2]*fac*d[ioff+1+n0[0]+n0[0]*n0[1]];
	  ind[0]--;
	  ind[1]--;
	  ind[2]--;
	}
	else {
	  return;
	}
      }
      else {
	return;
      }
    }
  }
}

/* 17.01.09: abstracted from writetraces */
int assembletrace(tracegeom const * tg,
		  offsegy * otr,
		  int nb,
		  float dz,
		  float dx,
		  float dy,
		  float ozg,
		  float oxg,
		  float oyg,
		  int idz,
		  int idx,
		  int idy,
		  float * work,
		  int wl) {

  Value val;
  
  float tmp;
  float tmps;
  
  int err=0;
  int iend=1;

#ifdef IWAVE_USE_MPI
  MPI_Comm wcomm=retrieveGlobalComm();  /* global communicator */
#endif

  if ((0 > nb) || (nb>tg->ntraces-1)) {
    fprintf(stderr,"ERROR: [tracegeom::assembletraces]\n");
    fprintf(stderr,"trace index %d out of range [0,%d\n]\n",nb,tg->ntraces-1);
    err=E_OTHER;
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,err);
#else
    return err;
#endif
  }
 
  /* clear out header */
  memset(&(otr->tr),0,240);

  val.h=tg->scalco;
  puthdval(&(otr->tr),"scalco",&val);
  
  tmp =oxg + ((tg->ig)[nb][idx]+(tg->rg)[nb][idx])*dx;
  tmps=oxg +((tg->is)[idx]+(tg->rs)[idx])*dx;
  if (tg->scalco > 0) { tmp /=  tg->scalco; tmps /=  tg->scalco; }
  if (tg->scalco < 0) { tmp *= -tg->scalco; tmps *= -tg->scalco; }
  val.i=tmp;
  puthdval(&(otr->tr),"gx",&val);
  val.i=tmps;
  puthdval(&(otr->tr),"sx",&val);
  
  tmp =oyg + ((tg->ig)[nb][idy]+(tg->rg)[nb][idy])*dy;
  tmps=oyg + ((tg->is)[idy]+(tg->rs)[idy])*dy;
  if (tg->scalco > 0) { tmp /=  tg->scalco; tmps /=  tg->scalco; }
  if (tg->scalco < 0) { tmp *= -tg->scalco; tmps *= -tg->scalco; }
  val.i=tmp;
  puthdval(&(otr->tr),"gy",&val);
  val.i=tmps;
  puthdval(&(otr->tr),"sy",&val);
  
  val.h=tg->scalel;
  puthdval(&(otr->tr),"scalel",&val);
  
  tmp = -(ozg + ((tg->ig)[nb][idz] + (tg->rg)[nb][idz])*dz);
  tmps= -(ozg + ((tg->is)[idz]+(tg->rs)[idz])*dz);
  if (tg->scalel > 0) { tmp /=  tg->scalel; tmps /=  tg->scalel; }
  if (tg->scalel < 0) { tmp *= -tg->scalel; tmps *= -tg->scalel; }
  val.i=tmp;
  puthdval(&(otr->tr),"gelev",&val);
  val.i=tmps;
  puthdval(&(otr->tr),"selev",&val);
  
  /*  val.u=1000.0*(tg->dtout);
   01.10.10 - have stored int value in musec, use it */
  val.u=tg->dtmus;
  puthdval(&(otr->tr),"dt",&val);
  
  val.u=tg->ntout;
  puthdval(&(otr->tr),"ns",&val);
  
  val.h=tg->t0out;
  puthdval(&(otr->tr),"delrt",&val);
  
  val.i=(tg->tracl)[nb];
  puthdval(&(otr->tr),"tracl",&val);
  
  val.i=(tg->tracr)[nb];
  puthdval(&(otr->tr),"tracr",&val);
  

  val.i=(tg->fldr)[nb];
  puthdval(&(otr->tr),"fldr",&val);
  
  val.i=(tg->tracf)[nb];
  puthdval(&(otr->tr),"tracf",&val);
  
  /* interpolate based on flag. Note that if internal and 
     output grids are same, interp is a no-op, so interpolating
     in error is harmless. */
  if (tg->interp) {
    if ( (err=cubic_(&(tg->t0),&(tg->dt),&((tg->buf)[nb*tg->nt]),&(tg->nt),
		     &(tg->t0out),&(tg->dtout),(otr->tr).data,&(tg->ntout),
		     &iend,work,&wl)) ) { 
      fprintf(stderr,"ERROR: tracegeom::assembletraces from cubic\n");
#ifdef IWAVE_USE_MPI
      MPI_Abort(wcomm,err);
#else
      return err;
#endif
    }
  }
  /* otherwise simple copy, with obvious sanity check */
  else {
    if (tg->nt != tg->ntout) return E_OTHER;
    memcpy((otr->tr).data,&((tg->buf)[nb*tg->nt]),tg->nt*sizeof(float));
  } 
  
  otr->m=(tg->troff)[nb];

  return err;
}


/* 10.07.08: added og back in - note: GLOBAL */
/* 15.01.09: pass pointer to tg, avoid copy */
int writetraces(tracegeom const * tg, 
		RPNT d,
		RPNT og,
		FILE * stream) {
  
  /* BEGIN DECLARE WORKSPACE **************************************************/
  
  offsegy otr;                  /* segy+off_t struct for trace comm */
#ifdef IWAVE_USE_MPI
  segy tr;                      /* trace workspace for accumulation */
#endif
  int nb;                       /* trace number counter */

  int * init;                   /* flag array - need to re-read */

  float dz = 1.0;               /* grid step in z */
  float dx = 1.0;               /* grid step in x */
  float dy = 1.0;               /* grid step in y */
  float ozg = 0.0;              /* global grid origin - z */
  float oxg = 0.0;              /* global grid origin - x */
  float oyg = 0.0;              /* global grid origin - y */
  int idz = 0;                  /* z axis index */
  int idx = 0;                  /* x axis index */
  int idy = 0;                  /* y axis index */

  int j;                        /* general purpose counter */

  int ndim;                     /* problem dimension (should be 1, 2, or 3!) */
  
  int err=0;                    /* error flag */
  float * work;                 /* interpolation workspace */
  int wl;                       /* length of interpolation workspace */

#ifdef IWAVE_USE_MPI
  int tleft;                    /* total number of traces left to record */
  MPI_Status stat;       

  MPI_Comm wcomm;               /* global communicator */
  MPI_Comm lcomm;               /* local communicator */
  int lsize;                    /* size of the local communicator */
#endif

  int lrank;                    /* rank in the local communicator  */

  /* END DECLARE WORKSPACE **************************************************/

  /* local variables for mpi info */
#ifdef IWAVE_USE_MPI
  wcomm=retrieveGlobalComm();
  lcomm=retrieveComm();
  lsize=retrieveSize();
#endif
  lrank=retrieveRank();

  ndim = tg->ndim;
  if (ndim < 1 || ndim > RARR_MAX_NDIM ) {
    fprintf(stream,"ERROR: writetraces\n");
    fprintf(stream,"ndim=%d out of range [1,%d]\n",ndim, RARR_MAX_NDIM);
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_OTHER);
#else
    return E_OTHER;
#endif
  }

  if (tg->ntr[tg->irec] > MAX_TRACES-1) {
    fprintf(stream,"ERROR: writetraces\n");
    fprintf(stream,"number of traces = %d in record %d exceeds hard-coded limit = %d\n",
	    tg->ntr[tg->irec],
	    tg->irec,
	    MAX_TRACES);
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_OTHER);
#else
    return E_OTHER;
#endif
  }

  if (tg->ntraces && (!(tg->buf))) {
    fprintf(stream,"ERROR: writetraces\n");
    fprintf(stream,"number of traces buffered %d > 0 but buffer not allocated\n",
	    tg->ntraces);
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_OTHER);
#else
    return E_ALLOC;
#endif
  }

  /* allocate init array */
  if (!(init=(int *)usermalloc_(MAX_TRACES*sizeof(int)))) {
    fprintf(stream,"Error: init_tracegeom\n");
    fprintf(stream,"failed to allocate %d ints in init array\n",MAX_TRACES);
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_ALLOC);
#else
    return E_ALLOC;
#endif
  }

  for (j=0;j<MAX_TRACES;j++) init[j]=0;  

  if (ndim > 0)  { idz=(tg->axord)[0]; dz = d[(tg->axord)[0]]; ozg = og[(tg->axord[0])]; }
  if (ndim > 0)  { idx=(tg->axord)[1]; dx = d[(tg->axord[1])]; oxg = og[(tg->axord[1])]; }
  if (ndim > 0)  { idy=(tg->axord)[2]; dy = d[(tg->axord)[2]]; oyg = og[(tg->axord)[2]]; }

  /* allocate workspace for interpolation, even if not needed */
  wl=cubic_getworksize(tg->nt);
  work=(float *)usermalloc_(sizeof(float)*wl);
  if (!work) {
    fprintf(stream,"ERROR: tracegeom::writetraces\n");
    fprintf(stream,"failed to allocate workspace for interpolation\n");
    fflush(stream);
#ifdef IWAVE_USE_MPI
    MPI_Abort(wcomm,E_OTHER);
#else
    return E_OTHER;
#endif
  }

  /* first pass: the only pass in serial, rank 0 in parallel */
#ifdef IWAVE_USE_MPI
  if (lrank==0) {
#endif

    /* seek to file begin */
    err=fseeko(tg->fpout,0L,SEEK_SET);
    if (err) {
      fprintf(stream,"ERROR: writetraces from fseeko 1, err=%d\n",err);
      fflush(stream);
#ifdef IWAVE_USE_MPI
      MPI_Abort(wcomm,err);
#else
      return err;
#endif
    }

    //##########
#ifdef IWAVE_VERBOSE
    fprintf(stream,"\nwritetraces -> ");
    iwave_fprintall(stream);

    fprintf(stream,"\nwritetraces: first pass trace loop rk 0\n");
#endif
    for (nb=0;nb<tg->ntraces;nb++) {
#ifdef IWAVE_VERBOSE 
      fprintf(stream,"writetraces: assemble trace %d\n",nb);
#endif  
      err=assembletrace(tg,&otr,nb,
			dz,dx,dy,ozg,oxg,oyg,idz,idx,idy,
			work,wl);
      if (err) {
	fprintf(stream,"ERROR: writetraces from assembletrace, err=%d\n",err);
	fflush(stream);
#ifdef IWAVE_USE_MPI
	MPI_Abort(wcomm,err);
#else
	return err;
#endif	
      }
	  
#ifdef IWAVE_VERBOSE
      fprintf(stream,"writetraces: rk 0 write trace %d tracr=%d offset=%ld\n",nb,otr.tr.tracr,otr.m);
#endif
      err=fseeko(tg->fpout,otr.m,SEEK_SET);
      if (err) {
	fprintf(stream,"ERROR: writetraces from fseeko 2, err=%d\n",err);
	fflush(stream);
#ifdef IWAVE_USE_MPI
	MPI_Abort(wcomm,err);
#else
	return err;
#endif	
      }      
      fputtr(tg->fpout,&(otr.tr));
      fflush(tg->fpout);
      //      for (j=0;j<otr.tr.ns;j++) fprintf(stream,"it=%d data=%g\n",j,otr.tr.data[j]);
      init[otr.tr.tracr]=1;
    }
    /*
    scanf("%d",&err);
    */
#ifdef IWAVE_USE_MPI
  }

  /* second pass - from here on rank 0 acts as server */

#ifdef IWAVE_VERBOSE 
  fprintf(stream,"writetraces: collect traces info \n");
#endif

  /* collect ntraces info */
  if(lsize>1){
    tleft=0;
    MPI_Reduce((void *)(&(tg->ntraces)),(void *)(&tleft),1,MPI_INT,MPI_SUM,0,lcomm);
  }
  else
    tleft=tg->ntraces;

#ifdef IWAVE_VERBOSE
  fprintf(stream,"writetraces: total number of traces = %d\n",tleft);
#endif
 
  if (lrank==0) {

    /* account for traces already written */
    tleft-=tg->ntraces;

    while (tleft > 0) {
     
#ifdef IWAVE_VERBOSE
      fprintf(stream,"writetraces: number of traces yet to be recd = %d\n",tleft);
#endif
      if(lsize>1)
	err=err || MPI_Recv(&otr,1,tg->p,MPI_ANY_SOURCE,0,lcomm,&stat);

      if (err) {
	fprintf(stream,"ERROR: writetraces from MPI_Recv\n");
	fflush(stream);
#ifdef IWAVE_USE_MPI
	MPI_Abort(wcomm,err);
#else
	return err;
#endif	
      }
#ifdef IWAVE_VERBOSE
      fprintf(stream,"RECV from rank %d offset=%ld tracr=%d ns=%d\n",stat.MPI_SOURCE,otr.m,otr.tr.tracr,otr.tr.ns); 
#endif
      /* two different options, depending on whether any part of this
	 trace has been written already. If so, read it back into buffer,
	 add data samples to trace passed via comm. In either case, write 
	 out.
      */
      if (init[otr.tr.tracr]) {
#ifdef IWAVE_VERBOSE
	fprintf(stream,"writetraces: init flag set, read to update\n");
#endif
	err=fseeko(tg->fpout,otr.m,SEEK_SET);
	if (err) {
	  fprintf(stream,"ERROR: writetraces from fseeko 3, err=%d\n",err);
	  fflush(stream);
#ifdef IWAVE_USE_MPI
	  MPI_Abort(wcomm,err);
#else
	  return err;
#endif	
	}
	if (HDRBYTES+(tg->ntout)*sizeof(float) !=fgettr(tg->fpout,&tr)) {
	  fprintf(stream,"ERROR: writetraces from fgettr, failed to read %ld bytes\n", HDRBYTES+(tg->ntout)*sizeof(float));
	  fflush(stream);
#ifdef IWAVE_USE_MPI
	  MPI_Abort(wcomm,err);
#else
	  return err;
#endif	
	}
	for (j=0;j<tr.ns;j++) 
	  otr.tr.data[j]+=tr.data[j];  	 
      }
      err=fseeko(tg->fpout,otr.m,SEEK_SET);
      if (err) {
	fprintf(stream,"ERROR: writetraces from fseeko 4, err=%d\n",err);
	fflush(stream);
#ifdef IWAVE_USE_MPI
	MPI_Abort(wcomm,err);
#else
	return err;
#endif	
      }
      fputtr(tg->fpout,&(otr.tr));       
      fflush(tg->fpout);
      init[otr.tr.tracr]=1;
      tleft--;
#ifdef IWAVE_VERBOSE
      fprintf(stream,"trace %d written, %d left\n",otr.tr.tracr,tleft);
#endif
    
    }

  }

  else {

#ifdef IWAVE_VERBOSE
    fprintf(stream,"writetraces: total number of traces to send = %d\n",tg->ntraces);
#endif
  
    nb=0;
    while (nb<tg->ntraces && !err) {
      err=assembletrace(tg,&otr,nb,
			dz,dx,dy,ozg,oxg,oyg,idz,idx,idy,
			work,wl);
      if (err) {
	fprintf(stream,"ERROR: writetraces from assembletraces\n");
	fflush(stream);
#ifdef IWAVE_USE_MPI
	MPI_Abort(wcomm,err);
#else
	return err;
#endif	
      }
#ifdef IWAVE_VERBOSE
      fprintf(stream,"SEND trace %d of %d from rank=%d offset=%ld ns=%d\n",nb,tg->ntraces,retrieveRank(),otr.m,otr.tr.ns);
      fflush(stream);
#endif
      if(lsize>1)
	err=err || MPI_Send(&otr,1,tg->p,0,0,lcomm);
      
      if (err) {
	fprintf(stream,"ERROR; writetraces from MPI_Send\n");
	fflush(stream);
#ifdef IWAVE_USE_MPI
	MPI_Abort(wcomm,err);
#else
	return err;
#endif 
      }
      nb++;      
    }

  }
#endif 

  userfree_(work);
  userfree_(init);

  /* flush output unit */
  if (lrank==0) fflush(tg->fpout);
  
#ifdef IWAVE_VERBOSE
  fprintf(stream,"at end of writetraces rk=%d err=%d\n",retrieveRank(),err);
  fflush(stream);
#endif

#ifdef IWAVE_USE_MPI
  MPI_Barrier(lcomm);
#endif

  return err;
}


