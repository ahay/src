#include "giwave.h"

int giwave_zero_recv(IWAVE * pstate,
		     GEN_UPDATE_FUN ud,
		     FILE * stream) {

  // no-op unless MPI is defined

#ifdef IWAVE_USE_MPI

  int i;         /* neighbor index */
  int ia;        /* array index */
  RARR * aptr;

  for (ia=0;ia<RDOM_MAX_NARR;ia++) {
    if (ud(ia,(pstate->model).tsind.iv,&(pstate->model))) {
      for ( i = 0; i < (pstate->model).nnei; ++i ) {
	if ( (pstate->pinfo).reinfo[ia][i].type != MPI_DATATYPE_NULL ) {
	  aptr = &(((((pstate->model).ld_r)[i])._s)[ia]);
	  if (pstate->printact > 5) { 
	    fprintf(stream,"giwave_zero_recv: zeroing array %d recv buf %d\n",ia,i);
	    ra_dump(aptr,stream);
	  }
	  ra_zero(aptr);
	}
      }
    }
  }

#endif

  return 0;
}

int giwave_synch(IWAVE * pstate, 
		 GEN_UPDATE_FUN ud,
		 int fwd,
		 FILE * stream) {

  // no-op unless MPI is defined

#ifdef IWAVE_USE_MPI

  int err = 0;
  // step info for the state to be updated

  int it;       /* step index */
  int iv;       /* substep index */
  int ia;       /* array index */

  int i;
  MPI_Status status;  
  int tmpdest, tmpsource;
  MPI_Datatype tmpdest_dt, tmpsource_dt;
  int tmpdest_val, tmpsource_val;
  void *tmpdest_buf, *tmpsource_buf;
  double time;
  time = 0.0; /* To avoid "unitialized variable" warning */

  /* we use the time step internal index for the perturbed field because
     this is the field being updated - on call should be same as index for
     unperturbed field, but don't use this */
  it = (pstate->model).tsind.it;
  iv = (pstate->model).tsind.iv;
  
  if ( pstate->printact > 5 ) {
    fprintf(stream,"\n------ giwave_synch: before exchange, step %d substep %d\n",it,iv);
    for (ia=0;ia<RDOM_MAX_NARR;ia++) {
      if (ud(ia,iv,&(pstate->model))) {
	fprintf(stream,"------ iarr = %d\n",ia);
	rd_print(&((pstate->model).ld_a), ia, stream);
      }
    }
    fflush(stream); 
  }

  for (ia=0;ia<RDOM_MAX_NARR;ia++) {
    if (ud(ia,iv,&(pstate->model))) {

      if ( pstate->printact > 1 ) {
	fprintf(stream,"\n------ giwave_synch fwd=%d array=%d -------------\n",fwd,ia);
	fflush(stream); 
      }

      IPNT gs;
      IPNT ge;
      RARR rsave;
      ra_setnull(&rsave);

      for ( i = 0; i < (pstate->model).nnei; ++i ) {
	// create send data - corresponds to ld_s	
	if ( (pstate->pinfo).seinfo[ia][i].type != MPI_DATATYPE_NULL ) {  
	  tmpdest = (pstate->pinfo).sranks[i];
	  tmpdest_buf = (pstate->pinfo).seinfo[ia][i].buf;
	  tmpdest_dt = (pstate->pinfo).seinfo[ia][i].type;
	  // adjoint case - save a copy of send RARR, which becomes recv buffer
	  // for adjoint
	  if (!fwd) {
	    rd_gse(&(((pstate->model).ld_s)[i]),ia,gs,ge);
	    if (err=ra_create(&rsave,(pstate->model).g.dim,gs,ge)) {
	      fprintf(stream,"\nError: giwave_synch from ra_create err=%d\n",err);
	      return err;
	    }
	    ra_zero(&rsave);
	    if (err=ra_copy(&rsave,&(((((pstate->model).ld_s)[i])._s)[ia]))) {
	      fprintf(stream,"\nError: giwave_synch from ra_copy err=%d\n",err);
	      return err;
	    }
	  }
	}
	else {
	  tmpdest = MPI_PROC_NULL;
	  tmpdest_buf = &tmpdest_val;
	  tmpdest_dt = MPI_INT;
	} 
	// create receive data - corresponds to ld_r
	if ( (pstate->pinfo).reinfo[ia][i].type != MPI_DATATYPE_NULL ) {
	  tmpsource = (pstate->pinfo).rranks[i];
	  tmpsource_buf = (pstate->pinfo).reinfo[ia][i].buf;
	  tmpsource_dt = (pstate->pinfo).reinfo[ia][i].type;
	} 
	else {
	  tmpsource = MPI_PROC_NULL;
	  tmpsource_buf = &tmpsource_val;
	  tmpsource_dt = MPI_INT;
	}
	  
	if (pstate->printact > 1) {
	  fprintf(stream, "    i = %d, sending to wrk = %d, receiving from wrk = %d, [NULL = %d]\n", 
		  i, tmpdest, tmpsource, MPI_PROC_NULL); 
          fflush(stream);
	}
	
	/* 
	   "dest" is send buffer
	   "source" is receive buffer

	   if fwd: 
	   SEND data in tmpdest_buf to rk=tmpdest
	   RECEIVE data in tmpsource_buf from rk=tmpsource

	   else:
	   SEND data in tmpsource_buf to rk=tmpsource
	   RECEIVE data in tmpdest_buf from rk=tmpdest

	*/
	if (fwd) {
	  err = MPI_Sendrecv(tmpdest_buf, 1, tmpdest_dt, tmpdest, iv,
			     tmpsource_buf, 1, tmpsource_dt, tmpsource, iv,
			     (pstate->pinfo).ccomm, &status);
	  if ( err != MPI_SUCCESS ) {
	    fprintf(stream, 
		    "ERROR. Internal: MPI_Sendrecv error #%d, nei=%d, iv=%d, arr=%d. ABORT.\n", 
		    err, i, iv, ia);
	    return err;
	  }
	}
	else {
	  // first receive into send buffer, which has been copied into a tmp buffer
	  // reciprocally, send receive buffer

	  err = MPI_Sendrecv(  tmpsource_buf, 1, tmpsource_dt, tmpsource, iv,
			       tmpdest_buf, 1, tmpdest_dt, tmpdest, iv,
			       (pstate->pinfo).ccomm, &status);

	  if ( err != MPI_SUCCESS ) {
	    fprintf(stream, 
		    "ERROR. Internal: MPI_Sendrecv error #%d, nei=%d, iv=%d, arr=%d. ABORT.\n", 
		    err, i, iv, ia);
	    return err;
	  }

	  // then add tmp buffer back
	  if ( ( (pstate->pinfo).seinfo[ia][i].type != MPI_DATATYPE_NULL ) ) {
	    if (!(rsave._s0)) {
	      fprintf(stream,"\nError: giwave_synch before axpy: rsave not initialized\n");
	      return E_INTERNAL;
	    }
	    fflush(stream);
	    if (err=ra_axpy(&(((((pstate->model).ld_s)[i])._s)[ia]),&rsave,REAL_ONE)) {
	      fprintf(stream,"\nError: giwave_synch from ra_axpy err=%d\n",err);
	      ra_dump(&(((((pstate->model).ld_r)[i])._s)[ia]),stream);
	      ra_dump(&rsave,stream);
	      return err;
	    }
	    ra_destroy(&rsave);
	  }

	  // now comp subdom data is correct - synchronize
	  err = MPI_Sendrecv(tmpdest_buf, 1, tmpdest_dt, tmpdest, iv,
			     tmpsource_buf, 1, tmpsource_dt, tmpsource, iv,
			     (pstate->pinfo).ccomm, &status);
	  if ( err != MPI_SUCCESS ) {
	    fprintf(stream, 
		    "ERROR. Internal: MPI_Sendrecv error #%d, nei=%d, iv=%d, arr=%d. ABORT.\n", 
		    err, i, iv, ia);
	    return err;
	  }

	}
      }
    }
  }

  if ( pstate->printact > 5 ) {
    fprintf(stream,"\n------ giwave_synch: after exchange, step %d substep %d\n",it,iv);
    for (ia=0;ia<RDOM_MAX_NARR;ia++) {
      if (ud(ia,iv,&(pstate->model))) {
	fprintf(stream,"------ iarr = %d\n",ia);
	rd_print(&((pstate->model).ld_a), ia, stream);
      }
    }
    fflush(stream); 
  }
#endif

  return 0;
}

int giwave_dynamic_init(IWAVE * state,
			int _it, int _itoff) {
  
  char * filename = NULL; /* workspace for input data filename */
  int i;
  int err = 0;    /* error flag */
  int rk = 0;     /* default proc rank = 0 */
  FD_MODEL * fdm;

#ifdef IWAVE_USE_MPI
  rk = retrieveGlobalRank();
#endif
    
  //  fprintf(stderr,"giwave_dynamic_init: it=%d itoff=%d\n",_it,_itoff);
  /* initialize dynamic field at _it */

  fdm = (FD_MODEL *)((state->model).specs);

  if (_itoff==0) {
    iwave_dynamic_init(state,_it);
  }
  else if (_itoff>0) {
    /* find input data file and read data into dynamic arrays*/   
    for (i=0;i<(state->model).ld_a.narr; i++) {
      //      if (isdyn(fdm,i)) {
      if (fd_isdyn(i)) {
	filename = iwave_get_a_checkfilename(i,_itoff);
	if (!filename) { 
	  err = 101;
	  fprintf(stderr,"\nERROR. In giwave_dynamic_init, proc %d, filename is NULL for array %d.\n",rk,i);
	  return err;
	}
	err = ra_fread(&((state->model).ld_a._s[i]), filename);
	if(filename) free(filename);
	if (err){ 
	  fprintf(stderr,"\nERROR. In giwave_dynamic_init, proc %d, ra_fread for array %d return err=%d for read %s.\n",rk,i,err,filename);
	  return err;
	}
      }
    }
    /* transfer check time to state */
    (state->model).tsind.it=_it;
    (state->model).tsind.iv=0;  
  }
  
  else /*itoff<0*/ {
    err = 1;
    fprintf(stderr,"\n ERROR. In giwave_dynamic_init, proc %d, itoff = %d < 0.\n",rk,_itoff);
  }
  
  return err;
}

/*----------------------------------------------------------------------------*/

int iwave_dynamic_takeshot(IWAVE * state, 
			   int _itoff ) /* _itoff=_itcheck - _itstart */
{
 
  char * filename = NULL; /* workspace for output data filename */
  int err = 0;    /* error flag */
  int iarr;
  int rk = 0;   /* default proc rank = 0 */
  FD_MODEL * fdm;

#ifdef IWAVE_USE_MPI
  rk = retrieveGlobalRank();
#endif

  fdm = (FD_MODEL *)((state->model).specs);  

  /* write dynamic array at _it into checkfile */
  for (iarr=0; iarr<(state->model).ld_a.narr; iarr++) {
    //    if (isdyn(fdm,iarr)) {
    if (fd_isdyn(iarr)) {
      filename = iwave_get_a_checkfilename(iarr,_itoff);
      if (!filename) { 
	err = 101;
	fprintf(stderr,"\nERROR. In iwave_dynamic_takeshot, proc %d, filename is NULL for array %d.\n",rk,iarr);
	return err;
      }
      err = ra_fwrite(&((state->model).ld_a._s[iarr]),filename); 
      if(filename) free(filename); 
      if (err) {
	fprintf(stderr,"\nERROR. In iwave_dynamic_takeshot, proc %d, ra_fwrite for array %d return err=%d.\n",rk,iarr,err);
	return err;
      }
    }
  }
  return err;
}


/*----------------------------------------------------------------------------*/
int iwave_remove_checkfile(IWAVE * state, int _itoff){
 
  char * filename = NULL; /* workspace for data filename */
  int err = 0;    /* error flag */
  int iarr;
  FD_MODEL * fdm;
  /* default proc rank = 0 
     int rk = 0;
  */  
  /*
    #ifdef IWAVE_USE_MPI
    rk = retrieveGlobalRank();
    #endif
  */
  fdm = (FD_MODEL *)((state->model).specs);  
  for (iarr=0; iarr<(state->model).ld_a.narr; iarr++) {
    //    if (isdyn(fdm,iarr)) {
    if (fd_isdyn(iarr)) {
      filename = iwave_get_a_checkfilename(iarr,_itoff);
      if (!filename) { 
	/*
	  fprintf(stderr,"\n Warning. In iwave_remove_checkfile, proc %d, filename is NULL for array %d.\n",rk,iarr);
	*/
      }
      if( remove(filename)){    
	/*
	  fprintf(stderr,"\n Warning. In  iwave_remove_checkfile, file %s not deleted. Is it exist?\n",filename);
	  return err;
	*/
      }
      if(filename) free(filename); 
    }
  }
  return err;
}

/*----------------------------------------------------------------------------*/

char * iwave_get_a_checkfilename(int iarr, 
				 int _itoff ) /* _itoff=_itcheck - _itstart */
{ 
  char * filename = NULL; /* workspace for output data filename */
  char tmp[100];  
  int rk = 0;   /* default proc rank = 0 */  

#ifdef IWAVE_USE_MPI
  rk = retrieveGlobalRank();
#endif

  /* generate checkfile names */
  if (getenv("DATAPATH")) {
    filename = (char*)malloc((strlen(getenv("DATAPATH")) + 20) * sizeof(char));
    if ( filename == NULL ) {
      fprintf(stderr,"\n ERROR. In iwave_get_a_checkfilename, proc %d, fail to create file path. \n",rk);
      return filename;
    }
    strcpy(filename,getenv("DATAPATH"));
  }
  else {
    filename = (char*)malloc(MAXPATHLEN * sizeof(char));    
    if ( filename == NULL ) {
      fprintf(stderr,"\n ERROR. In iwave_get_a_checkfilename, proc %d, memory allocation for filename.\n",rk);
      return filename;
    }
    if (!getcwd(filename,MAXPATHLEN)) {
      fprintf(stderr,"\nERROR. In iwave_get_a_checkfilename, proc %d, getcwd failed for cout filename.\n",rk);
      return filename;
    }
  }
  
  if (strlen(filename)>1 && (filename[strlen(filename)-1] != '/')) 
    strcat(filename,"/");

  sprintf(tmp,"statecheckpoint_%d_proc_%d_array_%d.bin",_itoff,rk,iarr);
  strcat(filename,tmp);
  
  return filename;  

}

 
