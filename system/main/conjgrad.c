/* Generic conjugate-gradient solver for linear inversion */
/*
  Copyright (C) 2005 University of Texas at Austin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#include <sys/types.h>
#include <unistd.h>
/*^*/

#include <stdio.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <rsf.h>

#define DLOOP(a) for (dsiz=nd, dbuf=nbuf; dsiz > 0; dsiz -= dbuf) {	\
	if (dsiz < dbuf) dbuf=dsiz; {a} }

#define MLOOP(a) for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {	\
	if (msiz < mbuf) mbuf=msiz; {a} }

#define DREAD(a) if (dbuf != fread(buf,sizeof(float),dbuf,a))	\
	sf_error("write error")

#define MREAD(a) if (mbuf != fread(buf,sizeof(float),mbuf,a))	\
	sf_error("write error")

#define DREAD2(a) if (dbuf != fread(buf2,sizeof(float),dbuf,a)) \
	sf_error("write error")

#define MREAD2(a) if (mbuf != fread(buf2,sizeof(float),mbuf,a)) \
	sf_error("write error")

#define DWRITE(a) if (dbuf != fwrite(buf,sizeof(float),dbuf,a)) \
	sf_error("write error")

#define MWRITE(a) if (mbuf != fwrite(buf,sizeof(float),mbuf,a)) \
	sf_error("write error")

#define DUP(a) if (dup(a) < 0) sf_error("dup error:")

int main(int argc, char* argv[])
{
    int i, iter, niter, p[6][2], status, *mask;
    float *buf, *buf2, *wht;
    double rn, rnp, alpha, beta;
    pid_t pid[6]={1,1,1,1,1,1};
    off_t nm, nd, msiz, dsiz, pos;
    size_t nbuf, mbuf, dbuf;
    FILE *xfile, *Rfile, *gfile, *sfile, *Sfile;
    char *x, *R, *g, *s, *S, *prog;
    sf_file mod, dat, from, mwt, x0, known;  /* input */
    sf_file to, out; /* output */

    extern int fseeko(FILE *stream, off_t offset, int whence);
    extern off_t ftello (FILE *stream);

    sf_init(argc,argv);

    dat = sf_input("in");
    mod = sf_input("mod");

    if (SF_FLOAT != sf_gettype(mod) ||
	SF_FLOAT != sf_gettype(dat)) 
	sf_error("Need float type in mod and dat");

    for (i=0; i < argc-1; i++) {
	argv[i]=argv[i+1];
    }
    for (i=0; i < argc-1; i++) {	
	/* find the program to run */
	if (NULL == strchr(argv[i],'=')) {
	    /* first one without the '=' */
	    prog = argv[0];
	    argv[0] = argv[i];
	    argv[i] = prog;
	    break;
	}
    }

    argv[argc-1] = sf_charalloc(6);
    snprintf(argv[argc-1],6,"adj=X");

    if (!sf_getint("niter",&niter)) niter=1;
    /* number of iterations */

    Rfile = sf_tempfile(&R,"w+b"); 
    xfile = sf_tempfile(&x,"w+b"); 
    gfile = sf_tempfile(&g,"w+b");
    sfile = sf_tempfile(&s,"w+b");
    Sfile = sf_tempfile(&S,"w+b");

    fclose(Rfile);
    fclose(xfile);
    fclose(gfile);
    fclose(sfile);
    fclose(Sfile);

    nm = sf_filesize(mod);
    nd = sf_filesize(dat);

    /* I/O buffers */
    nbuf = BUFSIZ/sizeof(float);
    buf  = sf_floatalloc(nbuf);
    buf2 = sf_floatalloc(nbuf);

    if (NULL != sf_getstring("mwt")) {
	mwt = sf_input("mwt"); /* model weight */
	wht = sf_floatalloc(nbuf);
    } else {
	mwt = NULL;
	wht = NULL;
    }

    if (NULL != sf_getstring("known")) {
	known = sf_input("known"); /* known model mask */
	if (SF_INT != sf_gettype(known)) sf_error("Need int type in known");
	mask = sf_intalloc(nbuf);
    } else {
	known = NULL;
	mask = NULL;
    }

    if (NULL != sf_getstring("x0")) {
	x0 = sf_input("x0"); /* initial model */
    } else {
	x0 = NULL;
    }

    for (i=0; i < 6; i++) { /* make six pipes */
	if (pipe(p[i]) < 0) sf_error("pipe error:");
    }

    for (iter=0; iter < niter; iter++) {
	for (i=0; i < 6; i++) { /* fork six children */
	    if ((pid[i] = fork()) < 0) sf_error("fork error:");
	    if (0 == pid[i]) break;
	}
	
	if (0 == pid[0]) {	
	    /* feeds rr to p[0] */
	    close(p[0][0]);
	    close(STDOUT_FILENO);
	    DUP(p[0][1]);

	    to = sf_output("out");

	    if (0 == iter) {
		xfile = fopen(x,"wb");

		if (NULL == x0) {
		    for (i=0; i < nbuf; i++) { buf[i] = 0.0f; }
		}
		
		MLOOP( if (NULL != x0) sf_floatread(buf,mbuf,x0);
		       MWRITE(xfile); );
		
		fclose(xfile);

		Rfile = fopen(R,"wb");
		DLOOP( sf_floatread(buf,dbuf,dat); 
		       for (i=0; i < dbuf; i++) { buf[i] = -buf[i]; }
		       sf_floatwrite(buf,dbuf,to);
		       DWRITE (Rfile); );
	    } else {
		Rfile = fopen(R,"rb");
		DLOOP( DREAD(Rfile); sf_floatwrite(buf,dbuf,to); );
	    }
	    fclose(Rfile);
	    
	    sf_warning("iter %d of %d",iter+1,niter);

	    _exit(0);
	}
 
	
	if (0==pid[1]) {
	    /* reads from p[0], runs the adjoint, and writes to p[1] */

	    close(p[0][1]);
	    close(STDIN_FILENO);
	    DUP(p[0][0]);
	    
	    close(p[1][0]);
	    close(STDOUT_FILENO);
	    DUP(p[1][1]);

	    argv[argc-1][4]='1';

	    execvp(argv[0],argv);

	    _exit(1);
	}


	if (0==pid[2]) {
	    /* reads gradient from p[1], uses it, and writes to p[2] */
	
	    close(p[1][1]);
	    close(STDIN_FILENO);
	    DUP(p[1][0]);
	    from = sf_input("in");
	    
	    close(p[2][0]);
	    close(STDOUT_FILENO);
	    DUP(p[2][1]);
	    to = sf_output("out");
	    sf_fileflush(to,mod);

	    gfile = fopen(g,"w+b");

	    /* rn = ddot(g) */
	    
	    rn = 0.;

	    if (NULL != mwt) sf_seek(mwt,0,SEEK_SET);
	    if (NULL != known) sf_seek(known,0,SEEK_SET);

	    MLOOP( sf_floatread(buf,mbuf,from);

		   if (NULL != mwt) { 
		       sf_floatread(wht,mbuf,mwt);
		       for (i=0; i < mbuf; i++) { buf[i] *= wht[i]; }
		   }

		   if (NULL != known) { 
		       sf_intread(mask,mbuf,known);
		       for (i=0; i < mbuf; i++) { sf_warning("%d %d",i,mask[i]); if (mask[i]) buf[i] = 0.0f; }
		   }

		   for (i=0; i < mbuf; i++) { rn += (double) buf[i] * buf[i]; }
		   MWRITE(gfile); );
	    
	    sfile = fopen(s,"r+b");

	    if (0==iter) {
		alpha = 0.;
	    } else {
		if (1 != fread(&rnp,sizeof(double),1,sfile)) sf_error("read error");

		alpha = rn/rnp;

		if (0 > fseeko(sfile,0,SEEK_SET))
		    sf_error ("seek problem");
	    }
	    
	    if (sizeof(double) != write(p[4][1],&alpha,sizeof(double)))
		sf_error("write error");
	    
	    fwrite(&rn,sizeof(double),1,sfile);

	    if (0 > fseeko(gfile,0,SEEK_SET))
		sf_error ("seek problem");

	    /* s = g + alpha*s */

	    MLOOP( MREAD(gfile); 		 

		   if (iter > 0) {
		       pos = ftello(sfile);

		       MREAD2(sfile); 
		       
		       if (0 > fseeko(sfile,pos,SEEK_SET))
			   sf_error ("seek problem");

		       cblas_saxpy(mbuf,alpha,buf2,1,buf,1);
		   } 
		   
		   MWRITE(sfile); ); 

	    if (0 > fseeko(gfile,0,SEEK_SET))
		sf_error ("seek problem");

	    if (NULL != mwt) sf_seek(mwt,0,SEEK_SET);
	    
	    MLOOP( MREAD(gfile); 

		   if (NULL != mwt) { 
		       sf_floatread(wht,mbuf,mwt);
		       for (i=0; i < mbuf; i++) { buf[i] *= wht[i]; }
		   }
		   
		   sf_floatwrite(buf,mbuf,to); );
	    
	    fclose(gfile);
	    fclose(sfile);

	    sf_warning("grad=%lg",rn);

	    _exit(2);
	}

	if (0==pid[3]) {
	    /* reads from p[2], runs the forward, and writes to p[3] */

	    close(p[2][1]);
	    close(STDIN_FILENO);
	    DUP(p[2][0]);
	    
	    close(p[3][0]);
	    close(STDOUT_FILENO);
	    DUP(p[3][1]);

	    argv[argc-1][4]='0';
	    execvp(argv[0],argv);

	    _exit(3);
	}

	if (0==pid[4]) {
	    /* reads conj. grad from p[3], does computations with it */

	    close(p[3][1]);
	    close(STDIN_FILENO);
	    DUP(p[3][0]);
	    from = sf_input("in");
	
	    if (sizeof(double) != read(p[4][0],&alpha,sizeof(double))) 
		sf_error("read error");

	    Sfile = fopen(S,"r+b");

	    /* beta = ddot(ss) */

	    beta = 0.;

	    /* ss = gg + alpha * ss */

	    DLOOP( sf_floatread(buf,dbuf,from);

		   if (iter > 0) {
		       pos = ftello(Sfile);

		       DREAD2(Sfile);
		       
		       if (0 > fseeko(Sfile,pos,SEEK_SET))
			   sf_error ("seek problem");
		    
		       cblas_saxpy(dbuf,alpha,buf2,1,buf,1);
		   }

		   beta += cblas_dsdot(dbuf,buf,1,buf,1);
		
		   DWRITE(Sfile); );
	    fclose(Sfile);

	    sfile = fopen(s,"rb");

	    if (1 != fread(&rn,sizeof(double),1,sfile)) sf_error("read error");

	    fclose(sfile);

	    alpha = - rn/beta;

	    if (sizeof(double) != write(p[5][1],&alpha,sizeof(double)))
		sf_error("write error");
	    
	    _exit(4);
	}

	if (0==pid[5]) {
	    /* updates x and r */

	    if (sizeof(double) != read(p[5][0],&alpha,sizeof(double)))
		sf_error("read error");

	    sfile = fopen(s,"rb"); if (NULL == sfile) sf_error("Cannot open %s:",s);

	    if (0 > fseeko(sfile,sizeof(double),SEEK_SET))
		sf_error ("seek problem");

	    xfile = fopen(x,"r+b"); if (NULL == xfile) sf_error("Cannot open %s:",x);

	    MLOOP( pos = ftello(xfile); 

		   MREAD(xfile);
		   MREAD2(sfile); 
		   
		   cblas_saxpy(mbuf,alpha,buf2,1,buf,1);

		   if (0 > fseeko(xfile,pos,SEEK_SET))
		       sf_error ("seek problem");
	    
		   MWRITE(xfile); );

	    fclose(sfile);
	    fclose(xfile);

	    Sfile = fopen(S,"rb");
	    Rfile = fopen(R,"r+b");

	    DLOOP( pos = ftello(Rfile); 

		   DREAD(Rfile);
		   DREAD2 (Sfile);
		   
		   cblas_saxpy(dbuf,alpha,buf2,1,buf,1);

		   if (0 > fseeko(Rfile,pos,SEEK_SET))
		       sf_error ("seek problem");
		   
		   DWRITE(Rfile); );

	    fclose(Sfile);
	    fclose(Rfile);

	    _exit(5);
	}
	
        /* parent waits */
	for (i=0; i < 6; i++) { 
	    if (0 == pid[i]) sf_error("A child alive");
	    waitpid(pid[i],&status,0);
	}
    }
	
    /* write x to out */  
    out = sf_output("out");
    sf_fileflush(out,mod);

    xfile = fopen(x,"rb");

    if (NULL != mwt) sf_seek(mwt,0,SEEK_SET);

    MLOOP( MREAD(xfile); 

	   if (NULL != mwt) { 
	       sf_floatread(wht,mbuf,mwt);
	       for (i=0; i < mbuf; i++) { buf[i] *= wht[i]; }
	   } 

	   sf_floatwrite(buf,mbuf,out); );

    fclose(xfile);

    unlink(R);
    unlink(x);
    unlink(g);
    unlink(s);
    unlink(S);


    exit(0);
}
