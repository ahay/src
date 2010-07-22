/* Generic conjugate-gradient solver for linear inversion with complex data */
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

#define DLOOP(a) for (dsiz=nd, dbuf=nbuf; dsiz > 0; dsiz -= dbuf) { \
	         if (dsiz < dbuf) dbuf=dsiz; {a} }

#define MLOOP(a) for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) { \
	         if (msiz < mbuf) mbuf=msiz; {a} }

#define DREAD(a) if (dbuf != fread(buf,sizeof(sf_complex),dbuf,a)) \
                  sf_error("write error")

#define MREAD(a) if (mbuf != fread(buf,sizeof(sf_complex),mbuf,a)) \
                  sf_error("write error")

#define DREAD2(a) if (dbuf != fread(buf2,sizeof(sf_complex),dbuf,a)) \
                  sf_error("write error")

#define MREAD2(a) if (mbuf != fread(buf2,sizeof(sf_complex),mbuf,a)) \
                  sf_error("write error")

#define DWRITE(a) if (dbuf != fwrite(buf,sizeof(sf_complex),dbuf,a)) \
                  sf_error("write error")

#define MWRITE(a) if (mbuf != fwrite(buf,sizeof(sf_complex),mbuf,a)) \
                  sf_error("write error")

#define DUP(a) if (dup(a) < 0) sf_error("dup error:")

int main(int argc, char* argv[])
{
    int i, iter, niter, p[6][2], status;
    sf_complex *buf, *buf2;
    double rn, rnp, alpha, beta;
    pid_t pid[6]={1,1,1,1,1,1};
    off_t nm, nd, msiz, dsiz, pos;
    size_t nbuf, mbuf, dbuf;
    FILE *xfile, *Rfile, *gfile, *sfile, *Sfile;
    char *x, *R, *g, *s, *S;
    sf_file mod=NULL;
    sf_file dat=NULL;
    sf_file out=NULL;
    sf_file from=NULL;
    sf_file to=NULL;
    extern int fseeko(FILE *stream, off_t offset, int whence);
    extern off_t ftello (FILE *stream);

    sf_init(argc,argv);
    dat = sf_input("in");
    mod = sf_input("mod");

    if (SF_COMPLEX != sf_gettype(mod) ||
	SF_COMPLEX != sf_gettype(dat)) 
	sf_error("Need complex type in mod and dat");

    for (i=0; i < argc-1; i++) {
	argv[i]=argv[i+1];
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
    fclose(gfile);
    fclose(sfile);
    fclose(Sfile);

    nm = sf_filesize(mod);
    nd = sf_filesize(dat);

    /* I/O buffers */
    nbuf = BUFSIZ/sizeof(sf_complex);
    buf  = sf_complexalloc(nbuf);
    buf2 = sf_complexalloc(nbuf);


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
	    sf_fileflush(to,dat);
	    sf_settype(to,SF_COMPLEX);

	    if (0 == iter) {
		xfile = fopen(x,"wb");
		for (i=0; i < nbuf; i++) { buf[i] = sf_cmplx(0.,0.); }
		MLOOP( MWRITE(xfile); );
		fclose(xfile);

		Rfile = fopen(R,"wb");
#ifdef SF_HAS_COMPLEX_H 
		DLOOP( sf_complexread(buf,dbuf,dat); 
		       for (i=0; i < dbuf; i++) { buf[i] = -buf[i]; }
		       sf_complexwrite(buf,dbuf,to);
		       DWRITE (Rfile); );
#else
		DLOOP( sf_complexread(buf,dbuf,dat); 
		       for (i=0; i < dbuf; i++) { buf[i] = sf_cneg(buf[i]); }
		       sf_complexwrite(buf,dbuf,to);
		       DWRITE (Rfile); );
#endif
	    } else {
		Rfile = fopen(R,"rb");
		DLOOP( DREAD(Rfile); sf_complexwrite(buf,dbuf,to); );
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
	    sf_settype(to,SF_COMPLEX);

	    gfile = fopen(g,"w+b");

	    /* rn = ddot(g) */
	    
	    rn = 0.;

#ifdef SF_HAS_COMPLEX_H
	    MLOOP( sf_complexread(buf,mbuf,from);
		   for (i=0; i < mbuf; i++) { 
		       rn += creal(conj(buf[i]) * buf[i]); }
		   MWRITE(gfile); );
#else
	    MLOOP( sf_complexread(buf,mbuf,from);
		   for (i=0; i < mbuf; i++) { 
		       rn += 
			   (double) buf[i].r * buf[i].r +
			   (double) buf[i].i * buf[i].i; }
		   MWRITE(gfile); );
#endif

	    sfile = fopen(s,"r+b");

	    if (0==iter) {
		alpha = 0.;
	    } else {
		if (1 != fread(&rnp,sizeof(double),1,sfile)) 
		    sf_error("read error");

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

#ifdef SF_HAS_COMPLEX_H	    
	    MLOOP( MREAD(gfile); sf_complexwrite(buf,mbuf,to);

		   if (iter > 0) {
		       pos = ftello(sfile);

		       MREAD2(sfile); 
		       
		       if (0 > fseeko(sfile,pos,SEEK_SET))
			   sf_error ("seek problem");

		       for (i=0; i < mbuf; i++) {
			   buf[i] += alpha * buf2[i];
		       }
		   } 
		   
		   MWRITE(sfile); );
#else
	    MLOOP( MREAD(gfile); sf_complexwrite(buf,mbuf,to);

		   if (iter > 0) {
		       pos = ftello(sfile);

		       MREAD2(sfile); 
		       
		       if (0 > fseeko(sfile,pos,SEEK_SET))
			   sf_error ("seek problem");

		       for (i=0; i < mbuf; i++) {
			   buf[i] = sf_cadd(buf[i],sf_crmul(buf2[i],alpha));
		       }
		   } 
		   
		   MWRITE(sfile); );
#endif

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

#ifdef SF_HAS_COMPLEX_H
	    DLOOP( sf_complexread(buf,dbuf,from);

		   if (iter > 0) {
		       pos = ftello(Sfile);

		       DREAD2(Sfile);
		       
		       if (0 > fseeko(Sfile,pos,SEEK_SET))
			   sf_error ("seek problem");
		    
		       for (i=0; i < dbuf; i++) {
			   buf[i] += alpha * buf2[i];
		       }
		   }

		   for (i=0; i < dbuf; i++) {
		       beta += creal(conjf(buf[i]) * buf[i]);
		   }
		
		   DWRITE(Sfile); );
#else
	    DLOOP( sf_complexread(buf,dbuf,from);
		   
		   if (iter > 0) {
		       pos = ftello(Sfile);
		       
		       DREAD2(Sfile);
		       
		       if (0 > fseeko(Sfile,pos,SEEK_SET))
			   sf_error ("seek problem");
		    
		       for (i=0; i < dbuf; i++) {
			   buf[i] = sf_cadd(buf[i],sf_crmul(buf2[i],alpha));
		       }
		   }

		   for (i=0; i < dbuf; i++) {
		       beta += 
			   (double) buf[i].r * buf[i].r +
			   (double) buf[i].i * buf[i].i;
		   }
		
		   DWRITE(Sfile); );
#endif

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

	    if (sizeof(double) !=  read(p[5][0],&alpha,sizeof(double)))
		sf_error("read error");

	    sfile = fopen(s,"rb"); 
	    if (NULL == sfile) sf_error("Cannot open %s:",s);

	    if (0 > fseeko(sfile,sizeof(double),SEEK_SET))
		sf_error ("seek problem");

	    xfile = fopen(x,"r+b"); 
	    if (NULL == xfile) sf_error("Cannot open %s:",x);

#ifdef SF_HAS_COMPLEX_H
	    MLOOP( pos = ftello(xfile); 

		   MREAD(xfile);
		   MREAD2(sfile); 
		   
		   for (i=0; i < mbuf; i++) {
		       buf[i] += alpha * buf2[i];
		   }

		   if (0 > fseeko(xfile,pos,SEEK_SET))
		   sf_error ("seek problem");
	    
		   MWRITE(xfile); );
#else
	    MLOOP( pos = ftello(xfile); 

		   MREAD(xfile);
		   MREAD2(sfile); 
		   
		   for (i=0; i < mbuf; i++) {
		       buf[i] = sf_cadd(buf[i],sf_crmul(buf2[i],alpha));
		   }

		   if (0 > fseeko(xfile,pos,SEEK_SET))
		   sf_error ("seek problem");
	    
		   MWRITE(xfile); );
#endif

	    fclose(sfile);
	    fclose(xfile);

	    Sfile = fopen(S,"rb");
	    Rfile = fopen(R,"r+b");

#ifdef SF_HAS_COMPLEX_H
	    DLOOP( pos = ftello(Rfile); 

		   DREAD(Rfile);
		   DREAD2 (Sfile);
		   
		   for (i=0; i < dbuf; i++) {
		       buf[i] += alpha * buf2[i];
		   }

		   if (0 > fseeko(Rfile,pos,SEEK_SET))
		   sf_error ("seek problem");
		   
		   DWRITE(Rfile); );
#else
	    DLOOP( pos = ftello(Rfile); 

		   DREAD(Rfile);
		   DREAD2 (Sfile);
		   
		   for (i=0; i < dbuf; i++) {
		       buf[i] = sf_cadd(buf[i],sf_crmul(buf2[i],alpha));
		   }

		   if (0 > fseeko(Rfile,pos,SEEK_SET))
		   sf_error ("seek problem");
		   
		   DWRITE(Rfile); );
#endif

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
    sf_settype(out,SF_COMPLEX);

    xfile = fopen(x,"rb");
    MLOOP( MREAD(xfile); sf_complexwrite(buf,mbuf,out); );
    fclose(xfile);

    unlink(R);
    unlink(x);
    unlink(g);
    unlink(s);
    unlink(S);


    exit(0);
}
