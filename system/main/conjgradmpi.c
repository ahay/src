/* Generic conjugate-gradient solver for linear inversion.

In this version, the linear operator program uses --input and --output instead of stdin and stdout.
*/
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
	if (dsiz < dbuf) { dbuf=dsiz; } {a} }

#define MLOOP(a) for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {	\
	if (msiz < mbuf) { mbuf=msiz; } {a} }

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
    int i, iter, niter;
    float *buf, *buf2, *wht;
    double rn, rnp=0.0, alpha, beta;
    off_t nm, nd, msiz, dsiz, pos;
    size_t nbuf, mbuf, dbuf, len, iolen, cmdlen;
    FILE *xfile, *Rfile, *gfile, *Gfile, *sfile, *Sfile;
    char *x, *R, *g, *G, *s, *S, *prog, *cmdline, *iostring, *arg;
    sf_file mod, dat, mwt, x0;  /* input */
    sf_file Rrsf, Grsf, grsf, out; /* output */

    extern int fseeko(FILE *stream, off_t offset, int whence);
    extern off_t ftello (FILE *stream);

    sf_init(argc,argv);

    dat = sf_input("--input");
    mod = sf_input("mod");

    out = sf_output("--output");
    sf_fileflush(out,mod);

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

    Rfile = sf_tempfile(&R,"w+"); 
    xfile = sf_tempfile(&x,"w+b"); 
    gfile = sf_tempfile(&g,"w+");
    Gfile = sf_tempfile(&G,"w+");
    sfile = sf_tempfile(&s,"w+b");
    Sfile = sf_tempfile(&S,"w+b");

    fclose(Rfile); Rrsf = sf_output(R); sf_readwrite(Rrsf,true); sf_fileflush(Rrsf,dat);
    fclose(xfile);
    fclose(gfile); 
    fclose(Gfile); 
    fclose(sfile);
    fclose(Sfile);

    nm = sf_filesize(mod);
    nd = sf_filesize(dat);

    /* I/O buffers */
    nbuf = BUFSIZ/sizeof(float);
    buf  = sf_floatalloc(nbuf);
    buf2 = sf_floatalloc(nbuf);

    cmdline = sf_charalloc(SF_CMDLEN);
    iostring = sf_charalloc(SF_CMDLEN);
    cmdlen = 0;
    for (i=0; i < argc; i++) {
	arg = argv[i];
	len = strlen(arg);
	if (cmdlen+len > SF_CMDLEN-2) sf_error("command line is too long");

	memcpy(cmdline+cmdlen,arg,len);
	cmdline[cmdlen+len]=' ';
	cmdlen += len+1;
    }

    if (NULL != sf_getstring("mwt")) {
	mwt = sf_input("mwt"); /* model weight */
	wht = sf_floatalloc(nbuf);
    } else {
	mwt = NULL;
	wht = NULL;
    }

    if (NULL != sf_getstring("x0")) {
	x0 = sf_input("x0"); /* initial model */
    } else {
	x0 = NULL;
    }

    for (iter=0; iter < niter; iter++) {
	if (0 == iter) {
	    xfile = fopen(x,"wb");
	    
	    if (NULL == x0) {
		for (i=0; i < nbuf; i++) { buf[i] = 0.0f; }
	    }
	    
	    MLOOP( if (NULL != x0) sf_floatread(buf,mbuf,x0);
		   MWRITE(xfile); );
	    
	    fclose(xfile);

	    DLOOP( sf_floatread(buf,dbuf,dat); 
		   for (i=0; i < dbuf; i++) { buf[i] = -buf[i]; }
		   sf_floatwrite(buf,dbuf,Rrsf); );
	    sf_fflush(Rrsf);
	} 
	
	sf_warning("iter %d of %d",iter+1,niter);
	
	/* reads from R, runs the adjoint, and writes to g */
	cmdline[cmdlen-2]='1';
	snprintf(iostring,SF_CMDLEN,"--input=%s --output=%s",R,g);
	iolen = strlen(iostring);
	if (cmdlen+iolen > SF_CMDLEN-1) sf_error("command line is too long");
	snprintf(cmdline+cmdlen,iolen+1,"%s",iostring);

	sf_warning(cmdline);
	sf_system(cmdline);
	
	/* rn = ddot(g) */
	    
	rn = 0.;

	if (NULL != mwt) sf_seek(mwt,0,SEEK_SET);

	grsf = sf_input(g); 
	
	MLOOP( sf_floatread(buf,mbuf,grsf);

	       if (NULL != mwt) { 
		   sf_floatread(wht,mbuf,mwt);
		   for (i=0; i < mbuf; i++) { buf[i] *= wht[i]; }
	       }

	       for (i=0; i < mbuf; i++) { rn += (double) buf[i] * buf[i]; }
	    );
	    
	sfile = fopen(s,"r+b");

	if (0==iter) {
	    alpha = 0.;
	} else {
	    alpha = rn/rnp;
	    
	    if (0 > fseeko(sfile,0,SEEK_SET))
		sf_error ("seek problem");
	}

	rnp = rn;
	    
	/* s = g + alpha*s */

	if (NULL != mwt) sf_seek(mwt,0,SEEK_SET);

	sf_seek(grsf,0,SEEK_SET);
	    
	MLOOP( sf_floatread(buf,mbuf,grsf);
		   
	       if (NULL != mwt) { 
		   sf_floatread(wht,mbuf,mwt);
		   for (i=0; i < mbuf; i++) { buf[i] *= wht[i]; }
	       }

	       if (iter > 0) {
		   pos = ftello(sfile);

		   MREAD2(sfile); 
		       
		   if (0 > fseeko(sfile,pos,SEEK_SET))
		       sf_error ("seek problem");

		   cblas_saxpy(mbuf,alpha,buf2,1,buf,1);
	       } 
		   
	       MWRITE(sfile); );
	
	fclose(sfile);
	sf_fileclose(grsf);

	sf_warning("grad=%lg",rn);

	/* reads from g, runs the forward, and writes to G */

	cmdline[cmdlen-2]='0';
	snprintf(iostring,SF_CMDLEN,"--input=%s --output=%s",g,G);
	iolen = strlen(iostring);
	if (cmdlen+iolen > SF_CMDLEN-1) sf_error("command line is too long");
	snprintf(cmdline+cmdlen,iolen+1,"%s",iostring);

	sf_warning(cmdline);
	sf_system(cmdline);
	
	/* reads G, does computations with it */

	Sfile = fopen(S,"r+b"); if (NULL == Sfile) sf_error("Cannot open %s:",S);

	/* beta = ddot(ss) */

	beta = 0.;

	/* ss = gg + alpha * ss */

	Grsf = sf_input(G); 

	DLOOP( sf_floatread(buf,dbuf,Grsf);

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

	sf_fileclose(Grsf);

	alpha = - rn/beta;

	/* updates x and r */
	
	sfile = fopen(s,"rb");  if (NULL == sfile) sf_error("Cannot open %s:",s);
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

	Sfile = fopen(S,"rb");  if (NULL == Sfile) sf_error("Cannot open %s:",S);

	sf_seek(Rrsf,0,SEEK_SET);

	DLOOP(
	    pos = sf_tell(Rrsf);
	    
	    sf_floatread(buf,dbuf,Rrsf);
	    
	    DREAD2 (Sfile);
	       
	    cblas_saxpy(dbuf,alpha,buf2,1,buf,1);
	       
	    sf_seek(Rrsf,pos,SEEK_SET);
	    
	    sf_floatwrite(buf,dbuf,Rrsf);
	    );

	fclose(Sfile);
    }
	
    /* write x to out */  

    xfile = fopen(x,"rb");
    MLOOP( MREAD(xfile); sf_floatwrite(buf,mbuf,out); );
    fclose(xfile);

    unlink(R);
    unlink(x);
    unlink(g);
    unlink(G);
    unlink(s);
    unlink(S);

    exit(0);
}
