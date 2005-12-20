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
#include <stdio.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int i, iter, niter, status, p[4][2];
    char *name;
    float *buf;
    double rn;
    pid_t pid[6];
    off_t nm, nd, msiz, dsiz;
    size_t nbuf, mbuf, dbuf;
    FILE *x, *rr, *s, *ss;
    sf_file mod, dat, out, from, to;

    sf_init(argc,argv);
    dat = sf_input("in");
    mod = sf_input("mod");

    if (SF_FLOAT != sf_gettype(mod) ||
	SF_FLOAT != sf_gettype(dat)) 
	sf_error("Need float type in mod and dat");
  
    free(argv[0]);
    for (i=0; i < argc-1; i++) {
	argv[i]=argv[i+1];
    }
    argv[argc-1] = sf_charalloc(6);
    snprintf(argv[argc-1],6,"adj=X");

    if (!sf_getint("niter",&niter)) niter=1;
    /* number of iterations */

    rr = sf_tempfile(&name,"w+b");
    x  = sf_tempfile(&name,"w+b");

    nm = sf_filesize(mod);
    nd = sf_filesize(dat);

    nbuf = BUFSIZ/sizeof(float);
    buf = sf_floatalloc(nbuf);

    /* r=-dat */

    for (dsiz=nd, dbuf=nbuf; dsiz > 0; dsiz -= dbuf) {
	if (dsiz < dbuf) dbuf=dsiz;

	sf_floatread(buf,dbuf,dat);
	for (i=0; i < dbuf; i++) {
	    buf[i]=-buf[i];
	}

	fwrite(buf,sizeof(float),dbuf,rr);
    }
    
    /* x=0 */

    for (i=0; i < nbuf; i++) {
	buf[i]=0.;
    }
    
    for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {
	if (msiz < mbuf) mbuf=msiz;
	fwrite(buf,sizeof(float),mbuf,x);
    }

    for (i=0; i < 4; i++) { /* make four pipes */
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
	    dup(p[0][1]);

	    to = sf_output("out");

	    rewind(rr);
	    for (dsiz=nd, dbuf=nbuf; dsiz > 0; dsiz -= dbuf) {
		if (dsiz < dbuf) dbuf=dsiz;
		fread(buf,sizeof(float),dbuf,rr);
		sf_floatwrite(buf,dbuf,to);
	    }

	    _exit(0);
	}
 
	
	if (0==pid[1]) {
	    /* reads from p[0], runs the adjoint, and writes to p[1] */

	    close(p[0][1]);
	    close(STDIN_FILENO);
	    dup(p[0][0]);
	    
	    close(p[1][0]);
	    close(STDOUT_FILENO);
	    dup(p[1][1]);
	    
	    argv[argc-1][4]='1';
	    execvp(argv[0],argv);
	    _exit(1);
	}


	if (0==pid[2]) {
	    /* reads gradient from p[1], uses it, and writes to p[2] */
	
	    close(p[1][1]);
	    close(STDIN_FILENO);
	    dup(p[1][0]);
	    from = sf_input("in");
	    
	    close(p[2][0]);
	    close(STDOUT_FILENO);
	    dup(p[2][1]);
	    to = sf_output("out");

	    rn = 0.;
	    for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {
		if (msiz < mbuf) mbuf=msiz;
		sf_floatread(buf,mbuf,from);
		sf_floatwrite(buf,mbuf,to);
	    }
	    _exit(2);
	}

	if (0==pid[3]) {
	    /* reads from p[2], runs the forward, and writes to p[3] */

	    close(p[2][1]);
	    close(STDIN_FILENO);
	    dup(p[2][0]);
	    
	    close(p[3][0]);
	    close(STDOUT_FILENO);
	    dup(p[3][1]);
	    
	    argv[argc-1][4]='0';
	    execvp(argv[0],argv);
	    _exit(3);
	}

	if (0==pid[4]) {
	    /* reads conj. grad from p[3], does computations with it */

	    close(p[3][1]);
	    close(STDIN_FILENO);
	    dup(p[3][0]);
	    from = sf_input("in");
	    
	    for (dsiz=nd, dbuf=nbuf; dsiz > 0; dsiz -= dbuf) {
		if (dsiz < dbuf) dbuf=dsiz;
		sf_floatread(buf,dbuf,from);
	    }

	}

	for (i=0; i < 6; i++) { 
	    if (0 == pid[i]) break;
	}    
	if (6==i) {
	    /* parent waits */
	    waitpid(pid[2],&status,0);
	    waitpid(pid[5],&status,0); 
	} 
	    
    }

    /* write x to out */  
    out = sf_output("out");
    sf_fileflush(out,mod);
  
    rewind(x);
    for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {
	if (msiz < mbuf) mbuf=msiz;
	fread(buf,sizeof(float),mbuf,x);
	sf_floatwrite(buf,mbuf,out);
    }

    exit(0);
}



#ifdef asdfcgf


	/* child1 feeds rr to A' */

	oper( T, F, g, rr);
	   
	/* child2 reads g, computes rn, alpha, and s and feeds g to A */
        
	oper( F, F, g, gg);       

	/* child 3 reads gg, computes ss and beta */

	/* parent: wait for children, update x and rr */

 
	rn = sum( dprod( g, g));

	if(iter=0) {
	    s  = g;
	    ss = gg;		
	} else {
	    alpha = rn / rnp;
	    s  =  g + alpha *  s;
	    ss = gg + alpha * ss;
	}

    beta = sum( dprod( ss, ss))
    if( beta > epsilon( beta)) {
           alpha = - rn / beta
           x =  x  + alpha *  s
           rr = rr + alpha * ss
           }
    rnp = rn;    forget = .false.;   stat = 0

#endif



