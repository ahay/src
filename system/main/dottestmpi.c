/* Generic dot-product test for linear operators with adjoints 

In this version, the linear operator program uses --input and --output instead of stdin and stdout.
*/
/*
  Copyright (C) 2004 University of Texas at Austin

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
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <time.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int i, im, id;
    unsigned long mseed, dseed;
    off_t nm, nd, msiz, dsiz;
    size_t nbuf, mbuf, dbuf, len, iolen, cmdlen;
    float *buf;
    double dp;
    sf_file mod, mrsf;
    sf_file dat, drsf;
    FILE *mfile, *dfile;
    char *cmdline, *iostring, *arg, *m, *d;

    sf_init(argc,argv);

    mod = sf_input("mod");
    dat = sf_input("dat");

    if (SF_FLOAT != sf_gettype(mod) ||
	SF_FLOAT != sf_gettype(dat))
	sf_error("Need float type in mod and dat");

    nm = sf_filesize(mod);
    nd = sf_filesize(dat);

    nbuf = BUFSIZ/sizeof(float);
    buf = sf_floatalloc(nbuf);

    mseed = (unsigned long) time(NULL);
    init_genrand(mseed);
    mseed = genrand_int32();
    dseed = genrand_int32();

    for (i=0; i < argc-1; i++) {
	argv[i]=argv[i+1];
    }
    argv[argc-1] = sf_charalloc(6);
    snprintf(argv[argc-1],6,"adj=X");

    cmdline = sf_charalloc(SF_CMDLEN);
    iostring = sf_charalloc(SF_CMDLEN);
    cmdlen = 0;
    for (i=0; i < argc; i++) {
	arg = argv[i];
	len = strlen(arg);
	if (cmdlen+len > SF_CMDLEN-2) sf_error("command line is too long");

	strncpy(cmdline+cmdlen,arg,len);
	cmdline[cmdlen+len]=' ';
	cmdlen += len+1;
    }

    mfile = sf_tempfile(&m,"w+b"); 
    fclose(mfile); 
    dfile = sf_tempfile(&d,"w+"); 
    fclose(dfile); 
    
    mrsf = sf_output(m);
    sf_fileflush(mrsf,mod);

    init_genrand(mseed);
    for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {
	if (msiz < mbuf) mbuf=msiz;

	sf_random(mbuf,buf);

	sf_floatwrite(buf,mbuf,mrsf);
    }
    sf_fileclose(mrsf);


    cmdline[cmdlen-2]='0';
    snprintf(iostring,SF_CMDLEN,"--input=%s --output=%s",m,d);
    iolen = strlen(iostring);
    if (cmdlen+iolen > SF_CMDLEN-1) sf_error("command line is too long");
    snprintf(cmdline+cmdlen,iolen+1,"%s",iostring);
    
    sf_warning(cmdline);
    sf_system(cmdline);

    drsf = sf_input(d);
    init_genrand(dseed);
    dp = 0.;
    for (dsiz=nd, dbuf=nbuf; dsiz > 0; dsiz -= dbuf) {
	if (dsiz < dbuf) dbuf=dsiz;
	
	sf_floatread(buf,dbuf,drsf);
	for (id=0; id < dbuf; id++) {
	    dp += buf[id]*genrand_real1 ();
	}	
    }
    sf_fileclose(drsf);

    sf_warning(" L[m]*d=%g",dp);

    drsf = sf_output(d);
    sf_fileflush(drsf,dat);
    
    init_genrand(dseed);
    for (dsiz=nd, dbuf=nbuf; dsiz > 0; dsiz -= dbuf) {
	if (dsiz < dbuf) dbuf=dsiz;

	sf_random(dbuf,buf);

	sf_floatwrite(buf,dbuf,drsf);
    }
    sf_fileclose(drsf);
 
    cmdline[cmdlen-2]='1';
    snprintf(iostring,SF_CMDLEN,"--input=%s --output=%s",d,m);
    iolen = strlen(iostring);
    if (cmdlen+iolen > SF_CMDLEN-1) sf_error("command line is too long");
    snprintf(cmdline+cmdlen,iolen+1,"%s",iostring);
    
    sf_warning(cmdline);
    sf_system(cmdline);

    mrsf = sf_input(m);
    init_genrand(mseed);
    dp = 0.;
    for (msiz=nm, mbuf=nbuf; msiz > 0; msiz -= mbuf) {
	if (msiz < mbuf) mbuf=msiz;
	
	sf_floatread(buf,mbuf,mrsf);
	
	for (im=0; im < mbuf; im++) {
	    dp += buf[im]*genrand_real1 ();
	}	
    }
    sf_fileclose(mrsf);

    sf_warning("L'[d]*m=%g",dp);
	
    exit(0);
}
