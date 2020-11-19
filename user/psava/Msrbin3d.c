/* 4-D data binning from traces at irregular coordinates */
/*
  Copyright (C) 2012 Colorado School of Mines
  
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

#include <rsf.h>

int main (int argc, char* argv[])
{
    bool verb;
    sf_axis a1,a2,b1,b2;

    float      o1,o2;
    float      d1,d2;
    int        n1,n2;
    int        i1,i2;

    sf_axis ak,ae,at;
    int     ik,ne,ie,it;
    sf_file Ftrc=NULL; /* data (traces) */
    sf_file Fbin=NULL; /* binned data   */
    sf_file Fkey=NULL; /* mapping keys  */

    int** tmp;
    int*  key;
    float*trc;

    int *e1, *e2;
    int j1=0,j2=0;
    int k1,k2;
    off_t iseek,start;

    /*------------------------------------------------------------*/
    sf_init (argc,argv);
    /*------------------------------------------------------------*/

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */

    /* I/O files */
    Ftrc = sf_input ("in");
    at = sf_iaxa(Ftrc,1); sf_setlabel(at,"t"); /* time axis */
    if(verb) sf_raxa(at);
    trc=sf_floatalloc(sf_n(at));

    Fbin = sf_output("out");
    sf_settype(Fbin,SF_FLOAT);

    Fkey = sf_input ("key");
    if (SF_INT != sf_gettype(Fkey)) sf_error("Need int keys");
    ak = sf_iaxa(Fkey,2); sf_setlabel(ak,"k"); sf_setunit(ak,"");/* key axis */
    if(verb) sf_raxa(ak);

    /*------------------------------------------------------------*/
    /* experiment space */
    if(! sf_getint  ("n1",&n1)) sf_error("need n1");
    if(! sf_getfloat("o1",&o1)) sf_error("need o1");
    if(! sf_getfloat("d1",&d1)) sf_error("need d1");
    a1 = sf_maxa(n1,o1,d1);

    if(! sf_getint  ("n2",&n2)) sf_error("need n2");
    if(! sf_getfloat("o2",&o2)) sf_error("need o2");
    if(! sf_getfloat("d2",&d2)) sf_error("need d2");
    a2 = sf_maxa(n2,o2,d2);

    tmp=sf_intalloc2(sf_n(a1),sf_n(a2));
    for    (i2=0; i2<sf_n(a2); i2++) {
	    for(i1=0; i1<sf_n(a1); i1++) {
		    tmp[i2][i1]=0;
	    }
    }
    /* go through all keys... */
    key=sf_intalloc(4);
    for(ik=0; ik<sf_n(ak); ik++) {
	    sf_intread(key,4,Fkey);
	    tmp[ key[1] ][ key[0] ]=1;
    }
    sf_seek(Fkey,0,SEEK_SET); /* seek to start */

    /* ... and count the number of unique experiments ... */
    ne=0;
    for    (i2=0; i2<sf_n(a2); i2++) {
	    for(i1=0; i1<sf_n(a1); i1++) {
		    ne+=tmp[i2][i1];
	    }
    }
    sf_warning("%d experiment(s)",ne);

    /* ... and keep the experiment indices */
    e1=sf_intalloc(ne);
    e2=sf_intalloc(ne);
    ie=0;
    for    (i2=0; i2<sf_n(a2); i2++) {
	    for(i1=0; i1<sf_n(a1); i1++) {
		    if(tmp[i2][i1]==1) {
			    e1[ie]=i1;
			    e2[ie]=i2;
			    ie++;
		    }
	    }
    }
    free(*tmp); free(tmp);

    /*------------------------------------------------------------*/

    /* make experiment axis */
    ae=sf_maxa(ne,0,1); 
    sf_setlabel(ae,"e");
    sf_setunit (ae,"");
    if(verb) sf_raxa(ae);

    /* input offset axes */
    if(! sf_getint  ("on1",&n1)) sf_error("need on1");
    if(! sf_getfloat("oo1",&o1)) sf_error("need oo1");
    if(! sf_getfloat("od1",&d1)) sf_error("need od1");
    b1 = sf_maxa(n1,o1,d1);
    k1 = o1/d1;

    if(! sf_getint  ("on2",&n2)) sf_error("need on2");
    if(! sf_getfloat("oo2",&o2)) sf_error("need oo2");
    if(! sf_getfloat("od2",&d2)) sf_error("need od2");
    b2 = sf_maxa(n2,o2,d2);
    k2 = o2/d2;

    /* output header */
    sf_oaxa(Fbin,at,1);
    sf_oaxa(Fbin,b1,2);
    sf_oaxa(Fbin,b2,3);
    sf_oaxa(Fbin,ae,4);

    /*------------------------------------------------------------*/

    /* create zero output array */
    for(it=0;it<sf_n(at);it++) {
	    trc[it]=0;
    }
    for(ie=0;ie<ne;ie++) {
	    for    (i2=0;i2<sf_n(b2);i2++) {
		    for(i1=0;i1<sf_n(b1);i1++) {
			    sf_floatwrite(trc,sf_n(at),Fbin);
		    }
	    }
    }
    sf_seek(Fbin,0,SEEK_SET);
    start = sf_tell(Fbin);

    /*------------------------------------------------------------*/
    /* loop over traces */
    if(verb) fprintf(stderr,"\n");
    for(ik=0; ik<sf_n(ak); ik++) {
	    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b%d",ik);

	    /* read keys for the trace */
	    sf_intread(key,4,Fkey);

	    /* find experiment index */
	    for(ie=0;ie<ne;ie++) {
		    if(key[0]==e1[ie] && key[1]==e2[ie]) {
			    j1=key[2];
			    j2=key[3];
			    break;
		    }
	    }

	    /* seek in the output file */
	    iseek = ie*sf_n(b1)*sf_n(b2)+
		    (j2-k2)*sf_n(b1)+
	    (j1-k1);
	sf_seek(Fbin,start+iseek*sf_n(at)*sizeof(float),SEEK_SET);
	    
	/* read a trace */
	sf_floatread(trc,sf_n(at),Ftrc);
	/* write trace */
	sf_floatwrite(trc,sf_n(at),Fbin);

    }
    if(verb) fprintf(stderr,"\n");

    exit(0);
}


