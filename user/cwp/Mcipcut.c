/* cut at CIPs */

/*
  Copyright (C) 2011 Colorado School of Mines

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

#define BOUND(i,n) (i<0) ? 0 : ( (i>n-1) ? n-1 : i );

int main(int argc, char* argv[])
{
    bool verb;

    sf_file Fcub=NULL; /* CUB */
    sf_file Fcip=NULL; /* CIP */
    sf_file Fcut=NULL; /* CUT */

    sf_axis ax,ay,az,ac,aj;
    int     ix,iy,iz,ic;
    float   fx,fy,fz;
    float w000,w001,w010,w011,w100,w101,w110,w111;


    /*  arrays  */
    float ***cub;
    vc3d    *cip;
    float   *cut;

    /*-----------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);    
    
    if(! sf_getbool("verb",&verb)) verb=false;	/* verbosity flag */

    Fcub=sf_input ( "in");
    Fcip=sf_input ("cip");
    Fcut=sf_output("out");

    /* input axes */
    ax = sf_iaxa(Fcub,1); sf_setlabel(ax,"x");
    ay = sf_iaxa(Fcub,2); sf_setlabel(ay,"y");
    az = sf_iaxa(Fcub,3); sf_setlabel(az,"z");
    ac = sf_iaxa(Fcip,2); sf_setlabel(ac,"c");
		
    aj  = sf_maxa(1,0,1);

    /* output axes */
    sf_oaxa(Fcut,ac,1);
    sf_oaxa(Fcut,aj,2);
    sf_oaxa(Fcut,aj,3);
    sf_oaxa(Fcut,aj,4);
    sf_oaxa(Fcut,aj,5);

    if (verb){
	sf_raxa(ax);
	sf_raxa(ay);
	sf_raxa(az);
	sf_raxa(ac);
    }

    /*------------------------------------------------------------*/
    /* read cube */
    cub = sf_floatalloc3(sf_n(ax),sf_n(ay),sf_n(az));
    sf_floatread(cub[0][0],sf_n(ax)*sf_n(ay)*sf_n(az),Fcub);

    /*------------------------------------------------------------*/
    /* read cips */
    cip = (vc3d*) sf_alloc(sf_n(ac),sizeof(*cip));
    vc3dread1(Fcip,cip,sf_n(ac));

    cut = sf_floatalloc(sf_n(ac));

    /*------------------------------------------------------------*/
    /* loop over CIPs */
    if(verb) fprintf(stderr,"ic\n");
    for(ic=0;ic<sf_n(ac);ic++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",ic);

	ix=BOUND((cip[ic].dx-sf_o(ax))/sf_d(ax),sf_n(ax));
	iy=BOUND((cip[ic].dy-sf_o(ay))/sf_d(ay),sf_n(ay));	
	iz=BOUND((cip[ic].dz-sf_o(az))/sf_d(az),sf_n(az));

	fx=(cip[ic].dx-sf_o(ax))/sf_d(ax)-ix;
	fy=(cip[ic].dy-sf_o(ay))/sf_d(ay)-iy;
	fz=(cip[ic].dz-sf_o(az))/sf_d(az)-iz;

        w000 = (1-fz)*(1-fy)*(1-fx);
        w001 = (1-fz)*(  fy)*(1-fx);
        w010 = (1-fz)*(1-fy)*(  fx);
        w011 = (1-fz)*(  fy)*(  fx);
        w100 = (  fz)*(1-fy)*(1-fx);
        w101 = (  fz)*(  fy)*(1-fx);
        w110 = (  fz)*(1-fy)*(  fx);
        w111 = (  fz)*(  fy)*(  fx);

	cut[ic]=0;
	cut[ic]+=w000*cub[iz  ][iy  ][ix  ];
	cut[ic]+=w001*cub[iz  ][iy+1][ix  ];
	cut[ic]+=w010*cub[iz  ][iy  ][ix+1];
	cut[ic]+=w011*cub[iz  ][iy+1][ix+1];
	cut[ic]+=w100*cub[iz+1][iy  ][ix  ];
	cut[ic]+=w101*cub[iz+1][iy+1][ix  ];
	cut[ic]+=w110*cub[iz+1][iy  ][ix+1];
	cut[ic]+=w111*cub[iz+1][iy+1][ix+1];

    }
    if(verb) fprintf(stderr,"\n");
    /*------------------------------------------------------------*/
    sf_floatwrite(cut,sf_n(ac),Fcut);


    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"free memory...");
    free(**cub);free(*cub);free(cub);
    ;                      free(cut);
    ;                      free(cip);
    if(verb) fprintf(stderr,"OK\n");
    /*------------------------------------------------------------*/

    exit(0);
}		
