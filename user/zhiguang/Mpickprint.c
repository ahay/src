/* Write predictive painting result into a txt file */
/*
  Copyright (C) 2016 University of Texas at Austin

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

int main(int argc, char* argv[])
{
    bool verb;
    int i, nz, iz, nh, ih, ncdp, icdp, np, ip, jtau, *itau, *mute, *npicks;
    float dz, dh, dcdp, cut, s, z, h;
    float *tau, **ztau, *semb;
    sf_file in, pick, npick, semblance;
    FILE *fp;

    sf_init(argc, argv);
    in = sf_input("in");
    pick = sf_input("pick");
    npick = sf_input("npick");
    semblance = sf_input("semblance");
    fp=fopen("sogpicks.txt","w");

    if (!sf_histint(in,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    ncdp = sf_leftsize(in,2);
    if (!sf_histfloat(in,"d3",&dcdp)) sf_error("No d3= in input");
    if (!sf_histint(pick,"n1",&np)) sf_error("No n1= in pick");

    if (!sf_getbool("verb",&verb)) verb=true;
    /* if y, print icdp/ncdp during operation */
    if (!sf_getfloat("cut",&cut)) cut=0.;
    /* muting value in boundary */

	itau = sf_intalloc(np);
	tau = sf_floatalloc(np);
	mute = sf_intalloc(nz);
	ztau = sf_floatalloc2(nz, nh);
	semb = sf_floatalloc(nz);

	npicks = sf_intalloc(ncdp);
	sf_intread(npicks, ncdp, npick);

    for (icdp=0; icdp < ncdp; icdp++) {

		if(verb) sf_warning("icdp/ncdp=%d/%d;",icdp+1, ncdp);

		/* read data */
		sf_floatread(ztau[0], nz*nh, in);
		sf_floatread(tau, np, pick);
		sf_floatread(semb, nz, semblance);

		for (ip=0; ip<np; ip++){
			itau[ip] = tau[ip]/dz+0.5;
		}

		/* muting boundary */
		for (iz=0; iz<nz; iz++){
			if(ztau[0][iz]>=cut) break;
		}
		for (i=0; i<iz; i++) mute[i]=dz*(nh-1)*(iz-i)/cut+0.5;
		for (i=iz; i<nz; i++) mute[i]=0;

		fprintf(fp, "%4d %11.5f %2d\n", icdp+1, icdp*dcdp, npicks[icdp]);

		for (ip=0; ip<npicks[icdp]; ip++){
			jtau=itau[npicks[icdp]-1-ip];
			fprintf(fp, "%2d %11.5f %3d\n", ip+1, jtau*dz, nh-mute[jtau]);

			s=semb[jtau];
			for (ih=nh-1; ih>=mute[jtau]; ih--){
				h=dh*(nh-1-ih);
				z=ztau[ih][jtau];
				fprintf(fp, "%11.5f %11.5f %6.5f\n", h, z, s);
			} // end of ih
		} // end of ip
	} // end of icdp

	fclose(fp);

    exit(0);
}
