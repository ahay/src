/* Oriented zero-offset migration. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
#include <assert.h>

float ** t0, *m0, *minit, *m;

int   read4file_ref(char *fname,  int nx, int nz, float **s)
{
    int i, NM;
    FILE *fp;

    if((fp=fopen(fname, "rb"))==NULL) {
	printf("Cannot open file.\n");
    }

    /* slowness */
    if(fread(s[0], sizeof(float), nx*nz, fp) != nx*nz) {
	if (feof(fp))
	    printf("File read error - too small.");
    }

    /* nm */
    if(fread(&i, sizeof(int), 1, fp) != 1) {
	if (feof(fp))
	    printf("File read error - nm.");
    }

    NM = i;

    m0    = sf_floatalloc(NM);
    minit = sf_floatalloc(NM);
    m     = sf_floatalloc(NM);

    /* m0 */
    if(fread(m0, sizeof(float), NM, fp) != NM) {
	if (feof(fp))
	    printf("File read error - too small.");
    }
    if(fread(minit, sizeof(float), NM, fp) != NM) {
	if (feof(fp))
	    printf("File read error - too small.");
    }
    if(fread(m, sizeof(float), NM, fp) != NM) {
	if (feof(fp))
	    printf("File read error - too small.");
    }

    fclose(fp);

    return NM;
}
/*************************************************/

void putf(sf_file so1, int nx, int nz, float dx, float dz)
{
    sf_putint (so1, "n3", 1);
    sf_putint (so1, "n2", nx);
    sf_putint (so1, "n1", nz);
    sf_putfloat (so1, "d3", 0);
    sf_putfloat (so1, "d2", dx);
    sf_putfloat (so1, "d1", dz);
    sf_putfloat (so1, "o3", 0);
    sf_putfloat (so1, "o2", 0);
    sf_putfloat (so1, "o1", 0);
}
int main(int argc, char* argv[])
{
    int nm = 1e6;//1234567890;//1e10
    int n, nx, nz;
    float dx, dz;
    sf_file so, so1, so2, so3;
    char * fname = 0;

    sf_init(argc,argv);

    so = sf_output("out");

    so1 = sf_output("correct");
    so2 = sf_output("init");
    so3 = sf_output("final");

    if (!sf_getint("N",&n)) sf_error("No N= ");
    nx = nz = n;

    dx = dz = 1.f / (n - 1);

    fname = sf_getstring ("sample");

    t0    = sf_floatalloc2(nz,nx);

    // sprintf(fname,"sample%-3d",S.nx);
    nm=read4file_ref(fname, nx, nz, t0);

    putf(so, nx, nz, dx, dz);

    putf(so1, 1, nm, dx, dz);
    putf(so2, 1, nm, dx, dz);
    putf(so3, 1, nm, dx, dz);

    sf_floatwrite(t0[0],    nx*nz, so);

    sf_floatwrite(m0, nm, so1);
    sf_floatwrite(minit, nm, so2);
    sf_floatwrite(m,     nm, so3);

    sf_close();

    free(m0);
    free(m);
    free(minit);
    exit(0);
}
