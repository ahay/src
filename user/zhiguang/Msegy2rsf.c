/* Convert SEGY to RSF */
/*
 Copyright (C) 2014 University of Texas at Austin
 
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

int main(int argc, char *argv[])
{
	int i,j,ir, it, nshot,ntrace,nt, nh, nr;
	int start, max, ep, sx, gx, shotnum, dsint, drint, offset;
	float ds, dr;
	sf_axis at;

	int **hh, *number;
	float **dd;
	sf_file head, input, output;

	sf_init(argc, argv);

	head=sf_input("head"); /* header file */
	input=sf_input("in");
	output=sf_output("out");

	if(!sf_getint("ep", &ep)) ep=4; /* keyword shot, starting from 0 */
	if(!sf_getint("sx", &sx)) sx=21; /* keyword source position */
	if(!sf_getint("gx", &gx)) gx=23; /* keyword receiver position */
	if(!sf_getint("nshot", &nshot)) sf_error("Need to know shot numbers!");

	if(!sf_getfloat("ds", &ds)) sf_error("Need input source interval!");
	/* shot spacing */
	if(!sf_getfloat("dr", &dr)) sf_error("Need input receiver interval!");
	/* receiver spacing */
	dsint=ds*100;
	drint=dr*100;
	sf_warning("ds=%g dr=%g dsint=%d drint=%d", ds, dr, dsint, drint);

	if(!sf_histint(head, "n1", &nh)) sf_error("No n1= in header file");
	if(!sf_histint(head, "n2", &ntrace)) sf_error("No n2= in header file");
	at=sf_iaxa(input, 1); nt=sf_n(at);

	hh=sf_intalloc2(nh, ntrace);
	sf_intread(hh[0], nh*ntrace, head);
	number=sf_intalloc(nshot);
	for(i=0; i<nshot; i++)
		number[i]=0;

	/* get trace number for each shot and the maximum offset */
	max=abs(hh[0][gx]-hh[0][sx]);
	shotnum=hh[0][ep];
	j=0;
	for(i=0; i<ntrace; i++){
		if(hh[i][ep] == shotnum){
			offset=abs(hh[i][gx]-hh[i][sx]);
			if(max<offset)
				max=offset;
			number[j]++;
		}else{
			j++;
			shotnum=hh[i][ep];
			number[j]++;
			offset=abs(hh[i][gx]-hh[i][sx]);
			if(max<offset)
				max=offset;
		}
	}
	sf_warning("The maximum offset is %fm", max/100.);

	/* check the trace numbers at each shot */
	shotnum=0;
	for(i=0; i<nshot; i++){
		sf_warning("Shot #%d, Trace number: %d",i+1, number[i]);
		shotnum+=number[i];
	}
	if(shotnum != ntrace) sf_error("Calculation error and check the program please");

	max=max/drint;

	nr=2*max+1;
	sf_warning("The maximum trace number is %d", nr);

	sf_oaxa(output, at, 1);
	sf_putint(output, "n2", nr);
	sf_putfloat(output, "d2", dr);
	sf_putfloat(output, "o2", -max*dr);
	sf_putstring(output, "label2", "Offset");
	sf_putstring(output, "unit2", "m");

	sf_putint(output, "n3", nshot);
	sf_putfloat(output, "d3", ds);
	sf_putfloat(output, "o3", hh[0][sx]/100.);
	sf_putstring(output, "label3", "Shot");
	sf_putstring(output, "unit3", "m");

	dd=sf_floatalloc2(nt, nr);

	j=0;
	for(i=0; i<nshot; i++){
		for(ir=0; ir<nr; ir++)
			for(it=0; it<nt; it++)
				dd[ir][it]=0.;
		start=(hh[j][gx]-hh[j][sx])/drint+max;
		sf_floatread(dd[0]+start*nt, number[i]*nt, input);
		sf_floatwrite(dd[0], nr*nt, output);
		j+=number[i];
	}

	/* Another check */
	if(j != ntrace) sf_error("Calculation error and check the program please");

	sf_fileclose(head);
	sf_fileclose(input);
	sf_fileclose(output);

	return 0;
}
