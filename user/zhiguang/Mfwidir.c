/* Update the conjugate direction in full waveform inversion */
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

float direction_cg_polak(float **g0,
                         float **d0,
                         float **g1,
                         float **d1,
                         const int n1, 
                         const int n2)
/*< Polak and Ribiere 1969 >*/ 
{
    int i,j; 
    float diff, beta_top, beta_bot, beta;

    beta_top=0.0;
    beta_bot=0.0;

    for ( j = 0; j < n2; j++ ) {
    for ( i = 0; i < n1; i++ ) {
        diff = g1[j][i] - g0[j][i];
        beta_top += g1[j][i] * diff;
        beta_bot += g0[j][i] * g0[j][i];
    }
    }

    beta=beta_top / beta_bot;
    
    if (beta < 0.0 ) { 
        beta =0.0; 
    }

    for ( j = 0; j < n2; j++ ){
    for ( i = 0; i < n1; i++ ) { 
        d1[j][i] = -g1[j][i] + beta * d0[j][i];
    }
    }
    return beta; 
}

float direction_cg_fletcher(float **g0,
                            float **d0,
                            float **g1,
                            float **d1,
                            const int n1, 
                            const int n2)
/*< Fletcher and Reeves 1964 >*/ 
{
    int i,j; 
    float beta_top, beta_bot, beta;

    beta_top=0.0;
    beta_bot=0.0;
    for ( j = 0; j < n2; j++ ) {
    for ( i = 0; i < n1; i++ ) {
        beta_top += g1[j][i] * g1[j][i];
        beta_bot += g0[j][i] * g0[j][i];
    }
    }

    beta= beta_top / beta_bot;
    
    if (beta < 0.0 ) { 
        beta =0.0; 
    }

    for ( j = 0; j< n2; j++ ){
    for ( i = 0; i < n1; i++ ) { 
        d1[j][i] = g1[j][i] + beta * d0[j][i];
    }
    }
    return beta; 
}

float direction_cg_fletcher2(float **g0,
                            float **d0,
                            float **g1,
                            float **d1,
                            const int n1, 
                            const int n2)
/*< Fletcher 1987 >*/ 
{
    int i,j; 
    float beta_top, beta_bot, beta;

    beta_top=0.0;
    beta_bot=0.0;
    for ( j = 0; j < n2; j++ ) {
    for ( i = 0; i < n1; i++ ) {
        beta_top += g1[j][i] * g1[j][i];
        beta_bot += d0[j][i] * g0[j][i];
    }
    }

    beta= - beta_top / beta_bot;
    
    if (beta < 0.0 ) { 
        beta =0.0; 
    }

    for ( j = 0; j< n2; j++ ){
    for ( i = 0; i < n1; i++ ) { 
        d1[j][i] = g1[j][i] + beta * d0[j][i];
    }
    }
    return beta; 
}

float direction_cg_dai(float **g0,
                         float **d0,
                         float **g1,
                         float **d1,
                         const int n1, 
                         const int n2)
/*< Dai and Yuan 1999 >*/ 
{
    int i,j; 
    float diff, beta_top, beta_bot, beta;

    beta_top=0.0;
    beta_bot=0.0;

    for ( j = 0; j < n2; j++ ) {
    for ( i = 0; i < n1; i++ ) {
        diff = g1[j][i] - g0[j][i];
        beta_top += g1[j][i] * g1[j][i];
        beta_bot += d0[j][i] * diff;
    }
    }

    beta=beta_top / beta_bot;
    
    if (beta < 0.0 ) { 
        beta =0.0; 
    }

    for ( j = 0; j < n2; j++ ){
    for ( i = 0; i < n1; i++ ) { 
        d1[j][i] = -g1[j][i] + beta * d0[j][i];
    }
    }
    return beta; 
}

int main(int argc, char* argv[])
{
	int n1, n2;
	float beta0;
	char *option;
	float **g0, **g1, **d0, **d1;
	sf_file in, out, grad0, grad1, beta;

	sf_init(argc, argv);

	in=sf_input("in");
	/* old search direction */
	out=sf_output("out");
	/* new search direction */
	grad0=sf_input("grad0");
	/* old gradient */
	grad1=sf_input("grad1");
	/* new gradient */

	option=sf_getstring("option");
	/* CG update parameter (Polak, Fletcher-Reeves, Fletcher, Dai-Yuan) */
	if(NULL == option) sf_error("Need option=");

	if(!sf_histint(in, "n1", &n1)) sf_error("No n1= in input.");
	if(!sf_histint(in, "n2", &n2)) sf_error("No n2= in input.");

	g0=sf_floatalloc2(n1, n2);
	g1=sf_floatalloc2(n1, n2);
	d0=sf_floatalloc2(n1, n2);
	d1=sf_floatalloc2(n1, n2);

	sf_floatread(g0[0], n1*n2, grad0);
	sf_floatread(g1[0], n1*n2, grad1);
	sf_floatread(d0[0], n1*n2, in);

	switch(option[0]){
		case 'p':
			beta0=direction_cg_polak(g0, d0, g1, d1, n1, n2);
			break;
		case 'r':
			beta0=direction_cg_fletcher(g0, d0, g1, d1, n1, n2);
			break;
		case 'f':
			beta0=direction_cg_fletcher2(g0, d0, g1, d1, n1, n2);
			break;
		case 'd':
			beta0=direction_cg_dai(g0, d0, g1, d1, n1, n2);
			break;
	}

	sf_floatwrite(d1[0], n1*n2, out);
	if(NULL != sf_getstring("beta")) {
		beta=sf_output("beta");
		sf_putint(beta,"n1",1);
		sf_putint(beta,"n2",1);
		sf_floatwrite(&beta0, 1, beta);
	}

	exit(0);
}

