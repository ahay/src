/* Ps to vplot */
/*
Copyright (C) 2013 Zhonghuan Chen, Tsinghua University

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
#include <rsfplot.h>

#define BUFSZ	100
#define PARSZ	10
#define ARRAYSZ	10

void box(float x1, float x2, float a1, float a2, float th);
void ell(float x1, float x2, float a1, float a2, float th);

int main (int argc, char **argv)
{
	int np, nl;
	char buf[BUFSZ*PARSZ], para[PARSZ][BUFSZ];
	char *sep = " \t", *p;
	float x1, x2;

	sf_init(argc, argv);

    vp_init();

	nl = 0;
	while(fgets(buf, BUFSZ*PARSZ, stdin))
	{
		nl++;
		np = 0;
		for(p=strtok(buf,sep);
			p;
			p=strtok(NULL,sep))
		{
			strcpy(para[np], p);
			if(np++ >= PARSZ) sf_error("L%d: too many parameters", nl);
		}
		if(np<=0) continue;
		if(strcmp(para[0],"dash")==0 && np==2)
			vp_set_dash(atof(para[1]));
		else if(strcmp(para[0],"color")==0 && np==2)
			vp_color(atoi(para[1]));
		else if(strcmp(para[0],"style")==0 && np==2)
			vp_style(atoi(para[1]));
		else if(strcmp(para[0],"fat")==0 && np==2)
			vp_fat(atoi(para[1]));
		else if(strcmp(para[0],"move")==0 && np==3)
		{
			x1 = atof(para[1]);
			x2 = atof(para[2]);
			vp_move(x1, x2);
		}
		else if(strcmp(para[0],"line")==0 && np==3)
			vp_draw(atof(para[1]), atof(para[2]));
		else if(strcmp(para[0],"text")==0 && np==4)
			vp_text(x1, x2, atoi(para[1]), atof(para[2]), para[3]);
		else if(strcmp(para[0],"box")==0 && np==4)
			box(x1, x2, atof(para[1]), atof(para[2]), atof(para[3])*SF_PI/180);
		else if(strcmp(para[0],"ell")==0 && np==4)
			ell(x1, x2, atof(para[1]), atof(para[2]), atof(para[3])*SF_PI/180);
		else if(strcmp(para[0],"arrow")==0 && np==4)
			vp_arrow(x1, x2, atof(para[1]), atof(para[2]), atof(para[3]));
		else sf_error("L%d: unknown cmd %s", nl, para[0]);
	}
    return 0;
}

void rotate(float *r1, float *r2, int np, float theta)
{
	float r, sinth, costh;
	int i;
	sinth = sin(theta);
	costh = cos(theta);
	for(i=0; i<np; i++)
	{
		r = r1[i]*costh+r2[i]*sinth;
		r2[i] = r2[i]*costh+r1[i]*sinth;
		r1[i] = r;
	}
}


void box(float x1, float x2, float a1, float a2, float theta)
{
	float r1[5], r2[5];
	int i;
	r1[0] = a1/2;
	r2[0] = a2/2;
	r1[1] = -r1[0];
	r2[1] = r2[0];
	r1[2] = -r1[0];
	r2[2] = -r2[0];
	r1[3] = r1[0];
	r2[3] = -r2[0];
	r1[4] = r1[0];
	r2[4] = r2[0];

	rotate(r1, r2, 5, theta);
	for(i=0; i<5; i++)
	{	r1[i]+=x1; 	r2[i]+=x2;	}

	vp_upendn(r1[0], r2[0]);
	vp_pline(r1,r2,5);
	vp_penup();
}


void ell(float x1, float x2, float a1, float a2, float theta)
{
#define ELLNUM	200
	int i;
	float r1[ELLNUM], r2[ELLNUM];
	
	for(i=0; i<ELLNUM; i++)
	{
		r1[i] = a1*cos(2.0*SF_PI*i/(ELLNUM-1));
		r2[i] = a2*sin(2.0*SF_PI*i/(ELLNUM-1));
	}

	rotate(r1, r2, ELLNUM, theta);
	for(i=0; i<ELLNUM; i++)
	{	r1[i]+=x1; 	r2[i]+=x2;	}

	vp_upendn(r1[0], r2[0]);
	vp_pline(r1,r2,ELLNUM);
	vp_penup();
}

