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

int parsplit(char*buf, char**par)
{
	char *p1, *p2;
	int n2=0;
	bool inquote, loop;

	p1 = buf;
	p2 = par[n2];

	inquote=false;
	loop=true;
	while(loop)
	{
		switch(*p1)
		{
		case '"':
		case '\'':
			inquote=!inquote;
			p1++;
			break;
		case ' ':
		case '\t':
			if(inquote==false){
				*p2 = '\0';
				n2++;
				p2 = par[n2];
				p1++;
			}else{
				*p2++=*p1++;
			}
			break;
		case '#':
		case '\0':
			*p2 = '\0';
			n2++;
			loop=false;
			break;
		case '\\':
			p1++;
			*p2++=*p1++;
			break;
		default:
			*p2++=*p1++;
			break;
		} 
	}
	
	return n2;
}

int cmdtree(char **para, int np)
{
	static float x1, x2, w, h;
	static int fontsz=7;
	if(strcmp(para[0],"dash")==0 && np==2)
		vp_set_dash(atof(para[1]));
	else if(strcmp(para[0],"color")==0 && np==2)
		vp_color(atoi(para[1]));
	else if(strcmp(para[0],"style")==0 && np==2)
		vp_style(atoi(para[1]));
	else if(strcmp(para[0],"fat")==0 && np==2)
		vp_fat(atoi(para[1]));
	else if(strcmp(para[0],"fontsz")==0 && np==2)
		fontsz=atoi(para[1]);
	else if(strcmp(para[0],"move")==0 && np==3)
	{
		x1 = atof(para[1]);
		x2 = atof(para[2]);
		vp_move(x1, x2);
	}
	else if(strcmp(para[0],"line")==0 && np==3)
		vp_draw(atof(para[1]), atof(para[2]));
	else if(strcmp(para[0],"text")==0)	
		switch(np)
		{
		case 4:
		vp_text(x1, x2, atoi(para[1]), atof(para[2]), para[3]);
		break;
		case 3:
		vp_text(x1, x2, fontsz, atof(para[1]), para[2]);
		break;
		case 2:
		vp_text(x1, x2, fontsz, 0.0, para[1]);
		break;
		default: break;
		}
	else if(strcmp(para[0],"box")==0)
		switch(np)
		{
		case 4:
		box(x1, x2, atof(para[1]), atof(para[2]), atof(para[3])*SF_PI/180);
		break;
		case 3:
		box(x1, x2, atof(para[1]), atof(para[2]), 0.0);
		break;
		default: break;
		}
	else if(strcmp(para[0],"frame")==0 && np==4 )
	{
		w = atof(para[1]);
		h = atof(para[2]);
		box(x1+w/2, x2+h/2, w, h, 0.0);
		vp_text(x1+0.1, x2+0.1, fontsz, 0.0, para[3]);
	}else if(strcmp(para[0],"ell")==0 && np==4)
		ell(x1, x2, atof(para[1]), atof(para[2]), atof(para[3])*SF_PI/180);
	else if(strcmp(para[0],"arrow")==0 && np==4)
		vp_arrow(x1, x2, atof(para[1]), atof(para[2]), atof(para[3]));
	else return -1;
	return 0;
}

char* cmdline(char *buf, char *p, char**para)
{
	static int nl=0;
	int np=0;
	do{ *p--='\0'; }while(*p==' ' || *p=='\t');
	p=buf;
	while(p[0]==' ' || p[0] == '\t') p++;
/*		for(p=strtok(buf,sep);
			p;
			p=strtok(NULL,sep))
		{
			strcpy(para[np], p);
			if(np++ >= PARSZ) sf_error("L%d: too many parameters", nl);
		}
*/
	nl++;
	if(p[0]=='#' || p[0]=='\0') return buf;
	np = parsplit(p, para);
	if(np <= 0)  return buf;
	if(cmdtree(para, np)<0)
		sf_error("L%d: unknown cmd %s", nl, para[0]);
	return buf;
}

bool issep(char c, char*sep)
{
	char *p;
	if(c=='\n') return true;
	p = sep;
	while(*p!='\0')
	{
		if(c == *p) return true;
		p++;
	}
	return false;
}

int main (int argc, char **argv)
{
	int np;
	char buf[BUFSZ*PARSZ], **para, *p, *sep;

	sf_init(argc, argv);
	if ((sep=sf_getstring("sep"))==NULL) sep="";
	/* cmd separater */

	para = (char **)sf_alloc(PARSZ, sizeof(char*));
	*para = (char *)sf_alloc(PARSZ*BUFSZ, sizeof(char));
	for(np=1;np<PARSZ; np++) para[np] = para[np-1] + BUFSZ;

    vp_init();

	p=buf;
	while(1)
	{
		*p=fgetc(stdin);
		if(issep(*p, sep))
		{
			p = cmdline(buf, p, para);
		}else if(*p==EOF) 
		{
			p = cmdline(buf, p, para);
			break;
		}else	p++;
	}
	free(*para);
	free(para);
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

