/* byte swapping (LP64) */
/*
  Copyright (C) 2011 KAUST

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

#define MB 1048576

bool verb;
int size;
static int b,mb;

int is_little_endian(void)
{
	int n = 0x1;
	char *c = (char *)&n;
	return (c[0] == 1);
}

void swap2(unsigned int *p)
{
	*p = ((*p >> 8) & 0x00ff) | ((*p << 8) & 0xff00);
}

void swap4(unsigned int *p)
{
	*p = ((*p >> 8 ) & 0x00ff00ff) | ((*p << 8 ) & 0xff00ff00);
	*p = ((*p >> 16) & 0x0000ffff) | ((*p << 16) & 0xffff0000);
}

void swap8(unsigned long *p)
{
	*p = ((*p >> 8 ) & 0x00ff00ff00ff00ff) | ((*p << 8 ) & 0xff00ff00ff00ff00);
	*p = ((*p >> 16) & 0x0000ffff0000ffff) | ((*p << 16) & 0xffff0000ffff0000);
	*p = ((*p >> 32) & 0x00000000ffffffff) | ((*p << 32) & 0xffffffff00000000);
}

void swap_int_4(int *p)
{
	unsigned int n = *(unsigned int *)p;
	swap4(&n);
	*p = *(int *)&n;
}

void swap_float_4(float *p)
{
	unsigned int n = *(unsigned int *)p;
	swap4(&n);
	*p = *(float *)&n;
}

void count()
{
	if ((b += size) >= MB) {
		b = 0;
		if (!(++mb % 100))
			sf_warning("%d MB written;",mb);
	}
}

int main(int argc, char *argv[])
{
	char *type;
	void *p;
	off_t filesize;
	sf_file Fin,Fout;

	sf_init(argc,argv);

	Fin = sf_input("in");
	Fout= sf_output("out");

	if(!sf_getbool("verb",&verb)) verb = false; /* verbosity */
	type = sf_getstring("type"); /* int, float */

	sf_seek(Fin,0,SEEK_END);
	filesize = sf_tell(Fin);
	sf_seek(Fin,0,SEEK_SET);

	b = mb = 0;

	switch (type[0]) {
	case 'i' :
		p = (void *)sf_alloc(1,size = sizeof(int));
		while(sf_tell(Fin) < filesize) {
			sf_intread(p,1,Fin);
			swap_int_4((int *)p);
			sf_intwrite(p,1,Fout);
			if (verb) count();
		}
		break;
	case 'f' :
	default  :
		p = (void *)sf_alloc(1,size = sizeof(float));
		while(sf_tell(Fin) < filesize) {
			sf_floatread(p,1,Fin);
			swap_float_4((float *)p);
			sf_floatwrite(p,1,Fout);
			if (verb) count();
		}
	}
	if (verb) sf_warning("\n");

	sf_fileclose(Fin);
	sf_fileclose(Fout);
	return 0;
}
