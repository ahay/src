/* endianness conversion */
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

size_t size;

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

void swap4_int(int *p)
{
	unsigned int n = *(unsigned int *)p;
	swap4(&n);
	*p = *(int *)&n;
}

void swap4_float(float *p)
{
	unsigned int n = *(unsigned int *)p;
	swap4(&n);
	*p = *(float *)&n;
}

void count()
{
	static int b=0, mb=0;

	if ((b += size) >= MB) {
		b = 0;
		sf_warning("%d MB written;",++mb);
	}
}

int main(int argc, char *argv[])
{
    int i,n,n1;
	off_t s,s1;
	void *p;
	bool verb;
	sf_file Fin,Fout;
    sf_datatype type;

	sf_init(argc,argv);

	Fin = sf_input("in");
	Fout= sf_output("out");

    n1 = sf_filesize(Fin);

    if (!sf_getlargeint("n",(off_t *)&n)) n = n1;
    /* num of elements to convert */
	if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity */

    if (n < n1) {
        sf_putlargeint(Fout,"n1",(off_t)n);

        if (sf_leftsize(Fin,1) > 1) {
		sf_warning("squeez");
			if (sf_histint(Fin,"n2",&i)) sf_putint(Fout,"n2",1);
			if (sf_histint(Fin,"n3",&i)) sf_putint(Fout,"n3",1);
			if (sf_histint(Fin,"n4",&i)) sf_putint(Fout,"n4",1);
		}
    }

    switch (type = sf_gettype(Fin)) {
    case SF_INT   : size = sizeof(int  ); break;
    case SF_FLOAT : size = sizeof(float); break;
    default : sf_error("need type={int,float} in input");
    }

    s = n * size;
    s1= n1* size;
    p = (void *)sf_alloc(1,size);

	sf_seek(Fin,0,SEEK_SET);

    while(sf_tell(Fin) < s) {
        switch (sf_gettype(Fin)) {
        case SF_INT :
			sf_intread(p,1,Fin);
			swap4_int((int *)p);
			sf_intwrite(p,1,Fout);
            break;
        case SF_FLOAT :
			sf_floatread(p,1,Fin);
			swap4_float((float *)p);
			sf_floatwrite(p,1,Fout);
            break;
        default  :
            sf_error("need type={int,float} in input");
		}
        if (verb) count();
	}
	if (verb) sf_warning("%\nd of %d elements converted\n",sf_tell(Fin) / size,n1);

	sf_fileclose(Fin);
	sf_fileclose(Fout);
	return 0;
}
