#include <rsf.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[])
{
	char *path, *name1, *name2, number[5];
	int i;

	sf_file in, out, tmp;
	sf_init(argc, argv);

	out=sf_output("out");
	in=sf_input("in");
	path=sf_getstring("path");
	if(path==NULL) sf_error("Need path!");
	sf_warning("path is %s",path);

	name1=sf_charalloc(strlen(path));
	name2=sf_charalloc(strlen(path));
	strcpy(name1,path);
	strcpy(name2,path);
	i=10;
	sprintf(number,"%d",i);
	strcat(name1,number);
	strcat(name1,"f.rsf");
	strcat(name2,number);
	strcat(name2,"b");
	sf_warning("name1 is %s",name1);
	sf_warning("name2 is %s", name2);
	tmp=sf_input(name1);

	return 0;
}
