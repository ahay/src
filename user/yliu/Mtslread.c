/* Convert a TSL (MT, V5-2000 of Phoenix Geophysics) dataset to RSF. */
/*
  Copyright (C) 2015 Jilin University
  
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

#define TSL_DATA_BYTE 3
/*^*/

union {
    char byte[TSL_DATA_BYTE];
    struct {
	int data:(TSL_DATA_BYTE*8);
    } t;
} s;
/*^*/

int main (int argc, char* argv[])
{
    long flen;
    unsigned char *data;  /* equipment info   */
    short int sid;        /* equipment number */
    short int num;        /* scan number      */
    int i, j, format, channel, TSL_HEADER_BYTE=16;
    float temp;
    char stu1, stu2, stu3;
    char *filename;
    FILE *file;    
    sf_file out;
    
    sf_init (argc, argv); 
    out = sf_output("out");

    if (NULL == (filename = sf_getstring("data"))) {
	/* input data */ 
	file = stdin;
    } else if (NULL == (file = fopen(filename,"rb"))) {
	sf_error("Cannot open \"%s\" for reading:",filename);
    }
    
    if (!sf_getint("format",&format)) format=0;
    /* data format: [0] (TSH,TSL) or [1] (TSn) */
    
    if (0==format) {
	TSL_HEADER_BYTE=16;
    } else if (1==format) {
	TSL_HEADER_BYTE=32;
    } else {
	sf_error("wrong format=%d",format);
    }

    data=(unsigned char*)malloc(TSL_HEADER_BYTE);
    /* allocate data space */
    
    for(i=0;i<TSL_HEADER_BYTE;i++) {
	
	if(i==0) {
	    fseek(file,i,SEEK_SET);
	    fread(&data[0],1,1,file);
	}
	if(i==1) {
	    fseek(file,i,SEEK_SET);
	    fread(&data[1],1,1,file);
	}
	if(i==2) {
	    fseek(file,i,SEEK_SET);
	    fread(&data[2],1,1,file);
	}
	if(i==3) {
	    fseek(file,i,SEEK_SET);
	    fread(&data[3],1,1,file);
	}
	if(i==4) {
	    fseek(file,i,SEEK_SET);
	    fread(&data[4],1,1,file);
	}
	if(i==5) {
	    fseek(file,i,SEEK_SET);
	    fread(&data[5],1,1,file);
	}
	if(i==6) {
	    fseek(file,i,SEEK_SET);
	    fread(&data[6],1,1,file);
	    /* week */
	}
	if(i==7) {
	    fseek(file,i,SEEK_SET);
	    fread(&data[7],1,1,file);
	}
	if(i==8) {
	    fseek(file,i,SEEK_SET);
	    fread(&sid,2,1,file);
	}
	if(i==10) {
	    fseek(file,i,SEEK_SET);
	    fread(&num,2,1,file);
	}
	if(i==12) {
	    fseek(file,i,SEEK_SET);
	    fread(&channel,1,1,file);
	}
	if(i==13) {
	    fseek(file,i,SEEK_SET);
	    fread(&data[11],1,1,file);
	    /* Data format: [0] (TSn) or [1] (TSH,TSL) */
	    if (format != data[11]) {
		sf_error("wrong input parameter [format=%d]!",format);
	    }
	}
	if(i==14) {
	    fseek(file,i,SEEK_SET);
	    fread(&data[12],1,1,file);
	    /* Status code: 0-8 */
	}
	if(i==15) {
	    fseek(file,i,SEEK_SET);
	    fread(&data[13],1,1,file);
	    /* Bit pattern sign */
	}	
    }
    sf_warning("Record starts at %d0%d. %d. %d, %d: %d: %d",data[7],
	       data[5],data[4],data[3],data[2],data[1],data[0]);

    fseek(file,0L,SEEK_END);             /*  seek end of file */
    flen=ftell(file);                    /*  file bit size    */
    flen=flen/(TSL_HEADER_BYTE+channel*TSL_DATA_BYTE*num);
    /*  record number   */
    
    if(data==NULL) {
	sf_warning("Data space is null");
	exit (0);
    }
 
    sf_putint(out,"n1",flen*num);
    sf_putfloat(out,"o1",0.);
    sf_putfloat(out,"d1",1./24.);
    sf_putstring(out,"label1","Time");
    sf_putstring(out,"unit1","s");
    sf_putint(out,"n2",channel);
    sf_putfloat(out,"o2",0.);
    sf_putfloat(out,"d2",1.);
    sf_putstring(out,"label2","Channel");
    sf_putstring(out,"unit2","Sample");
        
    for(j=0;j<flen;j++) {
        for(i=0;i<num;i++) {
	    fseek(file,(16+i*6)+160*j,SEEK_SET); 
	    fread(&stu1,1,1,file);
	    s.byte[0]=stu1;
	    fseek(file,(16+(i*6+1))+160*j,SEEK_SET);
	    fread(&stu2,1,1,file);
	    s.byte[1]=stu2;
	    fseek(file,(16+(i*6+2))+160*j,SEEK_SET);
	    fread(&stu3,1,1,file);
	    s.byte[2]=stu3;
	    temp=(float) s.t.data;

	    sf_floatwrite(&temp,1,out);	    
	}
    }
    for(j=0;j<flen;j++) {
        for(i=0;i<num;i++) {
	    fseek(file,(16+i*6+3)+160*j,SEEK_SET); 
	    fread(&stu1,1,1,file);
	    s.byte[0]=stu1;
	    fseek(file,(16+(i*6+4))+160*j,SEEK_SET);
	    fread(&stu2,1,1,file);
	    s.byte[1]=stu2;
	    fseek(file,(16+(i*6+5))+160*j,SEEK_SET);
	    fread(&stu3,1,1,file);
	    s.byte[2]=stu3;
	    temp=(float) s.t.data;

	    sf_floatwrite(&temp,1,out);	    
	}
    }    
    exit (0);
}
