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
    unsigned char *data, head, ch;
    short int sid, tsid;
    short int num, tnum; 
    int i, j, k, ic;
    int TSL_HEADER_BYTE=16, dbit, length, channel=2, thead;
    bool format;
    float temp;
    char stu[TSL_DATA_BYTE];
    char *filename;
    FILE *file;    
    sf_file out, tfile;
    
    sf_init (argc, argv); 
    out = sf_output("out");

    if (NULL == (filename = sf_getstring("data"))) {
	/* input data */ 
	file = stdin;
    } else if (NULL == (file = fopen(filename,"rb"))) {
	sf_error("Cannot open \"%s\" for reading:",filename);
    }
    
    if (!sf_getbool("format",&format)) format=false;
    /* data format: [false] (TSL,TSH: 16) or [true] (TSn: 32) */
    
    if (!format) {
	TSL_HEADER_BYTE=16;
    } else {
	TSL_HEADER_BYTE=32;
    }
    
    if (NULL != sf_getstring("tfile")) {
	tfile = sf_output("tfile");
	sf_settype(tfile,SF_INT);	
    } else {
	tfile = NULL;
    }
    
    data=(unsigned char*)malloc(TSL_HEADER_BYTE);
    /* allocate data space */

    if(data==NULL) {
	sf_warning("Data space is null");
	exit (0);
    }
    
    for(i=0;i<TSL_HEADER_BYTE;i++) {
	
	if(i==0) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&data[0],1,1,file))
			abort();
	}
	if(i==1) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&data[1],1,1,file))
			abort();
	}
	if(i==2) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&data[2],1,1,file))
			abort();
	}
	if(i==3) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&data[3],1,1,file))
			abort();
	}
	if(i==4) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&data[4],1,1,file))
			abort();
	}
	if(i==5) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&data[5],1,1,file))
			abort();
	}
	if(i==6) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&data[6],1,1,file))
			abort();
		/* week */
	}
	if(i==7) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&data[7],1,1,file))
			abort();
	}
	if(i==8) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&sid,2,1,file))
			abort();
	}
	if(i==10) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&num,2,1,file))
			abort();
	}
	if(i==12) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&ch,1,1,file))
			abort();
		channel = (int)ch;
	}
	if(i==13) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&data[11],1,1,file))
			abort();
		/* tag: [32] (TSn) or [0] (TSH,TSL) */
		if (!format) {
		if(data[11]!=0) {
		    sf_error("wrong input parameter [format=false]!");
		}
	    } else {
		if(data[11]!=32) {
		    sf_error("wrong input parameter [format=true]!");
		}
	    }
	}
	if(i==14) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&data[12],1,1,file))
			abort();
		/* Status code: 0-8 */
	}
	if(i==15) {
	    fseek(file,i,SEEK_SET);
	    if (!fread(&data[13],1,1,file))
			abort();
		/* Bit pattern sign */
	}	
    }

    dbit = TSL_DATA_BYTE*channel;
    length = TSL_HEADER_BYTE+channel*TSL_DATA_BYTE*num;

    sf_warning("Record starts at %d0%d. %d. %d, %d: %d: %d",data[7],
	       data[5],data[4],data[3],data[2],data[1],data[0]);

    fseek(file,0L,SEEK_END);
    /*  seek end of file  */
    flen=ftell(file);
    /*  bit size of file  */
    flen=flen/(length);
    /*  record number     */
    
    sf_putint(out,"n1",flen*num);
    sf_putfloat(out,"o1",0.);
    sf_putfloat(out,"d1",1./num);
    sf_putstring(out,"label1","Time");
    sf_putstring(out,"unit1","s");
    sf_putint(out,"n2",channel);
    sf_putfloat(out,"o2",0.);
    sf_putfloat(out,"d2",1.);
    sf_putstring(out,"label2","Channel");
    sf_putstring(out,"unit2","Sample");

    if (NULL != tfile) {

	sf_putint(tfile,"n2",flen*num);
	sf_putfloat(tfile,"o2",0.);
	sf_putfloat(tfile,"d2",1./num);
	sf_putstring(tfile,"label2","Time");
	sf_putstring(tfile,"unit2","s");
	if (!format) {
	    sf_putint(tfile,"n1",TSL_HEADER_BYTE-2);
	} else {
	    sf_putint(tfile,"n1",TSL_HEADER_BYTE-6);
	}
	sf_putint(tfile,"o1",0);
	sf_putint(tfile,"d1",1);
	sf_putstring(tfile,"label1","Tag");
	sf_putstring(tfile,"unit1","Sample");
	
	
	/* write header file */
	for(i=0;i<flen;i++) {
	    sf_warning("header slice %d of %d;",i+1,flen);	    
	    for(j=0;j<num;j++) {
		for(k=0;k<TSL_HEADER_BYTE;k++) {
		    if(k==0) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&head,1,1,file))
				abort();
			thead = (int)head;
			sf_intwrite(&thead,1,tfile);
		    }
		    if(k==1) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&head,1,1,file))
				abort();
			thead = (int)head;
			sf_intwrite(&thead,1,tfile);	
		    }
		    if(k==2) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&head,1,1,file))
				abort();
			thead = (int)head;
			sf_intwrite(&thead,1,tfile);	
		    }
		    if(k==3) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&head,1,1,file))
				abort();
			thead = (int)head;
			sf_intwrite(&thead,1,tfile);	
		    }
		    if(k==4) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&head,1,1,file))
				abort();
			thead = (int)head;
			sf_intwrite(&thead,1,tfile);	
		    }
		    if(k==5) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&head,1,1,file))
				abort();
			thead = (int)head;
			sf_intwrite(&thead,1,tfile);	
		    }
		    if(k==6) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&head,1,1,file))
				abort();
			thead = (int)head;
			sf_intwrite(&thead,1,tfile);	
		    }
		    if(k==7) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&head,1,1,file))
				abort();
			thead = (int)head;
			sf_intwrite(&thead,1,tfile);	
		    }
		    if(k==8) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&tsid,2,1,file))
				abort();
			thead = (int)tsid;
			sf_intwrite(&thead,1,tfile);	
		    }
		    if(k==10) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&tnum,2,1,file))
				abort();
			thead = (int)tnum;
			sf_intwrite(&thead,1,tfile);	
		    }
		    if(k==12) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&head,1,1,file))
				abort();
			thead = (int)head;
			sf_intwrite(&thead,1,tfile);	
		    }
		    if(k==13) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&head,1,1,file))
				abort();
			thead = (int)head;
			sf_intwrite(&thead,1,tfile);	
		    }
		    if(k==14) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&head,1,1,file))
				abort();
			thead = (int)head;
			sf_intwrite(&thead,1,tfile);	
		    }
		    if(k==15) {
			fseek(file,i*length+k,SEEK_SET);
			if (!fread(&head,1,1,file))
				abort();
			thead = (int)head;
			sf_intwrite(&thead,1,tfile);	
		    }

		    if(format) {
			if(k==16) {
			    fseek(file,i*length+k,SEEK_SET);
			    if (!fread(&head,1,1,file))
					abort();
				thead = (int)head;
				sf_intwrite(&thead,1,tfile);	
			}
			if(k==17) {
			    fseek(file,i*length+k,SEEK_SET);
			    if (!fread(&head,1,1,file))
					abort();
				thead = (int)head;
				sf_intwrite(&thead,1,tfile);	
			}
			if(k==18) {
			    fseek(file,i*length+k,SEEK_SET);
			    if (!fread(&tnum,2,1,file))
					abort();
				thead = (int)tnum;
				sf_intwrite(&thead,1,tfile);	
			}
			if(k==20) {
			    fseek(file,i*length+k,SEEK_SET);
			    if (!fread(&head,1,1,file))
					abort();
				thead = (int)head;
				sf_intwrite(&thead,1,tfile);	
			}
			if(k==21) {
			    fseek(file,i*length+k,SEEK_SET);
			    if (!fread(&head,1,1,file))
					abort();
				thead = (int)head;
				sf_intwrite(&thead,1,tfile);	
			}
			if(k==22) {
			    fseek(file,i*length+k,SEEK_SET);
			    if (!fread(&thead,4,1,file))
					abort();
				sf_intwrite(&thead, 1, tfile);
			}
			if(k==26) {
			    fseek(file,i*length+k,SEEK_SET);
			    if (!fread(&head,1,1,file))
					abort();
				thead = (int)head;
				sf_intwrite(&thead,1,tfile);	
			}
			if(k==27) {
			    fseek(file,i*length+k,SEEK_SET);
			    if (!fread(&head,1,1,file))
					abort();
				thead = (int)head;
				sf_intwrite(&thead,1,tfile);	
			}			
			if(k==28) {
			    fseek(file,i*length+k,SEEK_SET);
			    if (!fread(&head,1,1,file))
					abort();
				thead = (int)head;
				sf_intwrite(&thead,1,tfile);	
			}
			if(k==29) {
			    fseek(file,i*length+k,SEEK_SET);
			    if (!fread(&head,1,1,file))
					abort();
				thead = (int)head;
				sf_intwrite(&thead,1,tfile);	
			}
			if(k==30) {
			    fseek(file,i*length+k,SEEK_SET);
			    if (!fread(&head,1,1,file))
					abort();
				thead = (int)head;
				sf_intwrite(&thead,1,tfile);	
			}
			if(k==31) {
			    fseek(file,i*length+k,SEEK_SET);
			    if (!fread(&head,1,1,file))
					abort();
				thead = (int)head;
				sf_intwrite(&thead,1,tfile);	
			}
		    } 
		}
	    }   
	}
    }

    /* write data file */
    for(ic=0;ic<channel;ic++) {
	for(i=0;i<flen;i++) {
	    for(j=0;j<num;j++) {
		for(k=0;k<TSL_DATA_BYTE;k++) {
		    fseek(file,(TSL_HEADER_BYTE+j*dbit+ic*TSL_DATA_BYTE+k)+
			  length*i,SEEK_SET); 
		    if (!fread(&stu[k],1,1,file))
				abort();
			s.byte[k] = stu[k];
		}
		temp=(float) s.t.data;
		sf_floatwrite(&temp,1,out);	    
	    }
	}
    }

    exit (0);
}
