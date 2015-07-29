/* Convert 2D cjb-SU data to RSF format.

   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng
     
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

/* prepared head files by myself */
#include "_cjb.h"
#include "_cjbsegy.h"

/* head files aumatically produced from *.c */
#include "zero.h"
#include "puthead.h"

/*****************************************************************************************/
int main(int argc, char* argv[])
{
	int   i;
        float *data;
        cjbsegy  *tr;
        char *fn;
        FILE *Fi;
        sf_file Fo;

        int   nz, nx;
        float fz;
        float fx;
        float dz;
        float dx;

        sf_init(argc,argv);

 
        if (NULL==(fn=sf_getstring("fn"))) fn="hess.rtm.absorb.vti.image.su.pc_LF";

        /* setup I/O files */
        Fo = sf_output("out"); /* pseudo-pure P-wave x-component */

        if((Fi=fopen(fn,"rb"))==NULL)
        {
           printf("File %s open error!\n",fn);
           exit(0);
        }

        tr = calloc(sizeof(cjbsegy),1);
        if (1 != fread(tr,sizeof(cjbsegy),1,Fi)) 
	    sf_error("fread error:");

        nz=(int)tr->ns;
        fz=tr->f1;
        fx=-5000.0*6.096*2/10;
        dz=6.096*2;
        dx=6.096*2;

        if(fseek(Fi, 0L, 2) ==-1)
          printf("input file size unknown; Please specify n2\n");
        nx=(int) (ftell(Fi)/((nz+60)*sizeof(float)));

        sf_warning("nx=%d nz=%d dx=%f dz=%f",nx,nz,dx,dz);

        sf_putint(Fo,"n1",nz);
        sf_putint(Fo,"n2",nx);
        sf_putfloat(Fo,"d1",dz);
        sf_putfloat(Fo,"o1",fz);
        sf_putfloat(Fo,"d2",dx);
        sf_putfloat(Fo,"o2",fx);
        sf_putstring(Fo,"label1","z");
        sf_putstring(Fo,"label2","x");
        sf_putstring(Fo,"unit1","m");
        sf_putstring(Fo,"unit2","m");

        data = sf_floatalloc(nz);

        rewind(Fi);
        for(i=0;;i++)
        {
	    if (1 != fread(tr,sizeof(cjbsegy),1,Fi))
		sf_error("fread error:");
          if(feof(Fi))break;
          if (nz != fread(data,sizeof(float),nz,Fi))
	      sf_error("fread error:");

          sf_floatwrite(data, nz, Fo);
        }

        fclose(Fi);
        free(tr);
        free(data);
        exit(0);
}
