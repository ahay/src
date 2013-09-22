/* Convert 3D Azimuth-IncidentAngle-Domain CIGs (x,y,azimuth,angle,tau) data to RSF format.

   Copyright (C) 2013 Tongji University, Shanghai, China 
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
        char *fn;

        sf_init(argc,argv);

        FILE *Fi;
        sf_file Fo;
        
        cjbsegy *tr;
 
        tr = calloc(sizeof(cjbsegy), 1);

        int nx, ny, nazim, nang, ntau;
        float dazim, dang, dtau, fazim=0.0, fang=0.0, ftau=0;
        int dx, dy, fx=0, fy=0;

        if (!sf_getint("nx",&nx)) nx=101;
        if (!sf_getint("ny",&ny)) ny=101;
        if (!sf_getint("nazim",&nazim)) nazim=8;
        if (!sf_getint("nang",&nang)) nang=21;
        if (!sf_getint("ntau",&ntau)) ntau=101;
        if (!sf_getint("dx",&dx)) dx=1;
        if (!sf_getint("dy",&dy)) dy=1;
        if (!sf_getfloat("dazim",&dazim)) dazim=22.5;
        if (!sf_getfloat("dang",&dang)) dang=2.0;
        if (!sf_getfloat("dtau",&dtau)) dtau=0.002;
        if (!sf_getint("fx",&fx)) fx=1;
        if (!sf_getint("fy",&fy)) fy=1;
        if (!sf_getfloat("ftau",&ftau)) ftau=0;
        if (NULL==(fn=sf_getstring("fn"))) fn="kpstm.ladcig.su.agc";

        /* setup I/O files */
        Fo = sf_output("out");

        if((Fi=fopen(fn,"rb"))==NULL)
        {
           printf("File %s open error!\n",fn);
           exit(0);
        }

        fread(tr,sizeof(cjbsegy),1,Fi);
        int iline0=tr->ep;
        sf_warning("ns=%d dt=%f iLineNo=%d ",tr->ns, tr->dt,iline0);


        if(fseek(Fi, 0L, 2) ==-1)
          printf("input file size unknown; Please specify n2\n");
        int nxy=(int) (ftell(Fi)/((60+ntau)*sizeof(float)));
       
        sf_warning("nxy=%d nx=%d ny=%d ",nxy, nx,ny);
        sf_warning("nazim=%d nang=%d ntau=%d",nazim,nang,ntau);
        sf_warning("dx=%d dy=%d dazim=%f dang=%f dtau=%f",dx,dy,dazim,dang,dtau);
        sf_warning("fx=%d fy=%d fazim=%f fang=%f ftau=%f",fx,fy,fazim,fang,ftau);

        if(nxy!=nx*ny*nazim*nang) {
          sf_warning("nx * ny * nazim * nang != nxy ");
          exit(0);  
         };

        sf_putint(Fo,"n1",ntau);
        sf_putint(Fo,"n2",nang);
        sf_putint(Fo,"n3",nazim);
        sf_putint(Fo,"n4",nx);
        sf_putint(Fo,"n5",ny);
        sf_putfloat(Fo,"d1",dtau);
        sf_putfloat(Fo,"o1",ftau);
        sf_putfloat(Fo,"d2",dang);
        sf_putfloat(Fo,"o2",fang);
        sf_putfloat(Fo,"d3",dazim);
        sf_putfloat(Fo,"o3",fazim);
        sf_putfloat(Fo,"d4",dx);
        sf_putfloat(Fo,"o4",fx);
        sf_putfloat(Fo,"d5",dy);
        sf_putfloat(Fo,"o5",fy);
        sf_putstring(Fo,"label1","z");
        sf_putstring(Fo,"label2","angle");
        sf_putstring(Fo,"label3","azimuth");
        sf_putstring(Fo,"label4","x");
        sf_putstring(Fo,"label5","y");
        sf_putstring(Fo,"unit1","ms");
        sf_putstring(Fo,"unit2","degree");
        sf_putstring(Fo,"unit3","degree");
        sf_putstring(Fo,"unit4","m");
        sf_putstring(Fo,"unit5","m");

        data = sf_floatalloc(ntau);

        rewind(Fi);
        for(i=0;;i++)
        {
          fread(tr,sizeof(cjbsegy),1,Fi);
          if(tr->ep != iline0){
            sf_warning("Read iLineNo=%d finished",iline0);
            iline0=tr->ep;
          }
          fread(data,sizeof(float),ntau,Fi);
          if(feof(Fi))break;

          sf_floatwrite(data, ntau, Fo);
        }
        sf_warning("Read iLineNo=%d finished",tr->ep);

        fclose(Fi);
        free(data);
        free(tr);
        exit(0);
}
